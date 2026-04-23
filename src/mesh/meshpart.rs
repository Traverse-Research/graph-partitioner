use crate::types::{Idx, Real, Result};
use crate::Mesh;
use crate::ctrl::Control;
use crate::graph::setup::setup_graph;
use crate::mesh::{create_graph_dual, build_node_element_csr};
use crate::partition;

/// Entry point for METIS_PartMeshDual logic.
pub fn part_dual(mesh: Mesh, epart: &mut [Idx], npart: &mut [Idx]) -> Result<Idx> {
    let ne = (mesh.element_offsets.len() - 1) as Idx;
    let nn = mesh.nn;
    let nparts = mesh.num_parts;
    let min_common_nodes = mesh.min_common_nodes;

    // Build dual graph
    let (xadj, adjacency) = create_graph_dual(ne, nn, mesh.element_offsets, mesh.element_indices, min_common_nodes);

    // Set up control (is_kway = true, matching C METIS which calls METIS_PartGraphKway)
    let mut ctrl = Control::new(&mesh.options, 1, nparts, true);
    if let Some(tp) = mesh.target_part_weights {
        ctrl.set_target_part_weights(tp);
    }

    // Create internal graph from dual
    let mut gdata = setup_graph(
        1,
        &xadj,
        &adjacency,
        mesh.vertex_weights,
        mesh.vertex_sizes,
        None,
    );

    // SetupKWayBalMultipliers (matches part_kway flow)
    let inv_total_vertex_weight = gdata.inv_total_vertex_weight.clone();
    ctrl.setup_kway_balance_multipliers(&inv_total_vertex_weight);

    // Compute CoarsenTo and nIparts (matches part_kway / METIS_PartGraphKway)
    let log2_nparts = (nparts as f64).log2().max(1.0) as Idx;
    let coarsen_to = ((gdata.num_vertices / (40 * log2_nparts)).max(30 * nparts)).max(1);
    ctrl.coarsen_to = coarsen_to;
    ctrl.num_init_parts = if coarsen_to == 30 * nparts { 4 } else { 5 };

    // Partition the dual graph
    let cut = partition::mlevel_kway_partitioning(&mut ctrl, &mut gdata, epart);

    // Build node-element CSR (matching C METIS METIS_PartMeshDual lines 136-150)
    let (nptr, nind) = build_node_element_csr(ne, nn, mesh.element_offsets, mesh.element_indices);

    // Derive node partition from element partition
    // (matching C METIS InduceRowPartFromColumnPart exactly)
    induce_row_part_from_column_part(
        nn as usize,
        &nptr,
        &nind,
        npart,
        epart,
        nparts as usize,
        None, // target_part_weights — NULL in C METIS when not provided
    );

    Ok(cut)
}

/// Derive node partition from element partition.
/// Matches C METIS InduceRowPartFromColumnPart from meshpart.c exactly.
///
/// Two-phase algorithm:
/// Phase 1: Assign nodes where ALL connected elements are in the same partition.
/// Phase 2: Assign remaining nodes (multi-partition) using iargmax with load balancing.
fn induce_row_part_from_column_part(
    nrows: usize,
    rowptr: &[Idx],
    rowind: &[Idx],
    rpart: &mut [Idx],
    cpart: &[Idx],
    nparts: usize,
    target_part_weights: Option<&[Real]>,
) {
    let mut part_weights = vec![0 as Idx; nparts];
    let mut neighbor_partitions = vec![0 as Idx; nparts];
    let mut neighbor_weights = vec![0 as Idx; nparts];
    let mut neighbor_mark = vec![-1 as Idx; nparts];

    // Initialize rpart to -1 (unassigned)
    for i in 0..nrows {
        rpart[i] = -1;
    }

    // Setup integer target partition weights
    let itarget_part_weights: Vec<Idx> = if let Some(tp) = target_part_weights {
        (0..nparts).map(|i| 1 + (nrows as Real * tp[i]) as Idx).collect()
    } else {
        vec![1 + (nrows / nparts) as Idx; nparts]
    };

    // Phase 1: Assign nodes where all connected elements are in the same partition.
    // Empty nodes (no elements) get -2 (permanently unassigned).
    for i in 0..nrows {
        let start = rowptr[i] as usize;
        let end = rowptr[i + 1] as usize;

        if start == end {
            rpart[i] = -2;
            continue;
        }

        let me = cpart[rowind[start] as usize];
        let mut all_same = true;
        for j in (start + 1)..end {
            if cpart[rowind[j] as usize] != me {
                all_same = false;
                break;
            }
        }

        if all_same {
            rpart[i] = me;
            part_weights[me as usize] += 1;
        }
    }

    // Phase 2: Assign remaining nodes (those touching multiple partitions)
    // in a balanced way.
    for i in 0..nrows {
        if rpart[i] != -1 {
            continue;
        }

        let start = rowptr[i] as usize;
        let end = rowptr[i + 1] as usize;

        // Count neighbor partitions and their element counts
        let mut nnbrs: usize = 0;
        for j in start..end {
            let me = cpart[rowind[j] as usize] as usize;
            if neighbor_mark[me] == -1 {
                neighbor_partitions[nnbrs] = me as Idx;
                neighbor_weights[nnbrs] = 1;
                neighbor_mark[me] = nnbrs as Idx;
                nnbrs += 1;
            } else {
                neighbor_weights[neighbor_mark[me] as usize] += 1;
            }
        }

        // Assign to partition with most connected elements (iargmax)
        let best_idx = iargmax(&neighbor_weights[..nnbrs]);
        rpart[i] = neighbor_partitions[best_idx];

        // If overweight, reassign to the lightest partition among neighbors
        let pi = rpart[i] as usize;
        if part_weights[pi] > itarget_part_weights[pi] {
            for j in 0..nnbrs {
                let nd = neighbor_partitions[j] as usize;
                if part_weights[nd] < itarget_part_weights[nd]
                    || (part_weights[nd] - itarget_part_weights[nd]) < (part_weights[pi] - itarget_part_weights[pi])
                {
                    rpart[i] = neighbor_partitions[j];
                    break;
                }
            }
        }
        part_weights[rpart[i] as usize] += 1;

        // Reset neighbor_mark array
        for j in 0..nnbrs {
            neighbor_mark[neighbor_partitions[j] as usize] = -1;
        }
    }
}

/// Return the index of the maximum element (first occurrence for ties).
/// Matches GKlib iargmax with incx=1.
fn iargmax(x: &[Idx]) -> usize {
    let mut max_i = 0;
    for i in 1..x.len() {
        if x[i] > x[max_i] {
            max_i = i;
        }
    }
    max_i
}
