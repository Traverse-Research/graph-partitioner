use crate::types::{Idx, Real, Result};
use crate::Mesh;
use crate::ctrl::Ctrl;
use crate::graph::setup::setup_graph;
use crate::mesh::{create_graph_dual, build_node_element_csr};
use crate::partition;

/// Entry point for METIS_PartMeshDual logic.
pub fn part_dual(mesh: Mesh, epart: &mut [Idx], npart: &mut [Idx]) -> Result<Idx> {
    let ne = (mesh.eptr.len() - 1) as Idx;
    let nn = mesh.nn;
    let nparts = mesh.nparts;
    let ncommon = mesh.ncommon;

    // Build dual graph
    let (xadj, adjncy) = create_graph_dual(ne, nn, mesh.eptr, mesh.eind, ncommon);

    // Set up control (is_kway = true, matching C METIS which calls METIS_PartGraphKway)
    let mut ctrl = Ctrl::new(&mesh.options, 1, nparts, true);
    if let Some(tp) = mesh.tpwgts {
        ctrl.set_tpwgts(tp);
    }

    // Create internal graph from dual
    let mut gdata = setup_graph(
        1,
        &xadj,
        &adjncy,
        mesh.vwgt,
        mesh.vsize,
        None,
    );

    // SetupKWayBalMultipliers (matches part_kway flow)
    let invtvwgt = gdata.invtvwgt.clone();
    ctrl.setup_kway_bal_multipliers(&invtvwgt);

    // Compute CoarsenTo and nIparts (matches part_kway / METIS_PartGraphKway)
    let log2_nparts = (nparts as f64).log2().max(1.0) as Idx;
    let coarsen_to = ((gdata.nvtxs / (40 * log2_nparts)).max(30 * nparts)).max(1);
    ctrl.coarsen_to = coarsen_to;
    ctrl.niparts = if coarsen_to == 30 * nparts { 4 } else { 5 };

    // Partition the dual graph
    let cut = partition::mlevel_kway_partitioning(&mut ctrl, &mut gdata, epart);

    // Build node-element CSR (matching C METIS METIS_PartMeshDual lines 136-150)
    let (nptr, nind) = build_node_element_csr(ne, nn, mesh.eptr, mesh.eind);

    // Derive node partition from element partition
    // (matching C METIS InduceRowPartFromColumnPart exactly)
    induce_row_part_from_column_part(
        nn as usize,
        &nptr,
        &nind,
        npart,
        epart,
        nparts as usize,
        None, // tpwgts — NULL in C METIS when not provided
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
    tpwgts: Option<&[Real]>,
) {
    let mut pwgts = vec![0 as Idx; nparts];
    let mut nbrdom = vec![0 as Idx; nparts];
    let mut nbrwgt = vec![0 as Idx; nparts];
    let mut nbrmrk = vec![-1 as Idx; nparts];

    // Initialize rpart to -1 (unassigned)
    for i in 0..nrows {
        rpart[i] = -1;
    }

    // Setup integer target partition weights
    let itpwgts: Vec<Idx> = if let Some(tp) = tpwgts {
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
            pwgts[me as usize] += 1;
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
            if nbrmrk[me] == -1 {
                nbrdom[nnbrs] = me as Idx;
                nbrwgt[nnbrs] = 1;
                nbrmrk[me] = nnbrs as Idx;
                nnbrs += 1;
            } else {
                nbrwgt[nbrmrk[me] as usize] += 1;
            }
        }

        // Assign to partition with most connected elements (iargmax)
        let best_idx = iargmax(&nbrwgt[..nnbrs]);
        rpart[i] = nbrdom[best_idx];

        // If overweight, reassign to the lightest partition among neighbors
        let pi = rpart[i] as usize;
        if pwgts[pi] > itpwgts[pi] {
            for j in 0..nnbrs {
                let nd = nbrdom[j] as usize;
                if pwgts[nd] < itpwgts[nd]
                    || (pwgts[nd] - itpwgts[nd]) < (pwgts[pi] - itpwgts[pi])
                {
                    rpart[i] = nbrdom[j];
                    break;
                }
            }
        }
        pwgts[rpart[i] as usize] += 1;

        // Reset nbrmrk array
        for j in 0..nnbrs {
            nbrmrk[nbrdom[j] as usize] = -1;
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
