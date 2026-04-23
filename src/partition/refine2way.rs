use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::partition::fm::fm_2way_cut_refine;
use crate::partition::balance::balance_2way;
use crate::partition::initpart::compute_2way_partition_params;

/// Refine2Way: refine the 2-way partition through the entire coarsening chain.
///
/// `graph` is the finest level. `levels` contains the coarser levels where
/// `levels[0]` is the first coarser level and `levels[last]` is the coarsest
/// (where initial partitioning was performed).
///
/// Walks from coarsest to finest: Balance -> FM refine -> project to finer.
pub fn refine_2way(ctrl: &mut Control, graph: &mut GraphData, levels: &mut Vec<GraphData>, target_part_weights: &[Real]) {
    let niter = ctrl.num_iter;
    let nlevels = levels.len(); // number of coarser levels

    // Compute 2-way partition parameters at the coarsest level
    compute_2way_partition_params(ctrl, &mut levels[nlevels - 1]);

    // Balance and refine the coarsest level
    balance_2way(ctrl, &mut levels[nlevels - 1], target_part_weights);
    fm_2way_cut_refine(ctrl, &mut levels[nlevels - 1], target_part_weights, niter);

    // Walk from coarsest to finest within the arena
    for i in (0..nlevels - 1).rev() {
        // Project partition from levels[i+1] (coarser) to levels[i] (finer)
        {
            let (left, right) = levels.split_at_mut(i + 1);
            let finer = &mut left[i];
            let coarser = &right[0]; // = levels[i+1]
            project_2way_partition(finer, coarser);
        }

        // Balance and refine at this level
        balance_2way(ctrl, &mut levels[i], target_part_weights);
        fm_2way_cut_refine(ctrl, &mut levels[i], target_part_weights, niter);
    }

    // Final projection: from levels[0] (first coarser level) to graph (finest)
    project_2way_partition(graph, &levels[0]);

    // Balance and refine at the finest level
    balance_2way(ctrl, graph, target_part_weights);
    fm_2way_cut_refine(ctrl, graph, target_part_weights, niter);
}

/// Project2WayPartition: project the 2-way partition from the coarser graph to the
/// finer graph.
///
/// `finer` is the finer level; `coarser` is the coarser level.
/// Uses `finer.coarse_map` to map fine vertices to coarse vertices.
/// Reuses `coarse_map[i]` to cache `cboundary_map[coarse_map[i]]` for the interior-node optimization.
/// Computes id/ed for each vertex and builds the boundary list.
/// Copies edge_cut and part_weights from the coarser graph.
fn project_2way_partition(finer: &mut GraphData, coarser: &GraphData) {
    let ncon = coarser.num_constraints as usize;
    let cedge_cut = coarser.edge_cut;
    let cpart_weights = &coarser.part_weights[..2 * ncon];

    // Allocate 2-way partition memory for this (finer) level
    finer.alloc_2way();

    let num_vertices = finer.num_vertices as usize;

    // Project partition: where[i] = cwhere[coarse_map[i]]
    // Also cache cboundary_map for the interior-node optimization: coarse_map[i] = cboundary_map[coarse_map[i]]
    for i in 0..num_vertices {
        let j = finer.coarse_map[i] as usize;
        finer.partition[i] = coarser.partition[j];
        finer.coarse_map[i] = coarser.boundary_map[j]; // Reuse coarse_map to store boundary info
    }

    // Compute id/ed for each vertex and build the boundary list
    let mut num_boundary: Idx = 0;

    for i in 0..num_vertices {
        let istart = finer.xadj[i] as usize;
        let iend = finer.xadj[i + 1] as usize;
        let mut tid: Idx = 0;
        let mut ted: Idx = 0;

        if finer.coarse_map[i] == -1 {
            // Interior node: the coarser node was interior (cboundary_map == -1),
            // so all neighbors must be in the same partition.
            for j in istart..iend {
                tid += finer.edge_weights[j];
            }
        } else {
            // Potentially interface node: check each neighbor's partition
            let me = finer.partition[i];
            for j in istart..iend {
                if me == finer.partition[finer.adjacency[j] as usize] {
                    tid += finer.edge_weights[j];
                } else {
                    ted += finer.edge_weights[j];
                }
            }
        }

        finer.internal_degree[i] = tid;
        finer.external_degree[i] = ted;

        if ted > 0 || istart == iend {
            // BNDInsert(num_boundary, boundary_list, boundary_map, i)
            finer.boundary_list[num_boundary as usize] = i as Idx;
            finer.boundary_map[i] = num_boundary;
            num_boundary += 1;
        }
    }

    finer.edge_cut = cedge_cut;
    finer.num_boundary = num_boundary;

    // Copy part_weights from coarser
    finer.part_weights[..2 * ncon].copy_from_slice(cpart_weights);
}
