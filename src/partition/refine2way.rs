use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::partition::fm::fm_2way_cut_refine;
use crate::partition::balance::balance_2way;
use crate::partition::initpart::compute_2way_partition_params;

/// Refine2Way: refine the 2-way partition through the entire coarsening chain.
///
/// Matches the C METIS Refine2Way algorithm:
///   - `orggraph` is a raw pointer to the original (finest) graph
///   - `graph` is the coarsest graph where initial partitioning was performed
///   - Starts at the coarsest level, iterates: Balance -> FM refine -> project to finer
///   - Stops when we reach the original graph
pub fn refine_2way(ctrl: &mut Control, orggraph: *mut GraphData, graph: &mut GraphData, target_part_weights: &[Real]) {
    let niter = ctrl.num_iter;

    // Compute 2-way partition parameters at the coarsest level
    compute_2way_partition_params(ctrl, graph);

    // Use raw pointer to walk the chain from coarsest to finest
    let mut cur: *mut GraphData = graph as *mut GraphData;

    loop {
        let g = unsafe { &mut *cur };

        // Balance and refine at the current level
        balance_2way(ctrl, g, target_part_weights);
        fm_2way_cut_refine(ctrl, g, target_part_weights, niter);

        // If we have reached the original (finest) graph, we are done
        if cur == orggraph {
            break;
        }

        // Move to finer level
        let finer_ptr = g.finer;
        assert!(!finer_ptr.is_null(), "refine_2way: finer pointer is null before reaching orggraph");

        // Project the partition from the current (coarser) level to the finer level
        project_2way_partition(unsafe { &mut *finer_ptr });

        cur = finer_ptr;
    }
}

/// Project2WayPartition: project the 2-way partition from the coarser graph to the
/// finer graph.
///
/// Matches the C METIS Project2WayPartition algorithm:
///   - `graph` is the finer level; `graph.coarser` is the coarser level
///   - Uses `graph.coarse_map` to map fine vertices to coarse vertices
///   - Reuses `coarse_map[i]` to cache `cboundary_map[coarse_map[i]]` for the interior-node optimization
///   - Computes id/ed for each vertex and builds the boundary list
///   - Copies edge_cut and part_weights from the coarser graph
///   - Frees (drops) the coarser graph
fn project_2way_partition(graph: &mut GraphData) {
    // Extract all needed data from the coarser graph before mutating this graph.
    // This avoids borrow checker conflicts (coarser is owned by graph).
    let (cwhere, cboundary_map, cedge_cut, ncon, cpart_weights) = {
        let cgraph = graph.coarser.as_ref()
            .expect("project_2way_partition: coarser graph is None");
        let ncon = cgraph.num_constraints as usize;
        (
            cgraph.partition.clone(),
            cgraph.boundary_map.clone(),
            cgraph.edge_cut,
            ncon,
            cgraph.part_weights[..2 * ncon].to_vec(),
        )
    };

    // Allocate 2-way partition memory for this (finer) level
    graph.alloc_2way();

    let num_vertices = graph.num_vertices as usize;

    // Project partition: where[i] = cwhere[coarse_map[i]]
    // Also cache cboundary_map for the interior-node optimization: coarse_map[i] = cboundary_map[coarse_map[i]]
    for i in 0..num_vertices {
        let j = graph.coarse_map[i] as usize;
        graph.partition[i] = cwhere[j];
        graph.coarse_map[i] = cboundary_map[j]; // Reuse coarse_map to store boundary info
    }

    // Compute id/ed for each vertex and build the boundary list
    let mut num_boundary: Idx = 0;

    for i in 0..num_vertices {
        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;
        let mut tid: Idx = 0;
        let mut ted: Idx = 0;

        if graph.coarse_map[i] == -1 {
            // Interior node: the coarser node was interior (cboundary_map == -1),
            // so all neighbors must be in the same partition.
            for j in istart..iend {
                tid += graph.edge_weights[j];
            }
        } else {
            // Potentially interface node: check each neighbor's partition
            let me = graph.partition[i];
            for j in istart..iend {
                if me == graph.partition[graph.adjacency[j] as usize] {
                    tid += graph.edge_weights[j];
                } else {
                    ted += graph.edge_weights[j];
                }
            }
        }

        graph.internal_degree[i] = tid;
        graph.external_degree[i] = ted;

        if ted > 0 || istart == iend {
            // BNDInsert(num_boundary, boundary_list, boundary_map, i)
            graph.boundary_list[num_boundary as usize] = i as Idx;
            graph.boundary_map[i] = num_boundary;
            num_boundary += 1;
        }
    }

    graph.edge_cut = cedge_cut;
    graph.num_boundary = num_boundary;

    // Copy part_weights from coarser
    graph.part_weights[..2 * ncon].copy_from_slice(&cpart_weights);

    // Free coarser graph (FreeGraph(&graph->coarser))
    graph.coarser = None;
}
