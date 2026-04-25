pub mod balance;
pub mod fm;
pub mod initpart;
pub mod kway;
pub mod kwayfm;
pub mod kwayrefine;
pub mod pmetis;
pub mod recursive;
pub mod refine2way;
pub mod split;

use crate::ctrl::Control;
use crate::graph::coarsen::coarsen_graph;
use crate::graph::GraphData;
use crate::types::{Idx, Real};

/// MlevelKWayPartitioning: multilevel k-way partitioning orchestration.
///
/// Matches C METIS MlevelKWayPartitioning from kmetis.c:
/// - For ncuts iterations: coarsen -> alloc kway memory -> init partition -> refine -> keep best
/// - Best-cut selection uses ComputeLoadImbalanceDiff with curbal <= 0.0005 criterion
pub fn mlevel_kway_partitioning(
    ctrl: &mut Control,
    graph: &mut GraphData,
    part: &mut [Idx],
) -> Idx {
    let nparts = ctrl.num_parts;
    let ncuts = ctrl.num_cuts;
    let num_vertices = graph.num_vertices as usize;

    let mut bestobj: Idx = 0;
    let mut bestbal: Real = 0.0;
    let mut best_part = vec![0 as Idx; num_vertices];

    // Save graph dimensions for pool size calculation
    let graph_num_edges = graph.num_edges;
    let graph_num_vertices = graph.num_vertices;

    for i in 0..ncuts {
        // Coarsen the graph — returns arena of coarser levels
        let mut levels = coarsen_graph(ctrl, graph);

        // AllocateKWayPartitionMemory on the coarsest graph
        let coarsest = levels.last_mut().unwrap();
        let cnum_vertices = coarsest.num_vertices as usize;
        coarsest.partition = vec![0; cnum_vertices];
        coarsest.part_weights = vec![0; nparts as usize * coarsest.num_constraints as usize];
        coarsest.boundary_map = vec![-1; cnum_vertices];
        coarsest.boundary_list = vec![0; cnum_vertices];

        // InitKWayPartitioning: compute initial partition via recursive bisection
        pmetis::init_kway_partition(ctrl, levels.last_mut().unwrap(), nparts);

        // AllocateRefinementWorkSpace: initialize neighbor_pool
        let pool_size = (2 * graph_num_edges as usize + graph_num_vertices as usize).max(1024);
        ctrl.init_neighbor_pool(pool_size);

        // RefineKWay: walk from coarsest to finest with k-way refinement
        kwayrefine::refine_kway(ctrl, graph, &mut levels);

        // Evaluate the result on the finest graph
        let curobj = graph.edge_cut;
        let curbal = compute_load_imbalance_diff(
            graph,
            nparts as usize,
            &ctrl.partition_ij_balance_multipliers,
            &ctrl.imbalance_tols,
        );
        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            best_part.copy_from_slice(&graph.partition);
            bestobj = curobj;
            bestbal = curbal;
        }

        // FreeRData: clean up partition data from the graph for the next iteration
        free_rdata(graph);

        if bestobj == 0 {
            break;
        }
    }

    // Copy best result to output
    part.copy_from_slice(&best_part);
    bestobj
}

/// ComputeLoadImbalanceDiff for k-way.
fn compute_load_imbalance_diff(
    graph: &GraphData,
    nparts: usize,
    partition_ij_balance_multipliers: &[Real],
    imbalance_tols: &[Real],
) -> Real {
    let ncon = graph.num_constraints as usize;
    let mut max_diff = Real::NEG_INFINITY;
    for i in 0..nparts {
        for j in 0..ncon {
            let idx = i * ncon + j;
            if idx < graph.part_weights.len()
                && idx < partition_ij_balance_multipliers.len()
                && j < imbalance_tols.len()
            {
                let diff = graph.part_weights[idx] as Real * partition_ij_balance_multipliers[idx]
                    - imbalance_tols[j];
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }
    max_diff
}

/// Free partition refinement data from the graph.
/// In C METIS this is FreeRData which frees partition, part_weights, boundary_map, boundary_list, kway_refinement_info, etc.
fn free_rdata(graph: &mut GraphData) {
    graph.partition = Vec::new();
    graph.part_weights = Vec::new();
    graph.boundary_map = Vec::new();
    graph.boundary_list = Vec::new();
    graph.kway_refinement_info = Vec::new();
    graph.internal_degree = Vec::new();
    graph.external_degree = Vec::new();
    graph.num_boundary = 0;
    graph.edge_cut = 0;
}
