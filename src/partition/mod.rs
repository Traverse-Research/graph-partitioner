pub mod kway;
pub mod pmetis;
pub mod initpart;
pub mod fm;
pub mod balance;
pub mod refine2way;
pub mod split;
pub mod kwayrefine;
pub mod kwayfm;

use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::graph::coarsen::coarsen_graph;

/// MlevelKWayPartitioning: multilevel k-way partitioning orchestration.
///
/// Matches C METIS MlevelKWayPartitioning from kmetis.c:
/// - For ncuts iterations: coarsen -> alloc kway memory -> init partition -> refine -> keep best
/// - Best-cut selection uses ComputeLoadImbalanceDiff with curbal <= 0.0005 criterion
pub fn mlevel_kway_partitioning(ctrl: &mut Ctrl, graph: &mut GraphData, part: &mut [Idx]) -> Idx {
    let nparts = ctrl.nparts;
    let ncuts = ctrl.ncuts;
    let nvtxs = graph.nvtxs as usize;

    let mut bestobj: Idx = 0;
    let mut bestbal: Real = 0.0;
    let mut best_part = vec![0 as Idx; nvtxs];

    // Save graph dimensions for pool size calculation
    let graph_nedges = graph.nedges;
    let graph_nvtxs = graph.nvtxs;

    for i in 0..ncuts {
        // Coarsen the graph
        coarsen_graph(ctrl, graph);

        // Use raw pointers to handle the borrow checker (we need both orggraph and cgraph
        // from the same coarsening chain simultaneously)
        let orggraph_ptr: *mut GraphData = graph as *mut GraphData;
        let cgraph_ptr: *mut GraphData = get_coarsest_ptr(orggraph_ptr);
        let cgraph = unsafe { &mut *cgraph_ptr };

        // AllocateKWayPartitionMemory on the coarsest graph
        let cnvtxs = cgraph.nvtxs as usize;
        cgraph.where_ = vec![0; cnvtxs];
        cgraph.pwgts = vec![0; nparts as usize * cgraph.ncon as usize];
        cgraph.bndptr = vec![-1; cnvtxs];
        cgraph.bndind = vec![0; cnvtxs];

        // InitKWayPartitioning: compute initial partition via recursive bisection
        pmetis::init_kway_partition(ctrl, cgraph, nparts);
        // AllocateRefinementWorkSpace: initialize cnbrpool
        let pool_size = (2 * graph_nedges as usize + graph_nvtxs as usize).max(1024);
        ctrl.cnbrpool_init(pool_size);

        // RefineKWay: walk from coarsest to finest with k-way refinement
        kwayrefine::refine_kway(ctrl, orggraph_ptr, cgraph);

        // Evaluate the result on the finest graph (safe to use graph again since
        // RefineKWay has already dropped all coarser graphs)
        let curobj = graph.mincut;
        let curbal = compute_load_imbalance_diff(graph, nparts as usize, &ctrl.pijbm, &ctrl.ubfactors);
        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            best_part.copy_from_slice(&graph.where_);
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
    pijbm: &[Real],
    ubfactors: &[Real],
) -> Real {
    let ncon = graph.ncon as usize;
    let mut max_diff = Real::NEG_INFINITY;
    for i in 0..nparts {
        for j in 0..ncon {
            let idx = i * ncon + j;
            if idx < graph.pwgts.len() && idx < pijbm.len() && j < ubfactors.len() {
                let diff = graph.pwgts[idx] as Real * pijbm[idx] - ubfactors[j];
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }
    max_diff
}

/// Traverse coarsening chain to find coarsest graph (raw pointer version).
fn get_coarsest_ptr(g: *mut GraphData) -> *mut GraphData {
    let mut cur = g;
    unsafe {
        while (*cur).coarser.is_some() {
            cur = &mut **(*cur).coarser.as_mut().unwrap() as *mut GraphData;
        }
    }
    cur
}

/// Free partition refinement data from the graph (and its coarsening chain).
/// In C METIS this is FreeRData which frees where_, pwgts, bndptr, bndind, ckrinfo, etc.
fn free_rdata(graph: &mut GraphData) {
    graph.where_ = Vec::new();
    graph.pwgts = Vec::new();
    graph.bndptr = Vec::new();
    graph.bndind = Vec::new();
    graph.ckrinfo = Vec::new();
    graph.id = Vec::new();
    graph.ed = Vec::new();
    graph.nbnd = 0;
    graph.mincut = 0;

    // Also need to free coarsening chain graphs
    // After RefineKWay, the coarsening chain should already be freed
    // (ProjectKWayPartition drops each coarser graph as it goes)
    // But if coarser still exists, drop it
    graph.coarser = None;
}
