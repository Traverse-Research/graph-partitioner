use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::{GraphData, CKRInfo};
use crate::partition::kwayfm;

const BNDTYPE_REFINE: Idx = 1;
const BNDTYPE_BALANCE: Idx = 2;

/// RefineKWay: entry point of k-way cut-based refinement.
///
/// Matches C METIS RefineKWay from kwayrefine.c:
/// - orggraph is a raw pointer to the original (finest) graph
/// - graph is the coarsest graph (starting point)
/// - Walks from coarsest to finest: compute params, balance if needed, refine, project
pub fn refine_kway(ctrl: &mut Ctrl, orggraph: *mut GraphData, graph: &mut GraphData) {
    let niter = ctrl.niter;
    let nparts = ctrl.nparts;

    // Determine how many levels there are
    let mut nlevels = 0;
    {
        let mut ptr: *mut GraphData = graph as *mut GraphData;
        while ptr != orggraph {
            let g = unsafe { &*ptr };
            ptr = g.finer;
            nlevels += 1;
        }
    }

    // Compute the parameters of the coarsest graph
    compute_kway_partition_params(ctrl, graph, nparts);

    // Refine each successively finer graph
    let mut cur: *mut GraphData = graph as *mut GraphData;
    let mut i = 0;

    loop {
        let g = unsafe { &mut *cur };

        // Balance check at deeper levels
        if 2 * i >= nlevels && !is_balanced(ctrl, g, 0.02) {
            compute_kway_boundary(g, nparts, BNDTYPE_BALANCE);
            kwayfm::greedy_kway_cut_optimize(ctrl, g, 1, kwayfm::OMODE_BALANCE);
            compute_kway_boundary(g, nparts, BNDTYPE_REFINE);
        }

        kwayfm::greedy_kway_cut_optimize(ctrl, g, niter, kwayfm::OMODE_REFINE);

        if cur == orggraph {
            break;
        }

        let finer_ptr = g.finer;
        assert!(!finer_ptr.is_null(), "refine_kway: finer pointer is null before reaching orggraph");

        // Project partition to finer level
        project_kway_partition(ctrl, unsafe { &mut *finer_ptr });

        cur = finer_ptr;
        i += 1;
    }

    // Final balance check at the finest level
    let g = unsafe { &mut *cur };
    if !is_balanced(ctrl, g, 0.0) {
        compute_kway_boundary(g, nparts, BNDTYPE_BALANCE);
        kwayfm::greedy_kway_cut_optimize(ctrl, g, 10, kwayfm::OMODE_BALANCE);
        compute_kway_boundary(g, nparts, BNDTYPE_REFINE);
        kwayfm::greedy_kway_cut_optimize(ctrl, g, niter, kwayfm::OMODE_REFINE);
    }
}

/// ComputeKWayPartitionParams: compute initial id/ed, boundary, and ckrinfo.
///
/// Matches C METIS ComputeKWayPartitionParams for METIS_OBJTYPE_CUT.
pub fn compute_kway_partition_params(ctrl: &mut Ctrl, graph: &mut GraphData, nparts: Idx) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    let np = nparts as usize;

    // Allocate partition memory
    if graph.where_.len() != nvtxs {
        graph.where_ = vec![0; nvtxs];
    }
    graph.pwgts = vec![0; np * ncon];
    graph.bndptr = vec![-1; nvtxs];
    graph.bndind = vec![0; nvtxs];
    graph.ckrinfo = vec![CKRInfo::default(); nvtxs];

    // Reset cnbrpool
    ctrl.cnbrpool_reset();

    // Compute pwgts
    for i in 0..nvtxs {
        let me = graph.where_[i] as usize;
        if ncon == 1 {
            graph.pwgts[me] += graph.vwgt[i];
        } else {
            for j in 0..ncon {
                graph.pwgts[me * ncon + j] += graph.vwgt[i * ncon + j];
            }
        }
    }

    let mut nbnd: Idx = 0;
    let mut mincut: Idx = 0;

    for i in 0..nvtxs {
        let me = graph.where_[i] as usize;
        let mut id: Idx = 0;
        let mut ed: Idx = 0;

        for j in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let other = graph.where_[graph.adjncy[j] as usize] as usize;
            if me == other {
                id += graph.adjwgt[j];
            } else {
                ed += graph.adjwgt[j];
            }
        }

        graph.ckrinfo[i].id = id;
        graph.ckrinfo[i].ed = ed;

        if ed > 0 {
            mincut += ed;

            // Allocate neighbor info from pool
            let adj_count = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
            let inbr = ctrl.cnbrpool_get_next(adj_count);
            graph.ckrinfo[i].inbr = inbr;

            // Compute per-partition external degrees
            let mut nnbrs: Idx = 0;
            for j in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let other = graph.where_[graph.adjncy[j] as usize];
                if me as Idx != other {
                    // Search existing neighbors
                    let mut found = false;
                    for k in 0..nnbrs as usize {
                        if ctrl.cnbrpool[inbr as usize + k].pid == other {
                            ctrl.cnbrpool[inbr as usize + k].ed += graph.adjwgt[j];
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        ctrl.cnbrpool[inbr as usize + nnbrs as usize].pid = other;
                        ctrl.cnbrpool[inbr as usize + nnbrs as usize].ed = graph.adjwgt[j];
                        nnbrs += 1;
                    }
                }
            }
            graph.ckrinfo[i].nnbrs = nnbrs;

            // Only ed-id >= 0 nodes are considered boundary (BNDTYPE_REFINE)
            if ed - id >= 0 {
                graph.bndind[nbnd as usize] = i as Idx;
                graph.bndptr[i] = nbnd;
                nbnd += 1;
            }
        } else {
            graph.ckrinfo[i].inbr = -1;
        }
    }

    graph.mincut = mincut / 2;
    graph.nbnd = nbnd;
}

/// ProjectKWayPartition: project partition from coarser to finer graph.
///
/// Matches C METIS ProjectKWayPartition for METIS_OBJTYPE_CUT.
fn project_kway_partition(ctrl: &mut Ctrl, graph: &mut GraphData) {
    let nparts = ctrl.nparts as usize;

    // Extract data from coarser graph before mutating
    let (cwhere, ced_info, cmincut, cpwgts) = {
        let cgraph = graph.coarser.as_ref()
            .expect("project_kway_partition: coarser graph is None");
        // For optimization: cmap[i] will store ckrinfo[cmap[i]].ed
        let ced: Vec<Idx> = (0..cgraph.nvtxs as usize)
            .map(|v| cgraph.ckrinfo[v].ed)
            .collect();
        (
            cgraph.where_.clone(),
            ced,
            cgraph.mincut,
            cgraph.pwgts.clone(),
        )
    };

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Project partition and store ed optimization hint
    // Allocate kway partition memory
    graph.pwgts = vec![0; nparts * ncon];
    graph.bndptr = vec![-1; nvtxs];
    graph.bndind = vec![0; nvtxs];
    graph.ckrinfo = vec![CKRInfo::default(); nvtxs];

    // where_[i] = cwhere[cmap[i]], cmap[i] = ced[cmap[i]] (optimization)
    let mut cmap_ed = vec![0 as Idx; nvtxs];
    if graph.where_.len() != nvtxs {
        graph.where_ = vec![0; nvtxs];
    }
    for i in 0..nvtxs {
        let k = graph.cmap[i] as usize;
        graph.where_[i] = cwhere[k];
        cmap_ed[i] = ced_info[k]; // Store ckrinfo[cmap[i]].ed for optimization
    }

    // Reset cnbrpool
    ctrl.cnbrpool_reset();

    // Hash table for per-partition neighbor accumulation
    let mut htable = vec![-1 as Idx; nparts];

    let mut nbnd: Idx = 0;

    for i in 0..nvtxs {
        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;

        if cmap_ed[i] == 0 {
            // Interior node: coarse vertex had ed==0, all neighbors same partition
            let mut tid: Idx = 0;
            for j in istart..iend {
                tid += graph.adjwgt[j];
            }
            graph.ckrinfo[i].id = tid;
            graph.ckrinfo[i].inbr = -1;
        } else {
            // Potentially interface node
            let adj_count = iend - istart;
            let inbr = ctrl.cnbrpool_get_next(adj_count);
            graph.ckrinfo[i].inbr = inbr;

            let me = graph.where_[i];
            let mut tid: Idx = 0;
            let mut ted: Idx = 0;
            let mut nnbrs: Idx = 0;

            for j in istart..iend {
                let other = graph.where_[graph.adjncy[j] as usize];
                if me == other {
                    tid += graph.adjwgt[j];
                } else {
                    ted += graph.adjwgt[j];
                    let hk = htable[other as usize];
                    if hk == -1 {
                        htable[other as usize] = nnbrs;
                        ctrl.cnbrpool[inbr as usize + nnbrs as usize].pid = other;
                        ctrl.cnbrpool[inbr as usize + nnbrs as usize].ed = graph.adjwgt[j];
                        nnbrs += 1;
                    } else {
                        ctrl.cnbrpool[inbr as usize + hk as usize].ed += graph.adjwgt[j];
                    }
                }
            }

            graph.ckrinfo[i].id = tid;
            graph.ckrinfo[i].ed = ted;

            if ted == 0 {
                // Remove space: this was interior after all
                // In C METIS: ctrl->nbrpoolcpos -= min(nparts, iend-istart)
                graph.ckrinfo[i].inbr = -1;
            } else {
                graph.ckrinfo[i].nnbrs = nnbrs;

                if ted - tid >= 0 {
                    graph.bndind[nbnd as usize] = i as Idx;
                    graph.bndptr[i] = nbnd;
                    nbnd += 1;
                }

                // Reset htable
                for k in 0..nnbrs as usize {
                    htable[ctrl.cnbrpool[inbr as usize + k].pid as usize] = -1;
                }
            }
        }
    }

    graph.nbnd = nbnd;
    graph.mincut = cmincut;

    // Copy pwgts from coarser
    let len = cpwgts.len().min(graph.pwgts.len());
    graph.pwgts[..len].copy_from_slice(&cpwgts[..len]);

    // Free coarser graph
    graph.coarser = None;
}

/// ComputeKWayBoundary: recompute boundary based on bndtype.
///
/// Matches C METIS ComputeKWayBoundary for METIS_OBJTYPE_CUT.
fn compute_kway_boundary(graph: &mut GraphData, _nparts: Idx, bndtype: Idx) {
    let nvtxs = graph.nvtxs as usize;

    graph.bndptr = vec![-1; nvtxs];
    // bndind already allocated
    if graph.bndind.len() < nvtxs {
        graph.bndind.resize(nvtxs, 0);
    }

    let mut nbnd: Idx = 0;

    if bndtype == BNDTYPE_REFINE {
        for i in 0..nvtxs {
            if graph.ckrinfo[i].ed > 0 && graph.ckrinfo[i].ed - graph.ckrinfo[i].id >= 0 {
                graph.bndind[nbnd as usize] = i as Idx;
                graph.bndptr[i] = nbnd;
                nbnd += 1;
            }
        }
    } else {
        // BNDTYPE_BALANCE
        for i in 0..nvtxs {
            if graph.ckrinfo[i].ed > 0 {
                graph.bndind[nbnd as usize] = i as Idx;
                graph.bndptr[i] = nbnd;
                nbnd += 1;
            }
        }
    }

    graph.nbnd = nbnd;
}

/// IsBalanced: check if partition weights are within balance constraints.
///
/// Matches C METIS IsBalanced = (ComputeLoadImbalanceDiff <= ffactor).
fn is_balanced(ctrl: &Ctrl, graph: &GraphData, ffactor: Real) -> bool {
    compute_load_imbalance_diff_kway(graph, ctrl.nparts as usize, &ctrl.pijbm, &ctrl.ubfactors) <= ffactor
}

/// ComputeLoadImbalanceDiff for k-way.
fn compute_load_imbalance_diff_kway(
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
