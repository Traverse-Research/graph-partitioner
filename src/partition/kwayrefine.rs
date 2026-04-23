use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::{GraphData, KwayCutInfo};
use crate::partition::kwayfm;

const BOUNDARY_REFINE: Idx = 1;
const BOUNDARY_BALANCE: Idx = 2;

/// RefineKWay: entry point of k-way cut-based refinement.
///
/// Matches C METIS RefineKWay from kwayrefine.c:
/// - orggraph is a raw pointer to the original (finest) graph
/// - graph is the coarsest graph (starting point)
/// - Walks from coarsest to finest: compute params, balance if needed, refine, project
pub fn refine_kway(ctrl: &mut Control, orggraph: *mut GraphData, graph: &mut GraphData) {
    let niter = ctrl.num_iter;
    let nparts = ctrl.num_parts;

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
            compute_kway_boundary(g, nparts, BOUNDARY_BALANCE);
            kwayfm::greedy_kway_cut_optimize(ctrl, g, 1, kwayfm::MODE_BALANCE);
            compute_kway_boundary(g, nparts, BOUNDARY_REFINE);
        }

        kwayfm::greedy_kway_cut_optimize(ctrl, g, niter, kwayfm::MODE_REFINE);

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
        compute_kway_boundary(g, nparts, BOUNDARY_BALANCE);
        kwayfm::greedy_kway_cut_optimize(ctrl, g, 10, kwayfm::MODE_BALANCE);
        compute_kway_boundary(g, nparts, BOUNDARY_REFINE);
        kwayfm::greedy_kway_cut_optimize(ctrl, g, niter, kwayfm::MODE_REFINE);
    }
}

/// ComputeKWayPartitionParams: compute initial id/ed, boundary, and kway_refinement_info.
///
/// Matches C METIS ComputeKWayPartitionParams for METIS_OBJTYPE_CUT.
pub fn compute_kway_partition_params(ctrl: &mut Control, graph: &mut GraphData, nparts: Idx) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;
    let np = nparts as usize;

    // Allocate partition memory
    if graph.partition.len() != num_vertices {
        graph.partition = vec![0; num_vertices];
    }
    graph.part_weights = vec![0; np * ncon];
    graph.boundary_map = vec![-1; num_vertices];
    graph.boundary_list = vec![0; num_vertices];
    graph.kway_refinement_info = vec![KwayCutInfo::default(); num_vertices];

    // Reset neighbor_pool
    ctrl.reset_neighbor_pool();

    // Compute part_weights
    for i in 0..num_vertices {
        let me = graph.partition[i] as usize;
        if ncon == 1 {
            graph.part_weights[me] += graph.vertex_weights[i];
        } else {
            for j in 0..ncon {
                graph.part_weights[me * ncon + j] += graph.vertex_weights[i * ncon + j];
            }
        }
    }

    let mut num_boundary: Idx = 0;
    let mut edge_cut: Idx = 0;

    for i in 0..num_vertices {
        let me = graph.partition[i] as usize;
        let mut id: Idx = 0;
        let mut ed: Idx = 0;

        for j in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let other = graph.partition[graph.adjacency[j] as usize] as usize;
            if me == other {
                id += graph.edge_weights[j];
            } else {
                ed += graph.edge_weights[j];
            }
        }

        graph.kway_refinement_info[i].internal_degree = id;
        graph.kway_refinement_info[i].external_degree = ed;

        if ed > 0 {
            edge_cut += ed;

            // Allocate neighbor info from pool
            let adj_count = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
            let inbr = ctrl.alloc_neighbor_info(adj_count);
            graph.kway_refinement_info[i].neighbor_offset = inbr;

            // Compute per-partition external degrees
            let mut nnbrs: Idx = 0;
            for j in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let other = graph.partition[graph.adjacency[j] as usize];
                if me as Idx != other {
                    // Search existing neighbors
                    let mut found = false;
                    for k in 0..nnbrs as usize {
                        if ctrl.neighbor_pool[inbr as usize + k].part_id == other {
                            ctrl.neighbor_pool[inbr as usize + k].external_degree += graph.edge_weights[j];
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].part_id = other;
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].external_degree = graph.edge_weights[j];
                        nnbrs += 1;
                    }
                }
            }
            graph.kway_refinement_info[i].num_neighbors = nnbrs;

            // Only ed-id >= 0 nodes are considered boundary (BOUNDARY_REFINE)
            if ed - id >= 0 {
                graph.boundary_list[num_boundary as usize] = i as Idx;
                graph.boundary_map[i] = num_boundary;
                num_boundary += 1;
            }
        } else {
            graph.kway_refinement_info[i].neighbor_offset = -1;
        }
    }

    graph.edge_cut = edge_cut / 2;
    graph.num_boundary = num_boundary;
}

/// ProjectKWayPartition: project partition from coarser to finer graph.
///
/// Matches C METIS ProjectKWayPartition for METIS_OBJTYPE_CUT.
fn project_kway_partition(ctrl: &mut Control, graph: &mut GraphData) {
    let nparts = ctrl.num_parts as usize;

    // Extract data from coarser graph before mutating
    let (cwhere, ced_info, cedge_cut, cpart_weights) = {
        let cgraph = graph.coarser.as_ref()
            .expect("project_kway_partition: coarser graph is None");
        // For optimization: coarse_map[i] will store kway_refinement_info[coarse_map[i]].external_degree
        let ced: Vec<Idx> = (0..cgraph.num_vertices as usize)
            .map(|v| cgraph.kway_refinement_info[v].external_degree)
            .collect();
        (
            cgraph.partition.clone(),
            ced,
            cgraph.edge_cut,
            cgraph.part_weights.clone(),
        )
    };

    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Project partition and store ed optimization hint
    // Allocate kway partition memory
    graph.part_weights = vec![0; nparts * ncon];
    graph.boundary_map = vec![-1; num_vertices];
    graph.boundary_list = vec![0; num_vertices];
    graph.kway_refinement_info = vec![KwayCutInfo::default(); num_vertices];

    // partition[i] = cwhere[coarse_map[i]], coarse_map[i] = ced[coarse_map[i]] (optimization)
    let mut cmap_ed = vec![0 as Idx; num_vertices];
    if graph.partition.len() != num_vertices {
        graph.partition = vec![0; num_vertices];
    }
    for i in 0..num_vertices {
        let k = graph.coarse_map[i] as usize;
        graph.partition[i] = cwhere[k];
        cmap_ed[i] = ced_info[k]; // Store kway_refinement_info[coarse_map[i]].external_degree for optimization
    }

    // Reset neighbor_pool
    ctrl.reset_neighbor_pool();

    // Hash table for per-partition neighbor accumulation
    let mut htable = vec![-1 as Idx; nparts];

    let mut num_boundary: Idx = 0;

    for i in 0..num_vertices {
        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;

        if cmap_ed[i] == 0 {
            // Interior node: coarse vertex had ed==0, all neighbors same partition
            let mut tid: Idx = 0;
            for j in istart..iend {
                tid += graph.edge_weights[j];
            }
            graph.kway_refinement_info[i].internal_degree = tid;
            graph.kway_refinement_info[i].neighbor_offset = -1;
        } else {
            // Potentially interface node
            let adj_count = iend - istart;
            let inbr = ctrl.alloc_neighbor_info(adj_count);
            graph.kway_refinement_info[i].neighbor_offset = inbr;

            let me = graph.partition[i];
            let mut tid: Idx = 0;
            let mut ted: Idx = 0;
            let mut nnbrs: Idx = 0;

            for j in istart..iend {
                let other = graph.partition[graph.adjacency[j] as usize];
                if me == other {
                    tid += graph.edge_weights[j];
                } else {
                    ted += graph.edge_weights[j];
                    let hk = htable[other as usize];
                    if hk == -1 {
                        htable[other as usize] = nnbrs;
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].part_id = other;
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].external_degree = graph.edge_weights[j];
                        nnbrs += 1;
                    } else {
                        ctrl.neighbor_pool[inbr as usize + hk as usize].external_degree += graph.edge_weights[j];
                    }
                }
            }

            graph.kway_refinement_info[i].internal_degree = tid;
            graph.kway_refinement_info[i].external_degree = ted;

            if ted == 0 {
                // Remove space: this was interior after all
                // In C METIS: ctrl->nbrpoolcpos -= min(nparts, iend-istart)
                graph.kway_refinement_info[i].neighbor_offset = -1;
            } else {
                graph.kway_refinement_info[i].num_neighbors = nnbrs;

                if ted - tid >= 0 {
                    graph.boundary_list[num_boundary as usize] = i as Idx;
                    graph.boundary_map[i] = num_boundary;
                    num_boundary += 1;
                }

                // Reset htable
                for k in 0..nnbrs as usize {
                    htable[ctrl.neighbor_pool[inbr as usize + k].part_id as usize] = -1;
                }
            }
        }
    }

    graph.num_boundary = num_boundary;
    graph.edge_cut = cedge_cut;

    // Copy part_weights from coarser
    let len = cpart_weights.len().min(graph.part_weights.len());
    graph.part_weights[..len].copy_from_slice(&cpart_weights[..len]);

    // Free coarser graph
    graph.coarser = None;
}

/// ComputeKWayBoundary: recompute boundary based on bndtype.
///
/// Matches C METIS ComputeKWayBoundary for METIS_OBJTYPE_CUT.
fn compute_kway_boundary(graph: &mut GraphData, _nparts: Idx, bndtype: Idx) {
    let num_vertices = graph.num_vertices as usize;

    graph.boundary_map = vec![-1; num_vertices];
    // boundary_list already allocated
    if graph.boundary_list.len() < num_vertices {
        graph.boundary_list.resize(num_vertices, 0);
    }

    let mut num_boundary: Idx = 0;

    if bndtype == BOUNDARY_REFINE {
        for i in 0..num_vertices {
            if graph.kway_refinement_info[i].external_degree > 0 && graph.kway_refinement_info[i].external_degree - graph.kway_refinement_info[i].internal_degree >= 0 {
                graph.boundary_list[num_boundary as usize] = i as Idx;
                graph.boundary_map[i] = num_boundary;
                num_boundary += 1;
            }
        }
    } else {
        // BOUNDARY_BALANCE
        for i in 0..num_vertices {
            if graph.kway_refinement_info[i].external_degree > 0 {
                graph.boundary_list[num_boundary as usize] = i as Idx;
                graph.boundary_map[i] = num_boundary;
                num_boundary += 1;
            }
        }
    }

    graph.num_boundary = num_boundary;
}

/// IsBalanced: check if partition weights are within balance constraints.
///
/// Matches C METIS IsBalanced = (ComputeLoadImbalanceDiff <= ffactor).
fn is_balanced(ctrl: &Control, graph: &GraphData, ffactor: Real) -> bool {
    compute_load_imbalance_diff_kway(graph, ctrl.num_parts as usize, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols) <= ffactor
}

/// ComputeLoadImbalanceDiff for k-way.
fn compute_load_imbalance_diff_kway(
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
            if idx < graph.part_weights.len() && idx < partition_ij_balance_multipliers.len() && j < imbalance_tols.len() {
                let diff = graph.part_weights[idx] as Real * partition_ij_balance_multipliers[idx] - imbalance_tols[j];
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }
    max_diff
}
