use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::{GraphData, KwayCutInfo};
use crate::partition::kwayfm;

const BOUNDARY_REFINE: Idx = 1;
const BOUNDARY_BALANCE: Idx = 2;

/// RefineKWay: entry point of k-way cut-based refinement.
///
/// `graph` is the finest level. `levels` contains the coarser levels where
/// `levels[0]` is the first coarser level and `levels[last]` is the coarsest
/// (where initial partitioning was performed).
///
/// Walks from coarsest to finest: compute params, balance if needed, refine, project.
pub fn refine_kway(ctrl: &mut Control, graph: &mut GraphData, levels: &mut Vec<GraphData>) {
    let niter = ctrl.num_iter;
    let nparts = ctrl.num_parts;
    let nlevels = levels.len(); // number of coarser levels

    // Compute the parameters of the coarsest graph
    compute_kway_partition_params(ctrl, &mut levels[nlevels - 1], nparts);

    // Refine each level from coarsest to finest within the arena
    // i goes from 0 (coarsest = levels[nlevels-1]) to nlevels-1 (levels[0])
    for step in 0..nlevels {
        let level_idx = nlevels - 1 - step;

        // Balance check at deeper levels (closer to finest)
        if 2 * step >= nlevels && !is_balanced(ctrl, &levels[level_idx], 0.02) {
            compute_kway_boundary(&mut levels[level_idx], nparts, BOUNDARY_BALANCE);
            kwayfm::greedy_kway_cut_optimize(ctrl, &mut levels[level_idx], 1, kwayfm::MODE_BALANCE);
            compute_kway_boundary(&mut levels[level_idx], nparts, BOUNDARY_REFINE);
        }

        kwayfm::greedy_kway_cut_optimize(ctrl, &mut levels[level_idx], niter, kwayfm::MODE_REFINE);

        // Project to finer level if not at the finest arena level
        if level_idx > 0 {
            let (left, right) = levels.split_at_mut(level_idx);
            let finer = &mut left[level_idx - 1];
            let coarser = &right[0]; // = levels[level_idx]
            project_kway_partition(ctrl, finer, coarser);
        }
    }

    // Final projection: from levels[0] to graph (finest)
    project_kway_partition(ctrl, graph, &levels[0]);

    // Balance check at finest level
    if 2 * nlevels >= nlevels && !is_balanced(ctrl, graph, 0.02) {
        compute_kway_boundary(graph, nparts, BOUNDARY_BALANCE);
        kwayfm::greedy_kway_cut_optimize(ctrl, graph, 1, kwayfm::MODE_BALANCE);
        compute_kway_boundary(graph, nparts, BOUNDARY_REFINE);
    }

    kwayfm::greedy_kway_cut_optimize(ctrl, graph, niter, kwayfm::MODE_REFINE);

    // Final balance check at the finest level
    if !is_balanced(ctrl, graph, 0.0) {
        compute_kway_boundary(graph, nparts, BOUNDARY_BALANCE);
        kwayfm::greedy_kway_cut_optimize(ctrl, graph, 10, kwayfm::MODE_BALANCE);
        compute_kway_boundary(graph, nparts, BOUNDARY_REFINE);
        kwayfm::greedy_kway_cut_optimize(ctrl, graph, niter, kwayfm::MODE_REFINE);
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

        let ji_start = graph.xadj[i] as usize;
        let ji_end = graph.xadj[i + 1] as usize;
        let adj_slice = &graph.adjacency[ji_start..ji_end];
        let ewgt_slice = &graph.edge_weights[ji_start..ji_end];

        for (&neighbor, &ew) in adj_slice.iter().zip(ewgt_slice) {
            let other = graph.partition[neighbor as usize] as usize;
            if me == other {
                id += ew;
            } else {
                ed += ew;
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
/// `finer` is the finer level; `coarser` is the coarser level.
/// Uses `finer.coarse_map` to map fine vertices to coarse vertices.
fn project_kway_partition(ctrl: &mut Control, finer: &mut GraphData, coarser: &GraphData) {
    let nparts = ctrl.num_parts as usize;

    // Extract optimization hint: coarse vertex external degree
    let ced: Vec<Idx> = (0..coarser.num_vertices as usize)
        .map(|v| coarser.kway_refinement_info[v].external_degree)
        .collect();

    let num_vertices = finer.num_vertices as usize;
    let ncon = finer.num_constraints as usize;

    // Allocate kway partition memory
    finer.part_weights = vec![0; nparts * ncon];
    finer.boundary_map = vec![-1; num_vertices];
    finer.boundary_list = vec![0; num_vertices];
    finer.kway_refinement_info = vec![KwayCutInfo::default(); num_vertices];

    // Project partition and store ed optimization hint
    let mut cmap_ed = vec![0 as Idx; num_vertices];
    if finer.partition.len() != num_vertices {
        finer.partition = vec![0; num_vertices];
    }
    for i in 0..num_vertices {
        let k = finer.coarse_map[i] as usize;
        finer.partition[i] = coarser.partition[k];
        cmap_ed[i] = ced[k];
    }

    // Reset neighbor_pool
    ctrl.reset_neighbor_pool();

    // Hash table for per-partition neighbor accumulation
    let mut htable = vec![-1 as Idx; nparts];

    let mut num_boundary: Idx = 0;

    for i in 0..num_vertices {
        let istart = finer.xadj[i] as usize;
        let iend = finer.xadj[i + 1] as usize;

        if cmap_ed[i] == 0 {
            // Interior node: coarse vertex had ed==0, all neighbors same partition
            let mut tid: Idx = 0;
            for j in istart..iend {
                tid += finer.edge_weights[j];
            }
            finer.kway_refinement_info[i].internal_degree = tid;
            finer.kway_refinement_info[i].neighbor_offset = -1;
        } else {
            // Potentially interface node
            let adj_count = iend - istart;
            let inbr = ctrl.alloc_neighbor_info(adj_count);
            finer.kway_refinement_info[i].neighbor_offset = inbr;

            let me = finer.partition[i];
            let mut tid: Idx = 0;
            let mut ted: Idx = 0;
            let mut nnbrs: Idx = 0;

            for j in istart..iend {
                let other = finer.partition[finer.adjacency[j] as usize];
                if me == other {
                    tid += finer.edge_weights[j];
                } else {
                    ted += finer.edge_weights[j];
                    let hk = htable[other as usize];
                    if hk == -1 {
                        htable[other as usize] = nnbrs;
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].part_id = other;
                        ctrl.neighbor_pool[inbr as usize + nnbrs as usize].external_degree = finer.edge_weights[j];
                        nnbrs += 1;
                    } else {
                        ctrl.neighbor_pool[inbr as usize + hk as usize].external_degree += finer.edge_weights[j];
                    }
                }
            }

            finer.kway_refinement_info[i].internal_degree = tid;
            finer.kway_refinement_info[i].external_degree = ted;

            if ted == 0 {
                // Remove space: this was interior after all
                finer.kway_refinement_info[i].neighbor_offset = -1;
            } else {
                finer.kway_refinement_info[i].num_neighbors = nnbrs;

                if ted - tid >= 0 {
                    finer.boundary_list[num_boundary as usize] = i as Idx;
                    finer.boundary_map[i] = num_boundary;
                    num_boundary += 1;
                }

                // Reset htable
                for k in 0..nnbrs as usize {
                    htable[ctrl.neighbor_pool[inbr as usize + k].part_id as usize] = -1;
                }
            }
        }
    }

    finer.num_boundary = num_boundary;
    finer.edge_cut = coarser.edge_cut;

    // Copy part_weights from coarser
    let len = coarser.part_weights.len().min(finer.part_weights.len());
    finer.part_weights[..len].copy_from_slice(&coarser.part_weights[..len]);
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
