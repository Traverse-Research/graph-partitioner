use super::balance::{
    better_balance_2way, compute_load_imbalance_diff_vec, iargmax_nrm, select_queue,
};
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::types::{Idx, Real};
use crate::util::pqueue::PQueue;

/// FM 2-way cut refinement with rollback (single constraint).
/// Matches METIS FM_2WayCutRefine.
#[allow(unused_assignments)]
pub fn fm_2way_cut_refine(
    ctrl: &mut Control,
    graph: &mut GraphData,
    target_part_weights: &[Real],
    niter: Idx,
) {
    let num_vertices = graph.num_vertices as usize;
    if graph.num_boundary <= 0 {
        return;
    }

    // Integer target weights
    let itarget_part_weights = [
        (graph.total_vertex_weight[0] as Real * target_part_weights[0]) as Idx,
        graph.total_vertex_weight[0]
            - (graph.total_vertex_weight[0] as Real * target_part_weights[0]) as Idx,
    ];

    let limit = (0.01 * num_vertices as f64).max(15.0).min(100.0) as usize;
    let avgvertex_weights = ((graph.part_weights[0] + graph.part_weights[1]) / 20)
        .min(2 * (graph.part_weights[0] + graph.part_weights[1]) / num_vertices.max(1) as Idx);
    let origdiff = (itarget_part_weights[0] - graph.part_weights[0]).abs();

    let mut moved = vec![-1 as Idx; num_vertices];
    let mut swaps = vec![0 as Idx; num_vertices];
    let mut queues = [PQueue::new(num_vertices), PQueue::new(num_vertices)];

    for _pass in 0..niter {
        let initcut = graph.edge_cut;
        let mut edge_cut = initcut;
        let mut newcut = initcut;
        let mut edge_cutorder: i32 = -1;
        let mut mindiff = (itarget_part_weights[0] - graph.part_weights[0]).abs();

        queues[0].reset();
        queues[1].reset();

        // Insert boundary vertices into PQs
        let num_boundary = graph.num_boundary as usize;
        let mut perm = vec![0 as Idx; num_boundary];
        if num_boundary < 10 {
            ctrl.rng
                .rand_array_permute(num_boundary, &mut perm, 0, true);
        } else {
            // C METIS FM_2WayCutRefine: irandArrayPermute(num_boundary, perm, num_boundary, 1)
            ctrl.rng.rand_array_permute_with_nshuffles(
                num_boundary,
                &mut perm,
                0,
                num_boundary,
                true,
            );
        }

        for ii in 0..num_boundary {
            let v = graph.boundary_list[perm[ii] as usize] as usize;
            let gain = graph.external_degree[v] - graph.internal_degree[v];
            queues[graph.partition[v] as usize].insert(v as Idx, gain as f64);
        }

        let mut nswaps = 0usize;

        loop {
            // Select source partition (overweight side)
            let from = if itarget_part_weights[0] - graph.part_weights[0]
                < itarget_part_weights[1] - graph.part_weights[1]
            {
                0
            } else {
                1
            };
            let to = 1 - from;

            let best_vertex = if let Some((v, _)) = queues[from].get_top() {
                v as usize
            } else {
                break;
            };

            newcut -= graph.external_degree[best_vertex] - graph.internal_degree[best_vertex];

            // Update partition weights
            graph.part_weights[to] += graph.vertex_weights[best_vertex];
            graph.part_weights[from] -= graph.vertex_weights[best_vertex];

            // Check if this is a new best
            let newdiff = (itarget_part_weights[0] - graph.part_weights[0]).abs();
            if (newcut < edge_cut && newdiff <= origdiff + avgvertex_weights)
                || (newcut == edge_cut && newdiff < mindiff)
            {
                edge_cut = newcut;
                mindiff = newdiff;
                edge_cutorder = nswaps as i32;
            } else if nswaps as i32 - edge_cutorder > limit as i32 {
                // Undo this move and stop
                graph.part_weights[from] += graph.vertex_weights[best_vertex];
                graph.part_weights[to] -= graph.vertex_weights[best_vertex];
                newcut += graph.external_degree[best_vertex] - graph.internal_degree[best_vertex];
                break;
            }

            // Commit the move
            graph.partition[best_vertex] = to as Idx;
            moved[best_vertex] = nswaps as Idx;
            swaps[nswaps] = best_vertex as Idx;

            // Swap id and ed for the moved vertex
            let tmp = graph.internal_degree[best_vertex];
            graph.internal_degree[best_vertex] = graph.external_degree[best_vertex];
            graph.external_degree[best_vertex] = tmp;

            // Update boundary of moved vertex
            if graph.external_degree[best_vertex] == 0
                && graph.xadj[best_vertex] < graph.xadj[best_vertex + 1]
            {
                graph.remove_from_boundary(best_vertex);
            }

            // Update neighbors
            for kk in graph.xadj[best_vertex] as usize..graph.xadj[best_vertex + 1] as usize {
                let k = graph.adjacency[kk] as usize;
                let weight_delta = if to == graph.partition[k] as usize {
                    graph.edge_weights[kk]
                } else {
                    -graph.edge_weights[kk]
                };

                graph.internal_degree[k] += weight_delta;
                graph.external_degree[k] -= weight_delta;

                // Update boundary and PQ for k
                if graph.boundary_map[k] != -1 {
                    // k is currently on boundary
                    if graph.external_degree[k] == 0 && graph.xadj[k] < graph.xadj[k + 1] {
                        graph.remove_from_boundary(k);
                        if moved[k] == -1 {
                            queues[graph.partition[k] as usize].delete(k as Idx);
                        }
                    } else if moved[k] == -1 {
                        queues[graph.partition[k] as usize].update(
                            k as Idx,
                            (graph.external_degree[k] - graph.internal_degree[k]) as f64,
                        );
                    }
                } else {
                    // k is not on boundary
                    if graph.external_degree[k] > 0 {
                        graph.add_to_boundary(k);
                        if moved[k] == -1 {
                            queues[graph.partition[k] as usize].insert(
                                k as Idx,
                                (graph.external_degree[k] - graph.internal_degree[k]) as f64,
                            );
                        }
                    }
                }
            }

            nswaps += 1;
        }

        // Rollback to best position
        // Reset moved flags
        for i in 0..nswaps {
            moved[swaps[i] as usize] = -1;
        }

        let rollback_start = if edge_cutorder < 0 {
            0
        } else {
            edge_cutorder as usize + 1
        };
        for i in (rollback_start..nswaps).rev() {
            let best_vertex = swaps[i] as usize;
            let to = graph.partition[best_vertex] as usize;
            let from = 1 - to;

            // Undo the move
            graph.partition[best_vertex] = from as Idx;

            // Swap id/ed back
            let tmp = graph.internal_degree[best_vertex];
            graph.internal_degree[best_vertex] = graph.external_degree[best_vertex];
            graph.external_degree[best_vertex] = tmp;

            graph.part_weights[from] += graph.vertex_weights[best_vertex];
            graph.part_weights[to] -= graph.vertex_weights[best_vertex];

            // Fix boundary
            if graph.external_degree[best_vertex] == 0
                && graph.xadj[best_vertex] < graph.xadj[best_vertex + 1]
            {
                graph.remove_from_boundary(best_vertex);
            } else if graph.boundary_map[best_vertex] == -1 {
                graph.add_to_boundary(best_vertex);
            }

            // Update neighbors
            for kk in graph.xadj[best_vertex] as usize..graph.xadj[best_vertex + 1] as usize {
                let k = graph.adjacency[kk] as usize;
                let weight_delta = if from == graph.partition[k] as usize {
                    graph.edge_weights[kk]
                } else {
                    -graph.edge_weights[kk]
                };

                graph.internal_degree[k] += weight_delta;
                graph.external_degree[k] -= weight_delta;

                if graph.boundary_map[k] != -1
                    && graph.external_degree[k] == 0
                    && graph.xadj[k] < graph.xadj[k + 1]
                {
                    graph.remove_from_boundary(k);
                } else if graph.boundary_map[k] == -1 && graph.external_degree[k] > 0 {
                    graph.add_to_boundary(k);
                }
            }
        }

        graph.edge_cut = edge_cut;

        if edge_cutorder <= 0 || edge_cut == initcut {
            break;
        }
    }
}

/// FM_2WayRefine dispatcher: calls FM_2WayCutRefine for ncon==1,
/// FM_Mc2WayCutRefine for ncon > 1.
pub fn fm_2way_refine(
    ctrl: &mut Control,
    graph: &mut GraphData,
    target_part_weights: &[Real],
    niter: Idx,
) {
    if graph.num_constraints == 1 {
        fm_2way_cut_refine(ctrl, graph, target_part_weights, niter);
    } else {
        fm_mc_2way_cut_refine(ctrl, graph, target_part_weights, niter);
    }
}

/// FM_Mc2WayCutRefine: multi-constraint cut-focused FM refinement.
///
/// Matches C METIS FM_Mc2WayCutRefine from fm.c exactly:
/// - Uses 2*ncon priority queues (one per constraint per partition)
/// - qnum[i] = iargmax_nrm(...) for dominant constraint per vertex
/// - SelectQueue picks which queue to drain
/// - ffactor = 0.5 / max(20, nvtxs) for balance fudge factor
/// - ubfactors relaxed to current imbalance via ComputeLoadImbalanceDiffVec
/// - Best criteria: (newcut < mincut && newbal <= ffactor) ||
///                  (newcut == mincut && (newbal < minbal ||
///                   (newbal == minbal && BetterBalance2Way)))
/// - Rollback to best position
#[allow(unused_assignments)]
pub fn fm_mc_2way_cut_refine(
    ctrl: &mut Control,
    graph: &mut GraphData,
    _target_part_weights: &[Real],
    niter: Idx,
) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    if graph.num_boundary <= 0 {
        return;
    }

    let limit = (0.01 * num_vertices as f64).max(25.0).min(150.0) as usize;
    let ffactor: Real = 0.5 / (20.0f64.max(num_vertices as f64)) as Real;

    // Initialize 2*ncon queues
    let mut queues: Vec<PQueue> = (0..2 * ncon).map(|_| PQueue::new(num_vertices)).collect();

    // Compute qnum for each vertex
    let mut qnum = vec![0usize; num_vertices];
    for i in 0..num_vertices {
        qnum[i] = iargmax_nrm(
            ncon,
            &graph.vertex_weights[i * ncon..(i + 1) * ncon],
            &graph.inv_total_vertex_weight,
        );
    }

    // Compute ubfactors: relaxed tolerance = max(current_imbalance, ctrl.imbalance_tols)
    let mut ubfactors = vec![0.0 as Real; ncon];
    let origbal = compute_load_imbalance_diff_vec(
        &graph.part_weights,
        ncon,
        2,
        &ctrl.partition_ij_balance_multipliers,
        &ctrl.imbalance_tols,
        &mut ubfactors,
    );
    for i in 0..ncon {
        ubfactors[i] = if ubfactors[i] > 0.0 {
            ctrl.imbalance_tols[i] + ubfactors[i]
        } else {
            ctrl.imbalance_tols[i]
        };
    }

    let mut moved = vec![-1 as Idx; num_vertices];
    let mut swaps = vec![0 as Idx; num_vertices];

    for _pass in 0..niter {
        for q in queues.iter_mut() {
            q.reset();
        }

        let mut mincutorder: i32 = -1;
        let mut newcut = graph.edge_cut;
        let mut mincut = graph.edge_cut;
        let initcut = graph.edge_cut;

        let mut minbalv = vec![0.0 as Real; ncon];
        let mut minbal = compute_load_imbalance_diff_vec(
            &graph.part_weights,
            ncon,
            2,
            &ctrl.partition_ij_balance_multipliers,
            &ubfactors,
            &mut minbalv,
        );

        // Insert boundary vertices into PQs
        let num_boundary = graph.num_boundary as usize;
        let mut perm = vec![0 as Idx; num_boundary];
        if num_boundary < 10 {
            ctrl.rng
                .rand_array_permute(num_boundary, &mut perm, 0, true);
        } else {
            ctrl.rng.rand_array_permute_with_nshuffles(
                num_boundary,
                &mut perm,
                0,
                num_boundary / 5,
                true,
            );
        }

        for ii in 0..num_boundary {
            let i = graph.boundary_list[perm[ii] as usize] as usize;
            let rgain = graph.external_degree[i] - graph.internal_degree[i];
            queues[2 * qnum[i] + graph.partition[i] as usize].insert(i as Idx, rgain as f64);
        }

        let mut nswaps = 0usize;
        let mut newbalv = vec![0.0 as Real; ncon];

        loop {
            if nswaps >= num_vertices {
                break;
            }

            let (from, cnum) = select_queue(
                &graph.part_weights,
                ncon,
                &ctrl.partition_ij_balance_multipliers,
                &ubfactors,
                &queues,
            );
            let to = if from >= 0 { (from + 1) % 2 } else { -1 };

            if from == -1 {
                break;
            }
            let higain = if let Some((v, _)) = queues[2 * cnum as usize + from as usize].get_top() {
                v as usize
            } else {
                break;
            };

            newcut -= graph.external_degree[higain] - graph.internal_degree[higain];

            // Update part_weights
            for j in 0..ncon {
                graph.part_weights[to as usize * ncon + j] +=
                    graph.vertex_weights[higain * ncon + j];
                graph.part_weights[from as usize * ncon + j] -=
                    graph.vertex_weights[higain * ncon + j];
            }

            let newbal = compute_load_imbalance_diff_vec(
                &graph.part_weights,
                ncon,
                2,
                &ctrl.partition_ij_balance_multipliers,
                &ubfactors,
                &mut newbalv,
            );

            if (newcut < mincut && newbal <= ffactor)
                || (newcut == mincut
                    && (newbal < minbal
                        || (newbal == minbal && better_balance_2way(ncon, &minbalv, &newbalv))))
            {
                mincut = newcut;
                minbal = newbal;
                mincutorder = nswaps as i32;
                minbalv.copy_from_slice(&newbalv);
            } else if nswaps as i32 - mincutorder > limit as i32 {
                // Undo last move
                newcut += graph.external_degree[higain] - graph.internal_degree[higain];
                for j in 0..ncon {
                    graph.part_weights[from as usize * ncon + j] +=
                        graph.vertex_weights[higain * ncon + j];
                    graph.part_weights[to as usize * ncon + j] -=
                        graph.vertex_weights[higain * ncon + j];
                }
                break;
            }

            graph.partition[higain] = to as Idx;
            moved[higain] = nswaps as Idx;
            swaps[nswaps] = higain as Idx;

            // Swap id and ed
            let tmp = graph.internal_degree[higain];
            graph.internal_degree[higain] = graph.external_degree[higain];
            graph.external_degree[higain] = tmp;

            if graph.external_degree[higain] == 0 && graph.xadj[higain] < graph.xadj[higain + 1] {
                graph.remove_from_boundary(higain);
            }

            // Update neighbors
            for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
                let k = graph.adjacency[kk] as usize;
                let weight_delta = if to == graph.partition[k] as i32 {
                    graph.edge_weights[kk]
                } else {
                    -graph.edge_weights[kk]
                };

                graph.internal_degree[k] += weight_delta;
                graph.external_degree[k] -= weight_delta;

                if graph.boundary_map[k] != -1 {
                    if graph.external_degree[k] == 0 {
                        graph.remove_from_boundary(k);
                        if moved[k] == -1 {
                            queues[2 * qnum[k] + graph.partition[k] as usize].delete(k as Idx);
                        }
                    } else if moved[k] == -1 {
                        let rgain = graph.external_degree[k] - graph.internal_degree[k];
                        queues[2 * qnum[k] + graph.partition[k] as usize]
                            .update(k as Idx, rgain as f64);
                    }
                } else {
                    if graph.external_degree[k] > 0 {
                        graph.add_to_boundary(k);
                        if moved[k] == -1 {
                            let rgain = graph.external_degree[k] - graph.internal_degree[k];
                            queues[2 * qnum[k] + graph.partition[k] as usize]
                                .insert(k as Idx, rgain as f64);
                        }
                    }
                }
            }

            nswaps += 1;
        }

        // Rollback to best position
        for i in 0..nswaps {
            moved[swaps[i] as usize] = -1;
        }

        if nswaps > 0 {
            let mut i = nswaps as i32 - 1;
            while i > mincutorder {
                let higain = swaps[i as usize] as usize;
                let to = (graph.partition[higain] + 1) % 2;
                graph.partition[higain] = to;

                let tmp = graph.internal_degree[higain];
                graph.internal_degree[higain] = graph.external_degree[higain];
                graph.external_degree[higain] = tmp;

                if graph.external_degree[higain] == 0
                    && graph.boundary_map[higain] != -1
                    && graph.xadj[higain] < graph.xadj[higain + 1]
                {
                    graph.remove_from_boundary(higain);
                } else if graph.external_degree[higain] > 0 && graph.boundary_map[higain] == -1 {
                    graph.add_to_boundary(higain);
                }

                for j in 0..ncon {
                    graph.part_weights[to as usize * ncon + j] +=
                        graph.vertex_weights[higain * ncon + j];
                    graph.part_weights[((to as usize + 1) % 2) * ncon + j] -=
                        graph.vertex_weights[higain * ncon + j];
                }

                for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
                    let k = graph.adjacency[kk] as usize;
                    let weight_delta = if to == graph.partition[k] {
                        graph.edge_weights[kk]
                    } else {
                        -graph.edge_weights[kk]
                    };

                    graph.internal_degree[k] += weight_delta;
                    graph.external_degree[k] -= weight_delta;

                    if graph.boundary_map[k] != -1 && graph.external_degree[k] == 0 {
                        graph.remove_from_boundary(k);
                    }
                    if graph.boundary_map[k] == -1 && graph.external_degree[k] > 0 {
                        graph.add_to_boundary(k);
                    }
                }

                i -= 1;
            }
        }

        graph.edge_cut = mincut;

        if mincutorder <= 0 || mincut == initcut {
            break;
        }
    }

    let _ = origbal; // suppress unused warning
}
