use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::util::pqueue::PQueue;

/// Compute the load imbalance diff for a 2-way partition.
///
/// Returns max over i in 0..nparts, j in 0..ncon of:
///   part_weights[i*ncon+j] * partition_ij_balance_multipliers[i*ncon+j] - imbalance_tols[j]
///
/// If <= 0, the partition is balanced.
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
            let diff =
                graph.part_weights[i * ncon + j] as Real * partition_ij_balance_multipliers[i * ncon + j] - imbalance_tols[j];
            if diff > max_diff {
                max_diff = diff;
            }
        }
    }
    max_diff
}

/// Balance a 2-way partition.
///
/// Matches METIS Balance2Way from balance.c:
/// 1. Check ComputeLoadImbalanceDiff <= 0 (already balanced) -- return early.
/// 2. For ncon == 1: check if the weight difference is tiny (< 3*total_vertex_weight/num_vertices) -- return early.
/// 3. If num_boundary > 0, use Bnd2WayBalance (boundary-only PQ balancing).
/// 4. If num_boundary == 0, use General2WayBalance (all-vertex PQ balancing).
/// 5. For ncon > 1, use McGeneral2WayBalance (multi-constraint).
pub fn balance_2way(ctrl: &mut Control, graph: &mut GraphData, ntarget_part_weights: &[Real]) {
    // Step 1: already balanced?
    let load_imb = compute_load_imbalance_diff(graph, 2, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols);
    if load_imb <= 0.0 {
        return;
    }

    let ncon = graph.num_constraints as usize;

    if ncon == 1 {
        // Step 2: check if the difference is negligibly small
        // C METIS uses integer arithmetic: 3*tvwgt[0]/nvtxs (integer division)
        let tpwgt0 = (ntarget_part_weights[0] * graph.total_vertex_weight[0] as Real) as Real;
        let diff = (tpwgt0 - graph.part_weights[0] as Real).abs();
        let threshold = (3 * graph.total_vertex_weight[0] / graph.num_vertices) as Real;
        if diff < threshold {
            return;
        }

        // Step 3/4: boundary vs general balancing
        if graph.num_boundary > 0 {
            bnd_2way_balance(ctrl, graph, ntarget_part_weights);
        } else {
            general_2way_balance(ctrl, graph, ntarget_part_weights);
        }
    } else {
        mc_general_2way_balance(ctrl, graph, ntarget_part_weights);
    }
}

/// Bnd2WayBalance: greedy boundary-vertex balancing for single-constraint graphs.
///
/// Matches METIS Bnd2WayBalance from balance.c:
/// - Computes integer target weights: target_part_weights[0] = total_vertex_weight*ntarget_part_weights[0], target_part_weights[1] = total_vertex_weight - target_part_weights[0]
/// - Determines overweight (from) and underweight (to) partitions
/// - mindiff = |target_part_weights[0] - part_weights[0]|
/// - Permutes boundary indices, inserts eligible vertices into PQ (gain = ed - id)
/// - Greedily moves vertices: pop max-gain, move if part_weights[to]+vertex_weights <= target_part_weights[to], else break
/// - Incrementally updates where, part_weights, id/ed, boundary, edge_cut, and PQ
/// - NO rollback -- this is a greedy monotone pass
fn bnd_2way_balance(ctrl: &mut Control, graph: &mut GraphData, ntarget_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;
    let num_boundary = graph.num_boundary as usize;

    // Integer target weights
    let target_part_weights0 = (graph.total_vertex_weight[0] as Real * ntarget_part_weights[0]) as Idx;
    let target_part_weights1 = graph.total_vertex_weight[0] - target_part_weights0;
    let target_part_weights = [target_part_weights0, target_part_weights1];

    let mindiff = (target_part_weights[0] - graph.part_weights[0]).abs();

    // Determine from (overweight) and to (underweight)
    let from: usize = if graph.part_weights[0] < target_part_weights[0] { 1 } else { 0 };
    let to: usize = 1 - from;

    let mut queue = PQueue::new(num_vertices);
    let mut moved = vec![-1 as Idx; num_vertices];

    // Permute boundary indices: irandArrayPermute(num_boundary, perm, num_boundary/5, 1)
    let mut perm = vec![0 as Idx; num_boundary];
    if num_boundary < 10 {
        ctrl.rng.rand_array_permute(num_boundary, &mut perm, 0, true);
    } else {
        ctrl.rng.rand_array_permute_with_nshuffles(num_boundary, &mut perm, 0, num_boundary / 5, true);
    }

    // Insert boundary vertices of 'from' partition with vertex_weights <= mindiff
    for ii in 0..num_boundary {
        let i = perm[ii] as usize;
        let v = graph.boundary_list[i] as usize;
        if graph.partition[v] as usize == from && graph.vertex_weights[v] <= mindiff {
            queue.insert(v as Idx, (graph.external_degree[v] - graph.internal_degree[v]) as f64);
        }
    }

    let mut edge_cut = graph.edge_cut;
    let mut nswaps: Idx = 0;

    // Greedy move loop
    loop {
        let best_vertex = if let Some((v, _)) = queue.get_top() {
            v as usize
        } else {
            break;
        };

        // If moving this vertex would overload the 'to' partition, stop
        if graph.part_weights[to] + graph.vertex_weights[best_vertex] > target_part_weights[to] {
            break;
        }

        edge_cut -= graph.external_degree[best_vertex] - graph.internal_degree[best_vertex];
        graph.part_weights[to] += graph.vertex_weights[best_vertex];
        graph.part_weights[from] -= graph.vertex_weights[best_vertex];

        graph.partition[best_vertex] = to as Idx;
        moved[best_vertex] = nswaps;

        // Swap id and ed for the moved vertex
        let tmp = graph.internal_degree[best_vertex];
        graph.internal_degree[best_vertex] = graph.external_degree[best_vertex];
        graph.external_degree[best_vertex] = tmp;

        if graph.external_degree[best_vertex] == 0 && graph.xadj[best_vertex] < graph.xadj[best_vertex + 1] {
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

            if graph.boundary_map[k] != -1 {
                // k is on boundary
                if graph.external_degree[k] == 0 {
                    graph.remove_from_boundary(k);
                    if moved[k] == -1 && graph.partition[k] as usize == from && graph.vertex_weights[k] <= mindiff {
                        queue.delete(k as Idx);
                    }
                } else {
                    if moved[k] == -1 && graph.partition[k] as usize == from && graph.vertex_weights[k] <= mindiff {
                        queue.update(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
                    }
                }
            } else {
                // k is not on boundary
                if graph.external_degree[k] > 0 {
                    graph.add_to_boundary(k);
                    if moved[k] == -1 && graph.partition[k] as usize == from && graph.vertex_weights[k] <= mindiff {
                        queue.insert(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
                    }
                }
            }
        }

        nswaps += 1;
    }

    graph.edge_cut = edge_cut;
}

/// General2WayBalance: greedy all-vertex balancing for single-constraint graphs when num_boundary == 0.
///
/// Matches METIS General2WayBalance from balance.c:
/// - Inserts ALL vertices from 'from' partition with vertex_weights <= mindiff into PQ
/// - Greedily moves vertices; neighbor loop only does rpqUpdate (no insert/delete)
/// - Boundary updates are separate from PQ updates
/// - Used when there are no boundary vertices
fn general_2way_balance(ctrl: &mut Control, graph: &mut GraphData, ntarget_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;

    // Integer target weights
    let target_part_weights0 = (graph.total_vertex_weight[0] as Real * ntarget_part_weights[0]) as Idx;
    let target_part_weights1 = graph.total_vertex_weight[0] - target_part_weights0;
    let target_part_weights = [target_part_weights0, target_part_weights1];

    let mindiff = (target_part_weights[0] - graph.part_weights[0]).abs();

    // Determine from (overweight) and to (underweight)
    let from: usize = if graph.part_weights[0] < target_part_weights[0] { 1 } else { 0 };
    let to: usize = 1 - from;

    let mut queue = PQueue::new(num_vertices);
    let mut moved = vec![-1 as Idx; num_vertices];

    // Permute vertex indices: irandArrayPermute(num_vertices, perm, num_vertices/5, 1)
    let mut perm = vec![0 as Idx; num_vertices];
    if num_vertices < 10 {
        ctrl.rng.rand_array_permute(num_vertices, &mut perm, 0, true);
    } else {
        ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut perm, 0, num_vertices / 5, true);
    }

    // Insert ALL vertices from 'from' partition with small weight
    for ii in 0..num_vertices {
        let i = perm[ii] as usize;
        if graph.partition[i] as usize == from && graph.vertex_weights[i] <= mindiff {
            queue.insert(i as Idx, (graph.external_degree[i] - graph.internal_degree[i]) as f64);
        }
    }

    let mut edge_cut = graph.edge_cut;
    let mut nswaps: Idx = 0;

    // Greedy move loop
    loop {
        let best_vertex = if let Some((v, _)) = queue.get_top() {
            v as usize
        } else {
            break;
        };

        // If moving this vertex would overload the 'to' partition, stop
        if graph.part_weights[to] + graph.vertex_weights[best_vertex] > target_part_weights[to] {
            break;
        }

        edge_cut -= graph.external_degree[best_vertex] - graph.internal_degree[best_vertex];
        graph.part_weights[to] += graph.vertex_weights[best_vertex];
        graph.part_weights[from] -= graph.vertex_weights[best_vertex];

        graph.partition[best_vertex] = to as Idx;
        moved[best_vertex] = nswaps;

        // Swap id and ed for the moved vertex
        let tmp = graph.internal_degree[best_vertex];
        graph.internal_degree[best_vertex] = graph.external_degree[best_vertex];
        graph.external_degree[best_vertex] = tmp;

        // Update boundary status of moved vertex (two separate checks, matching C METIS)
        if graph.external_degree[best_vertex] == 0 && graph.boundary_map[best_vertex] != -1 && graph.xadj[best_vertex] < graph.xadj[best_vertex + 1] {
            graph.remove_from_boundary(best_vertex);
        }
        if graph.external_degree[best_vertex] > 0 && graph.boundary_map[best_vertex] == -1 {
            graph.add_to_boundary(best_vertex);
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

            // Update queue position (only rpqUpdate, no insert/delete)
            if moved[k] == -1 && graph.partition[k] as usize == from && graph.vertex_weights[k] <= mindiff {
                queue.update(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
            }

            // Update boundary information (separate from PQ)
            if graph.external_degree[k] == 0 && graph.boundary_map[k] != -1 {
                graph.remove_from_boundary(k);
            } else if graph.external_degree[k] > 0 && graph.boundary_map[k] == -1 {
                graph.add_to_boundary(k);
            }
        }

        nswaps += 1;
    }

    graph.edge_cut = edge_cut;
}

/// McGeneral2WayBalance: greedy multi-constraint 2-way balancing.
///
/// For ncon > 1: inserts all vertices from the overweight partition that could help
/// reduce the imbalance, greedily moves them while the target side can absorb them.
fn mc_general_2way_balance(_ctrl: &mut Control, graph: &mut GraphData, ntarget_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Integer target weights per constraint per partition
    let mut target_part_weights = vec![0 as Idx; 2 * ncon];
    for j in 0..ncon {
        target_part_weights[j] = (graph.total_vertex_weight[j] as Real * ntarget_part_weights[j]) as Idx;
        target_part_weights[ncon + j] = graph.total_vertex_weight[j] - target_part_weights[j];
    }

    // Determine which partition is overweight (using first constraint that is out of balance)
    let from: usize = if graph.part_weights[0] > target_part_weights[0] { 0 } else { 1 };
    let to: usize = 1 - from;

    // Build PQ from all vertices in the 'from' partition
    let mut queue = PQueue::new(num_vertices);

    for v in 0..num_vertices {
        if graph.partition[v] as usize == from {
            let gain = graph.external_degree[v] - graph.internal_degree[v];
            queue.insert(v as Idx, gain as f64);
        }
    }

    // Greedy move loop
    while let Some((higain_idx, _key)) = queue.get_top() {
        let best_vertex = higain_idx as usize;

        // Check if moving this vertex would overload the 'to' partition in any constraint
        let mut fits = true;
        for j in 0..ncon {
            let new_to_wgt = graph.part_weights[to * ncon + j] + graph.vertex_weights[best_vertex * ncon + j];
            if new_to_wgt > target_part_weights[to * ncon + j] {
                fits = false;
                break;
            }
        }
        if !fits {
            break;
        }

        // Move vertex from 'from' to 'to'
        graph.edge_cut -= graph.external_degree[best_vertex] - graph.internal_degree[best_vertex];

        graph.partition[best_vertex] = to as Idx;
        for j in 0..ncon {
            graph.part_weights[to * ncon + j] += graph.vertex_weights[best_vertex * ncon + j];
            graph.part_weights[from * ncon + j] -= graph.vertex_weights[best_vertex * ncon + j];
        }

        // Swap id and ed for the moved vertex
        let tmp = graph.internal_degree[best_vertex];
        graph.internal_degree[best_vertex] = graph.external_degree[best_vertex];
        graph.external_degree[best_vertex] = tmp;

        // Update boundary status of moved vertex
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
            let weight_delta = if to == graph.partition[k] as usize {
                graph.edge_weights[kk]
            } else {
                -graph.edge_weights[kk]
            };

            graph.internal_degree[k] += weight_delta;
            graph.external_degree[k] -= weight_delta;

            // Update boundary and PQ for neighbor k
            if graph.boundary_map[k] != -1 {
                // k is on boundary
                if graph.external_degree[k] == 0 && graph.xadj[k] < graph.xadj[k + 1] {
                    graph.remove_from_boundary(k);
                    if queue.contains(k as Idx) {
                        queue.delete(k as Idx);
                    }
                } else if queue.contains(k as Idx) {
                    queue.update(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
                }
            } else {
                // k is not on boundary
                if graph.external_degree[k] > 0 {
                    graph.add_to_boundary(k);
                    if graph.partition[k] as usize == from {
                        queue.insert(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
                    }
                }
            }
        }
    }
}
