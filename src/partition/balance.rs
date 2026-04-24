use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::util::pqueue::PQueue;

/// iargmax_nrm: returns the index of the maximum x[i]*y[i].
/// Used to find the dominant constraint for each vertex.
pub(crate) fn iargmax_nrm(n: usize, x: &[Idx], y: &[Real]) -> usize {
    let mut max = 0;
    for i in 1..n {
        if x[i] as Real * y[i] > x[max] as Real * y[max] {
            max = i;
        }
    }
    max
}

/// iargmax2_nrm: returns the index of the second-largest x[i]*y[i].
/// Used in McGeneral2WayBalance for empty queue redistribution.
pub(crate) fn iargmax2_nrm(n: usize, x: &[Idx], y: &[Real]) -> usize {
    let (mut max1, mut max2);
    if x[0] as Real * y[0] > x[1] as Real * y[1] {
        max1 = 0;
        max2 = 1;
    } else {
        max1 = 1;
        max2 = 0;
    }
    for i in 2..n {
        if x[i] as Real * y[i] > x[max1] as Real * y[max1] {
            max2 = max1;
            max1 = i;
        } else if x[i] as Real * y[i] > x[max2] as Real * y[max2] {
            max2 = i;
        }
    }
    max2
}

/// ComputeLoadImbalanceDiffVec: for each constraint, compute the max imbalance diff
/// across all partitions, and store it in diffvec[i].
/// Returns the overall max.
pub(crate) fn compute_load_imbalance_diff_vec(
    part_weights: &[Idx],
    ncon: usize,
    nparts: usize,
    pijbm: &[Real],
    ubfactors: &[Real],
    diffvec: &mut [Real],
) -> Real {
    let mut max = -1.0 as Real;
    for i in 0..ncon {
        diffvec[i] = part_weights[i] as Real * pijbm[i] - ubfactors[i];
        for j in 1..nparts {
            let cur = part_weights[j * ncon + i] as Real * pijbm[j * ncon + i] - ubfactors[i];
            if cur > diffvec[i] {
                diffvec[i] = cur;
            }
        }
        if max < diffvec[i] {
            max = diffvec[i];
        }
    }
    max
}

/// BetterBalance2Way: returns true if y is better balanced than x.
/// Both vectors are ubfactor-centered load imbalance vectors.
pub(crate) fn better_balance_2way(n: usize, x: &[Real], y: &[Real]) -> bool {
    let mut nrm1: Real = 0.0;
    let mut nrm2: Real = 0.0;
    for i in (0..n).rev() {
        if x[i] > 0.0 {
            nrm1 += x[i] * x[i];
        }
        if y[i] > 0.0 {
            nrm2 += y[i] * y[i];
        }
    }
    nrm2 < nrm1
}

/// SelectQueue: select the partition and constraint queue to move vertices from.
/// Returns (from_partition, constraint_num). Returns (-1, -1) if nothing to move.
pub(crate) fn select_queue(
    part_weights: &[Idx],
    ncon: usize,
    pijbm: &[Real],
    ubfactors: &[Real],
    queues: &[PQueue],
) -> (i32, i32) {
    let mut from: i32 = -1;
    let mut cnum: i32 = -1;
    let mut max: Real = 0.0;

    // First: find the most violated balancing constraint
    for part in 0..2i32 {
        for i in 0..ncon {
            let idx = part as usize * ncon + i;
            let tmp = part_weights[idx] as Real * pijbm[idx] - ubfactors[i];
            // the '=' ensures that under tight constraints the partition at max is selected
            if tmp >= max {
                max = tmp;
                from = part;
                cnum = i as i32;
            }
        }
    }

    if from != -1 {
        // If the desired queue is empty, select a queue from the same side
        let qi = 2 * cnum as usize + from as usize;
        if queues[qi].len() == 0 {
            let mut found = false;
            for i in 0..ncon {
                if queues[2 * i + from as usize].len() > 0 {
                    max = part_weights[from as usize * ncon + i] as Real * pijbm[from as usize * ncon + i] - ubfactors[i];
                    cnum = i as i32;
                    found = true;
                    break;
                }
            }

            if found {
                let start = cnum as usize + 1;
                for i in start..ncon {
                    let tmp = part_weights[from as usize * ncon + i] as Real * pijbm[from as usize * ncon + i] - ubfactors[i];
                    if tmp > max && queues[2 * i + from as usize].len() > 0 {
                        max = tmp;
                        cnum = i as i32;
                    }
                }
            }
        }
    } else {
        // Partition does not violate balancing constraints: select based on cut criteria
        let mut fmax: f64 = f64::NEG_INFINITY;
        for part in 0..2i32 {
            for i in 0..ncon {
                let qi = 2 * i + part as usize;
                if queues[qi].len() > 0 {
                    if let Some((_, key)) = queues[qi].peek_top() {
                        if from == -1 || key > fmax {
                            fmax = key;
                            from = part;
                            cnum = i as i32;
                        }
                    }
                }
            }
        }
    }

    (from, cnum)
}

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

/// McGeneral2WayBalance: multi-constraint 2-way balancing with 2*ncon queues.
///
/// Matches C METIS McGeneral2WayBalance from balance.c exactly:
/// - Uses 2*ncon priority queues (one per constraint per partition)
/// - qnum[i] = iargmax_nrm(...) for dominant constraint per vertex
/// - Redistributes vertices from empty queues using iargmax2_nrm
/// - SelectQueue picks which queue to drain
/// - Greedy moves with rollback to best position
/// - Stops when balance target is achieved (minbal <= 0)
fn mc_general_2way_balance(ctrl: &mut Control, graph: &mut GraphData, _ntarget_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    let limit = (0.01 * num_vertices as f64).max(15.0).min(100.0) as usize;

    // Initialize 2*ncon queues
    let mut queues: Vec<PQueue> = (0..2 * ncon).map(|_| PQueue::new(num_vertices)).collect();
    let mut qsizes = vec![0usize; 2 * ncon];

    // Compute qnum for each vertex
    let mut qnum = vec![0usize; num_vertices];
    for i in 0..num_vertices {
        qnum[i] = iargmax_nrm(ncon, &graph.vertex_weights[i * ncon..(i + 1) * ncon], &graph.inv_total_vertex_weight);
        qsizes[2 * qnum[i] + graph.partition[i] as usize] += 1;
    }

    // Redistribute vertices from empty queues
    for from in 0..2usize {
        for j in 0..ncon {
            if qsizes[2 * j + from] == 0 {
                for i in 0..num_vertices {
                    if graph.partition[i] as usize != from {
                        continue;
                    }
                    let k = iargmax2_nrm(ncon, &graph.vertex_weights[i * ncon..(i + 1) * ncon], &graph.inv_total_vertex_weight);
                    if k == j
                        && qsizes[2 * qnum[i] + from] > qsizes[2 * j + from]
                        && graph.vertex_weights[i * ncon + qnum[i]] as Real * graph.inv_total_vertex_weight[qnum[i]]
                            < 1.3 * graph.vertex_weights[i * ncon + j] as Real * graph.inv_total_vertex_weight[j]
                    {
                        qsizes[2 * qnum[i] + from] -= 1;
                        qsizes[2 * j + from] += 1;
                        qnum[i] = j;
                    }
                }
            }
        }
    }

    let mut minbalv = vec![0.0 as Real; ncon];
    let mut minbal = compute_load_imbalance_diff_vec(
        &graph.part_weights, ncon, 2, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols, &mut minbalv,
    );

    let mut newcut = graph.edge_cut;
    let mut mincut = graph.edge_cut;
    let mut mincutorder: i32 = -1;

    let mut moved = vec![-1 as Idx; num_vertices];
    let mut swaps = vec![0 as Idx; num_vertices];

    // Insert ALL vertices into priority queues
    let mut perm = vec![0 as Idx; num_vertices];
    if num_vertices < 10 {
        ctrl.rng.rand_array_permute(num_vertices, &mut perm, 0, true);
    } else {
        ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut perm, 0, num_vertices / 10, true);
    }

    for ii in 0..num_vertices {
        let i = perm[ii] as usize;
        let gain = graph.external_degree[i] - graph.internal_degree[i];
        queues[2 * qnum[i] + graph.partition[i] as usize].insert(i as Idx, gain as f64);
    }

    let mut nswaps = 0usize;
    let mut newbalv = vec![0.0 as Real; ncon];

    loop {
        if nswaps >= num_vertices {
            break;
        }
        if minbal <= 0.0 {
            break;
        }

        let (from, cnum) = select_queue(
            &graph.part_weights, ncon, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols, &queues,
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

        // Update part_weights: iaxpy
        for j in 0..ncon {
            graph.part_weights[to as usize * ncon + j] += graph.vertex_weights[higain * ncon + j];
            graph.part_weights[from as usize * ncon + j] -= graph.vertex_weights[higain * ncon + j];
        }

        let newbal = compute_load_imbalance_diff_vec(
            &graph.part_weights, ncon, 2, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols, &mut newbalv,
        );

        if newbal < minbal
            || (newbal == minbal
                && (newcut < mincut
                    || (newcut == mincut && better_balance_2way(ncon, &minbalv, &newbalv))))
        {
            mincut = newcut;
            minbal = newbal;
            mincutorder = nswaps as i32;
            minbalv.copy_from_slice(&newbalv);
        } else if nswaps as i32 - mincutorder > limit as i32 {
            // Undo last move
            newcut += graph.external_degree[higain] - graph.internal_degree[higain];
            for j in 0..ncon {
                graph.part_weights[from as usize * ncon + j] += graph.vertex_weights[higain * ncon + j];
                graph.part_weights[to as usize * ncon + j] -= graph.vertex_weights[higain * ncon + j];
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

        // Update boundary
        if graph.external_degree[higain] == 0 && graph.boundary_map[higain] != -1 && graph.xadj[higain] < graph.xadj[higain + 1] {
            graph.remove_from_boundary(higain);
        }
        if graph.external_degree[higain] > 0 && graph.boundary_map[higain] == -1 {
            graph.add_to_boundary(higain);
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

            // Update queue position
            if moved[k] == -1 {
                queues[2 * qnum[k] + graph.partition[k] as usize].update(k as Idx, (graph.external_degree[k] - graph.internal_degree[k]) as f64);
            }

            // Update boundary
            if graph.external_degree[k] == 0 && graph.boundary_map[k] != -1 {
                graph.remove_from_boundary(k);
            } else if graph.external_degree[k] > 0 && graph.boundary_map[k] == -1 {
                graph.add_to_boundary(k);
            }
        }

        nswaps += 1;
    }

    // Rollback to best position
    let rollback_start = mincutorder + 1;
    if nswaps > 0 {
        let mut i = nswaps as i32 - 1;
        while i > mincutorder {
            let higain = swaps[i as usize] as usize;
            let to = (graph.partition[higain] + 1) % 2;
            graph.partition[higain] = to;

            let tmp = graph.internal_degree[higain];
            graph.internal_degree[higain] = graph.external_degree[higain];
            graph.external_degree[higain] = tmp;

            if graph.external_degree[higain] == 0 && graph.boundary_map[higain] != -1 && graph.xadj[higain] < graph.xadj[higain + 1] {
                graph.remove_from_boundary(higain);
            } else if graph.external_degree[higain] > 0 && graph.boundary_map[higain] == -1 {
                graph.add_to_boundary(higain);
            }

            for j in 0..ncon {
                graph.part_weights[to as usize * ncon + j] += graph.vertex_weights[higain * ncon + j];
                graph.part_weights[((to as usize + 1) % 2) * ncon + j] -= graph.vertex_weights[higain * ncon + j];
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
    let _ = rollback_start; // suppress unused warning
}
