use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::util::pqueue::PQueue;

/// Compute the load imbalance diff for a 2-way partition.
///
/// Returns max over i in 0..nparts, j in 0..ncon of:
///   pwgts[i*ncon+j] * pijbm[i*ncon+j] - ubfactors[j]
///
/// If <= 0, the partition is balanced.
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
            let diff =
                graph.pwgts[i * ncon + j] as Real * pijbm[i * ncon + j] - ubfactors[j];
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
/// 2. For ncon == 1: check if the weight difference is tiny (< 3*tvwgt/nvtxs) -- return early.
/// 3. If nbnd > 0, use Bnd2WayBalance (boundary-only PQ balancing).
/// 4. If nbnd == 0, use General2WayBalance (all-vertex PQ balancing).
/// 5. For ncon > 1, use McGeneral2WayBalance (multi-constraint).
pub fn balance_2way(ctrl: &mut Ctrl, graph: &mut GraphData, ntpwgts: &[Real]) {
    // Step 1: already balanced?
    if compute_load_imbalance_diff(graph, 2, &ctrl.pijbm, &ctrl.ubfactors) <= 0.0 {
        return;
    }

    let ncon = graph.ncon as usize;

    if ncon == 1 {
        // Step 2: check if the difference is negligibly small
        let tpwgt0 = (ntpwgts[0] * graph.tvwgt[0] as Real) as Real;
        let diff = (tpwgt0 - graph.pwgts[0] as Real).abs();
        let threshold = 3.0 * graph.tvwgt[0] as Real / graph.nvtxs as Real;
        if diff < threshold {
            return;
        }

        // Step 3/4: boundary vs general balancing
        if graph.nbnd > 0 {
            bnd_2way_balance(ctrl, graph, ntpwgts);
        } else {
            general_2way_balance(ctrl, graph, ntpwgts);
        }
    } else {
        mc_general_2way_balance(ctrl, graph, ntpwgts);
    }
}

/// Bnd2WayBalance: greedy boundary-vertex balancing for single-constraint graphs.
///
/// Matches METIS Bnd2WayBalance from balance.c:
/// - Computes integer target weights: tpwgts[0] = tvwgt*ntpwgts[0], tpwgts[1] = tvwgt - tpwgts[0]
/// - Determines overweight (from) and underweight (to) partitions
/// - mindiff = |tpwgts[0] - pwgts[0]|
/// - Permutes boundary indices, inserts eligible vertices into PQ (gain = ed - id)
/// - Greedily moves vertices: pop max-gain, move if pwgts[to]+vwgt <= tpwgts[to], else break
/// - Incrementally updates where, pwgts, id/ed, boundary, mincut, and PQ
/// - NO rollback -- this is a greedy monotone pass
fn bnd_2way_balance(ctrl: &mut Ctrl, graph: &mut GraphData, ntpwgts: &[Real]) {
    let nvtxs = graph.nvtxs as usize;
    let nbnd = graph.nbnd as usize;

    // Integer target weights
    let tpwgts0 = (graph.tvwgt[0] as Real * ntpwgts[0]) as Idx;
    let tpwgts1 = graph.tvwgt[0] - tpwgts0;
    let tpwgts = [tpwgts0, tpwgts1];

    let mindiff = (tpwgts[0] - graph.pwgts[0]).abs();

    // Determine from (overweight) and to (underweight)
    let from: usize = if graph.pwgts[0] < tpwgts[0] { 1 } else { 0 };
    let to: usize = 1 - from;

    let mut queue = PQueue::new(nvtxs);
    let mut moved = vec![-1 as Idx; nvtxs];

    // Permute boundary indices: irandArrayPermute(nbnd, perm, nbnd/5, 1)
    let mut perm = vec![0 as Idx; nbnd];
    if nbnd < 10 {
        ctrl.rng.rand_array_permute(nbnd, &mut perm, 0, true);
    } else {
        ctrl.rng.rand_array_permute_with_nshuffles(nbnd, &mut perm, 0, nbnd / 5, true);
    }

    // Insert boundary vertices of 'from' partition with vwgt <= mindiff
    for ii in 0..nbnd {
        let i = perm[ii] as usize;
        let v = graph.bndind[i] as usize;
        if graph.where_[v] as usize == from && graph.vwgt[v] <= mindiff {
            queue.insert(v as Idx, (graph.ed[v] - graph.id[v]) as f64);
        }
    }

    let mut mincut = graph.mincut;
    let mut nswaps: Idx = 0;

    // Greedy move loop
    loop {
        let higain = if let Some((v, _)) = queue.get_top() {
            v as usize
        } else {
            break;
        };

        // If moving this vertex would overload the 'to' partition, stop
        if graph.pwgts[to] + graph.vwgt[higain] > tpwgts[to] {
            break;
        }

        mincut -= graph.ed[higain] - graph.id[higain];
        graph.pwgts[to] += graph.vwgt[higain];
        graph.pwgts[from] -= graph.vwgt[higain];

        graph.where_[higain] = to as Idx;
        moved[higain] = nswaps;

        // Swap id and ed for the moved vertex
        let tmp = graph.id[higain];
        graph.id[higain] = graph.ed[higain];
        graph.ed[higain] = tmp;

        if graph.ed[higain] == 0 && graph.xadj[higain] < graph.xadj[higain + 1] {
            graph.bnd_delete(higain);
        }

        // Update neighbors
        for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
            let k = graph.adjncy[kk] as usize;
            let kwgt = if to == graph.where_[k] as usize {
                graph.adjwgt[kk]
            } else {
                -graph.adjwgt[kk]
            };

            graph.id[k] += kwgt;
            graph.ed[k] -= kwgt;

            if graph.bndptr[k] != -1 {
                // k is on boundary
                if graph.ed[k] == 0 {
                    graph.bnd_delete(k);
                    if moved[k] == -1 && graph.where_[k] as usize == from && graph.vwgt[k] <= mindiff {
                        queue.delete(k as Idx);
                    }
                } else {
                    if moved[k] == -1 && graph.where_[k] as usize == from && graph.vwgt[k] <= mindiff {
                        queue.update(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
                    }
                }
            } else {
                // k is not on boundary
                if graph.ed[k] > 0 {
                    graph.bnd_insert(k);
                    if moved[k] == -1 && graph.where_[k] as usize == from && graph.vwgt[k] <= mindiff {
                        queue.insert(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
                    }
                }
            }
        }

        nswaps += 1;
    }

    graph.mincut = mincut;
}

/// General2WayBalance: greedy all-vertex balancing for single-constraint graphs when nbnd == 0.
///
/// Matches METIS General2WayBalance from balance.c:
/// - Inserts ALL vertices from 'from' partition with vwgt <= mindiff into PQ
/// - Greedily moves vertices; neighbor loop only does rpqUpdate (no insert/delete)
/// - Boundary updates are separate from PQ updates
/// - Used when there are no boundary vertices
fn general_2way_balance(ctrl: &mut Ctrl, graph: &mut GraphData, ntpwgts: &[Real]) {
    let nvtxs = graph.nvtxs as usize;

    // Integer target weights
    let tpwgts0 = (graph.tvwgt[0] as Real * ntpwgts[0]) as Idx;
    let tpwgts1 = graph.tvwgt[0] - tpwgts0;
    let tpwgts = [tpwgts0, tpwgts1];

    let mindiff = (tpwgts[0] - graph.pwgts[0]).abs();

    // Determine from (overweight) and to (underweight)
    let from: usize = if graph.pwgts[0] < tpwgts[0] { 1 } else { 0 };
    let to: usize = 1 - from;

    let mut queue = PQueue::new(nvtxs);
    let mut moved = vec![-1 as Idx; nvtxs];

    // Permute vertex indices: irandArrayPermute(nvtxs, perm, nvtxs/5, 1)
    let mut perm = vec![0 as Idx; nvtxs];
    if nvtxs < 10 {
        ctrl.rng.rand_array_permute(nvtxs, &mut perm, 0, true);
    } else {
        ctrl.rng.rand_array_permute_with_nshuffles(nvtxs, &mut perm, 0, nvtxs / 5, true);
    }

    // Insert ALL vertices from 'from' partition with small weight
    for ii in 0..nvtxs {
        let i = perm[ii] as usize;
        if graph.where_[i] as usize == from && graph.vwgt[i] <= mindiff {
            queue.insert(i as Idx, (graph.ed[i] - graph.id[i]) as f64);
        }
    }

    let mut mincut = graph.mincut;
    let mut nswaps: Idx = 0;

    // Greedy move loop
    loop {
        let higain = if let Some((v, _)) = queue.get_top() {
            v as usize
        } else {
            break;
        };

        // If moving this vertex would overload the 'to' partition, stop
        if graph.pwgts[to] + graph.vwgt[higain] > tpwgts[to] {
            break;
        }

        mincut -= graph.ed[higain] - graph.id[higain];
        graph.pwgts[to] += graph.vwgt[higain];
        graph.pwgts[from] -= graph.vwgt[higain];

        graph.where_[higain] = to as Idx;
        moved[higain] = nswaps;

        // Swap id and ed for the moved vertex
        let tmp = graph.id[higain];
        graph.id[higain] = graph.ed[higain];
        graph.ed[higain] = tmp;

        // Update boundary status of moved vertex (two separate checks, matching C METIS)
        if graph.ed[higain] == 0 && graph.bndptr[higain] != -1 && graph.xadj[higain] < graph.xadj[higain + 1] {
            graph.bnd_delete(higain);
        }
        if graph.ed[higain] > 0 && graph.bndptr[higain] == -1 {
            graph.bnd_insert(higain);
        }

        // Update neighbors
        for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
            let k = graph.adjncy[kk] as usize;
            let kwgt = if to == graph.where_[k] as usize {
                graph.adjwgt[kk]
            } else {
                -graph.adjwgt[kk]
            };

            graph.id[k] += kwgt;
            graph.ed[k] -= kwgt;

            // Update queue position (only rpqUpdate, no insert/delete)
            if moved[k] == -1 && graph.where_[k] as usize == from && graph.vwgt[k] <= mindiff {
                queue.update(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
            }

            // Update boundary information (separate from PQ)
            if graph.ed[k] == 0 && graph.bndptr[k] != -1 {
                graph.bnd_delete(k);
            } else if graph.ed[k] > 0 && graph.bndptr[k] == -1 {
                graph.bnd_insert(k);
            }
        }

        nswaps += 1;
    }

    graph.mincut = mincut;
}

/// McGeneral2WayBalance: greedy multi-constraint 2-way balancing.
///
/// For ncon > 1: inserts all vertices from the overweight partition that could help
/// reduce the imbalance, greedily moves them while the target side can absorb them.
fn mc_general_2way_balance(_ctrl: &mut Ctrl, graph: &mut GraphData, ntpwgts: &[Real]) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Integer target weights per constraint per partition
    let mut tpwgts = vec![0 as Idx; 2 * ncon];
    for j in 0..ncon {
        tpwgts[j] = (graph.tvwgt[j] as Real * ntpwgts[j]) as Idx;
        tpwgts[ncon + j] = graph.tvwgt[j] - tpwgts[j];
    }

    // Determine which partition is overweight (using first constraint that is out of balance)
    let from: usize = if graph.pwgts[0] > tpwgts[0] { 0 } else { 1 };
    let to: usize = 1 - from;

    // Build PQ from all vertices in the 'from' partition
    let mut queue = PQueue::new(nvtxs);

    for v in 0..nvtxs {
        if graph.where_[v] as usize == from {
            let gain = graph.ed[v] - graph.id[v];
            queue.insert(v as Idx, gain as f64);
        }
    }

    // Greedy move loop
    while let Some((higain_idx, _key)) = queue.get_top() {
        let higain = higain_idx as usize;

        // Check if moving this vertex would overload the 'to' partition in any constraint
        let mut fits = true;
        for j in 0..ncon {
            let new_to_wgt = graph.pwgts[to * ncon + j] + graph.vwgt[higain * ncon + j];
            if new_to_wgt > tpwgts[to * ncon + j] {
                fits = false;
                break;
            }
        }
        if !fits {
            break;
        }

        // Move vertex from 'from' to 'to'
        graph.mincut -= graph.ed[higain] - graph.id[higain];

        graph.where_[higain] = to as Idx;
        for j in 0..ncon {
            graph.pwgts[to * ncon + j] += graph.vwgt[higain * ncon + j];
            graph.pwgts[from * ncon + j] -= graph.vwgt[higain * ncon + j];
        }

        // Swap id and ed for the moved vertex
        let tmp = graph.id[higain];
        graph.id[higain] = graph.ed[higain];
        graph.ed[higain] = tmp;

        // Update boundary status of moved vertex
        if graph.ed[higain] == 0
            && graph.xadj[higain] < graph.xadj[higain + 1]
        {
            graph.bnd_delete(higain);
        } else if graph.bndptr[higain] == -1 {
            graph.bnd_insert(higain);
        }

        // Update neighbors
        for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
            let k = graph.adjncy[kk] as usize;
            let kwgt = if to == graph.where_[k] as usize {
                graph.adjwgt[kk]
            } else {
                -graph.adjwgt[kk]
            };

            graph.id[k] += kwgt;
            graph.ed[k] -= kwgt;

            // Update boundary and PQ for neighbor k
            if graph.bndptr[k] != -1 {
                // k is on boundary
                if graph.ed[k] == 0 && graph.xadj[k] < graph.xadj[k + 1] {
                    graph.bnd_delete(k);
                    if queue.contains(k as Idx) {
                        queue.delete(k as Idx);
                    }
                } else if queue.contains(k as Idx) {
                    queue.update(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
                }
            } else {
                // k is not on boundary
                if graph.ed[k] > 0 {
                    graph.bnd_insert(k);
                    if graph.where_[k] as usize == from {
                        queue.insert(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
                    }
                }
            }
        }
    }
}
