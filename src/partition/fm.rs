use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::util::pqueue::PQueue;

/// FM 2-way cut refinement with rollback (single constraint).
/// Matches METIS FM_2WayCutRefine.
#[allow(unused_assignments)]
pub fn fm_2way_cut_refine(ctrl: &mut Ctrl, graph: &mut GraphData, tpwgts: &[Real], niter: Idx) {
    let nvtxs = graph.nvtxs as usize;
    if graph.nbnd <= 0 {
        return;
    }

    // Integer target weights
    let itpwgts = [
        (graph.tvwgt[0] as Real * tpwgts[0]) as Idx,
        graph.tvwgt[0] - (graph.tvwgt[0] as Real * tpwgts[0]) as Idx,
    ];

    let limit = (0.01 * nvtxs as f64).max(15.0).min(100.0) as usize;
    let avgvwgt = ((graph.pwgts[0] + graph.pwgts[1]) / 20)
        .min(2 * (graph.pwgts[0] + graph.pwgts[1]) / nvtxs.max(1) as Idx);
    let origdiff = (itpwgts[0] - graph.pwgts[0]).abs();

    let mut moved = vec![-1 as Idx; nvtxs];
    let mut swaps = vec![0 as Idx; nvtxs];

    for _pass in 0..niter {
        let initcut = graph.mincut;
        let mut mincut = initcut;
        let mut newcut = initcut;
        let mut mincutorder: i32 = -1;
        let mut mindiff = (itpwgts[0] - graph.pwgts[0]).abs();

        let mut queues = [PQueue::new(nvtxs), PQueue::new(nvtxs)];

        // Insert boundary vertices into PQs
        let nbnd = graph.nbnd as usize;
        let mut perm = vec![0 as Idx; nbnd];
        if nbnd < 10 {
            ctrl.rng.rand_array_permute(nbnd, &mut perm, 0, true);
        } else {
            // C METIS FM_2WayCutRefine: irandArrayPermute(nbnd, perm, nbnd, 1)
            ctrl.rng.rand_array_permute_with_nshuffles(nbnd, &mut perm, 0, nbnd, true);
        }

        for ii in 0..nbnd {
            let v = graph.bndind[perm[ii] as usize] as usize;
            let gain = graph.ed[v] - graph.id[v];
            queues[graph.where_[v] as usize].insert(v as Idx, gain as f64);
        }

        let mut nswaps = 0usize;

        loop {
            // Select source partition (overweight side)
            let from = if itpwgts[0] - graph.pwgts[0] < itpwgts[1] - graph.pwgts[1] {
                0
            } else {
                1
            };
            let to = 1 - from;

            let higain = if let Some((v, _)) = queues[from].get_top() {
                v as usize
            } else {
                break;
            };

            newcut -= graph.ed[higain] - graph.id[higain];

            // Update partition weights
            graph.pwgts[to] += graph.vwgt[higain];
            graph.pwgts[from] -= graph.vwgt[higain];

            // Check if this is a new best
            let newdiff = (itpwgts[0] - graph.pwgts[0]).abs();
            if (newcut < mincut && newdiff <= origdiff + avgvwgt)
                || (newcut == mincut && newdiff < mindiff)
            {
                mincut = newcut;
                mindiff = newdiff;
                mincutorder = nswaps as i32;
            } else if nswaps as i32 - mincutorder > limit as i32 {
                // Undo this move and stop
                graph.pwgts[from] += graph.vwgt[higain];
                graph.pwgts[to] -= graph.vwgt[higain];
                newcut += graph.ed[higain] - graph.id[higain];
                break;
            }

            // Commit the move
            graph.where_[higain] = to as Idx;
            moved[higain] = nswaps as Idx;
            swaps[nswaps] = higain as Idx;

            // Swap id and ed for the moved vertex
            let tmp = graph.id[higain];
            graph.id[higain] = graph.ed[higain];
            graph.ed[higain] = tmp;

            // Update boundary of moved vertex
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

                // Update boundary and PQ for k
                if graph.bndptr[k] != -1 {
                    // k is currently on boundary
                    if graph.ed[k] == 0 && graph.xadj[k] < graph.xadj[k + 1] {
                        graph.bnd_delete(k);
                        if moved[k] == -1 {
                            queues[graph.where_[k] as usize].delete(k as Idx);
                        }
                    } else if moved[k] == -1 {
                        queues[graph.where_[k] as usize].update(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
                    }
                } else {
                    // k is not on boundary
                    if graph.ed[k] > 0 {
                        graph.bnd_insert(k);
                        if moved[k] == -1 {
                            queues[graph.where_[k] as usize].insert(k as Idx, (graph.ed[k] - graph.id[k]) as f64);
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

        let rollback_start = if mincutorder < 0 { 0 } else { mincutorder as usize + 1 };
        for i in (rollback_start..nswaps).rev() {
            let higain = swaps[i] as usize;
            let to = graph.where_[higain] as usize;
            let from = 1 - to;

            // Undo the move
            graph.where_[higain] = from as Idx;

            // Swap id/ed back
            let tmp = graph.id[higain];
            graph.id[higain] = graph.ed[higain];
            graph.ed[higain] = tmp;

            graph.pwgts[from] += graph.vwgt[higain];
            graph.pwgts[to] -= graph.vwgt[higain];

            // Fix boundary
            if graph.ed[higain] == 0 && graph.xadj[higain] < graph.xadj[higain + 1] {
                graph.bnd_delete(higain);
            } else if graph.bndptr[higain] == -1 {
                graph.bnd_insert(higain);
            }

            // Update neighbors
            for kk in graph.xadj[higain] as usize..graph.xadj[higain + 1] as usize {
                let k = graph.adjncy[kk] as usize;
                let kwgt = if from == graph.where_[k] as usize {
                    graph.adjwgt[kk]
                } else {
                    -graph.adjwgt[kk]
                };

                graph.id[k] += kwgt;
                graph.ed[k] -= kwgt;

                if graph.bndptr[k] != -1 && graph.ed[k] == 0 && graph.xadj[k] < graph.xadj[k + 1] {
                    graph.bnd_delete(k);
                } else if graph.bndptr[k] == -1 && graph.ed[k] > 0 {
                    graph.bnd_insert(k);
                }
            }
        }

        graph.mincut = mincut;

        if mincutorder <= 0 || mincut == initcut {
            break;
        }
    }
}
