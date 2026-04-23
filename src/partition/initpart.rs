use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;

/// Initialize a 2-way partition on the coarsest graph. Tries niparts attempts.
pub fn init_2way_partition(ctrl: &mut Ctrl, graph: &mut GraphData, tpwgts: &[Real], niparts: Idx) {
    graph.alloc_2way();
    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; graph.nvtxs as usize];

    for _inbfs in 0..niparts {
        if ctrl.iptype == 1 {
            random_bisection(ctrl, graph, tpwgts);
        } else {
            grow_bisection(ctrl, graph, tpwgts);
        }

        compute_2way_partition_params(ctrl, graph);
        super::balance::balance_2way(ctrl, graph, tpwgts);
        super::fm::fm_2way_cut_refine(ctrl, graph, tpwgts, ctrl.niter);

        if graph.mincut < bestcut {
            bestcut = graph.mincut;
            bestwhere.copy_from_slice(&graph.where_);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.where_.copy_from_slice(&bestwhere);
    compute_2way_partition_params(ctrl, graph);
}

/// GrowBisection: BFS from random seed, growing partition 0.
fn grow_bisection(ctrl: &mut Ctrl, graph: &mut GraphData, tpwgts: &[Real]) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Target weights for partition 1
    let onemaxpwgt = (ctrl.ubfactors[0] * graph.tvwgt[0] as Real * tpwgts[ncon]) as Idx;
    let oneminpwgt = ((1.0 / ctrl.ubfactors[0]) * graph.tvwgt[0] as Real * tpwgts[ncon]) as Idx;

    // Initialize all to partition 1
    for i in 0..nvtxs {
        graph.where_[i] = 1;
    }

    let mut touched = vec![0u8; nvtxs];
    let mut queue = vec![0usize; nvtxs];
    let mut pwgts = [0 as Idx; 2];
    pwgts[1] = graph.tvwgt[0];

    // Pick random seed
    let start = ctrl.rng.rand_in_range(graph.nvtxs) as usize;
    queue[0] = start;
    touched[start] = 1;
    let mut first = 0usize;
    let mut last = 1usize;
    let mut nleft = nvtxs - 1;
    let mut drain = false;

    loop {
        if first == last {
            if nleft == 0 || drain {
                break;
            }
            // Find random untouched vertex
            let k = ctrl.rng.rand_in_range(nleft as Idx) as usize;
            let mut cnt = 0;
            for i in 0..nvtxs {
                if touched[i] == 0 {
                    if cnt == k {
                        queue[0] = i;
                        touched[i] = 1;
                        first = 0;
                        last = 1;
                        nleft -= 1;
                        break;
                    }
                    cnt += 1;
                }
            }
        }

        let i = queue[first];
        first += 1;

        // Minimum weight guard
        if pwgts[0] > 0 && pwgts[1] - graph.vwgt[i] < oneminpwgt {
            drain = true;
            continue;
        }

        // Move to partition 0
        graph.where_[i] = 0;
        pwgts[0] += graph.vwgt[i];
        pwgts[1] -= graph.vwgt[i];

        if pwgts[1] <= onemaxpwgt {
            break;
        }

        drain = false;

        // Expand BFS
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let j = graph.adjncy[k] as usize;
            if touched[j] == 0 {
                touched[j] = 1;
                queue[last] = j;
                last += 1;
                nleft -= 1;
            }
        }
    }

    // Degenerate cases
    if pwgts[1] == 0 {
        let v = ctrl.rng.rand_in_range(graph.nvtxs) as usize;
        graph.where_[v] = 1;
    }
    if pwgts[0] == 0 {
        let v = ctrl.rng.rand_in_range(graph.nvtxs) as usize;
        graph.where_[v] = 0;
    }
}

/// RandomBisection: randomly assign vertices to partition 0.
fn random_bisection(ctrl: &mut Ctrl, graph: &mut GraphData, tpwgts: &[Real]) {
    let nvtxs = graph.nvtxs as usize;

    let zeromaxpwgt = (ctrl.ubfactors[0] * graph.tvwgt[0] as Real * tpwgts[0]) as Idx;

    for i in 0..nvtxs {
        graph.where_[i] = 1;
    }

    if nvtxs < 10 {
        let mut perm = vec![0 as Idx; nvtxs];
        ctrl.rng.rand_array_permute(nvtxs, &mut perm, 0, true);

        let mut pwgt0 = 0 as Idx;
        for &pi in &perm {
            let i = pi as usize;
            if pwgt0 + graph.vwgt[i] <= zeromaxpwgt {
                graph.where_[i] = 0;
                pwgt0 += graph.vwgt[i];
            }
        }
    } else {
        let mut perm = vec![0 as Idx; nvtxs];
        let nshuffles = nvtxs / 2;
        ctrl.rng.rand_array_permute_with_nshuffles(nvtxs, &mut perm, 0, nshuffles, true);

        let mut pwgt0 = 0 as Idx;
        for &pi in &perm {
            let i = pi as usize;
            if pwgt0 + graph.vwgt[i] <= zeromaxpwgt {
                graph.where_[i] = 0;
                pwgt0 += graph.vwgt[i];
            }
        }
    }
}

/// Compute 2-way partition parameters from scratch.
pub fn compute_2way_partition_params(_ctrl: &Ctrl, graph: &mut GraphData) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Reset
    graph.pwgts = vec![0; 2 * ncon];
    graph.bndptr = vec![-1; nvtxs];
    graph.bndind = vec![0; nvtxs];
    graph.nbnd = 0;
    graph.id = vec![0; nvtxs];
    graph.ed = vec![0; nvtxs];

    // Compute pwgts
    if ncon == 1 {
        for i in 0..nvtxs {
            graph.pwgts[graph.where_[i] as usize] += graph.vwgt[i];
        }
    } else {
        for i in 0..nvtxs {
            let p = graph.where_[i] as usize;
            for j in 0..ncon {
                graph.pwgts[p * ncon + j] += graph.vwgt[i * ncon + j];
            }
        }
    }

    // Compute id, ed, boundary, mincut
    let mut mincut: Idx = 0;

    for i in 0..nvtxs {
        let me = graph.where_[i];
        let mut tid: Idx = 0;
        let mut ted: Idx = 0;

        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;

        for k in istart..iend {
            if graph.where_[graph.adjncy[k] as usize] == me {
                tid += graph.adjwgt[k];
            } else {
                ted += graph.adjwgt[k];
            }
        }

        graph.id[i] = tid;
        graph.ed[i] = ted;

        if ted > 0 || istart == iend {
            graph.bnd_insert(i);
            mincut += ted;
        }
    }

    graph.mincut = mincut / 2;
}
