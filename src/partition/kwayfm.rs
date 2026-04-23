use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::util::pqueue::PQueue;

const BNDTYPE_REFINE: Idx = 1;
const BNDTYPE_BALANCE: Idx = 2;
pub const OMODE_REFINE: Idx = 1;
pub const OMODE_BALANCE: Idx = 2;

const VPQSTATUS_PRESENT: i8 = 1;
const VPQSTATUS_EXTRACTED: i8 = 2;
const VPQSTATUS_NOTPRESENT: i8 = 3;

#[allow(dead_code)]
/// Greedy k-way optimization entry point.
pub fn greedy_kway_optimize(ctrl: &mut Ctrl, graph: &mut GraphData, niter: Idx, _ffactor: Real, omode: Idx) {
    greedy_kway_cut_optimize(ctrl, graph, niter, omode);
}

/// Greedy k-way cut optimization matching METIS Greedy_KWayCutOptimize.
pub fn greedy_kway_cut_optimize(ctrl: &mut Ctrl, graph: &mut GraphData, niter: Idx, omode: Idx) {
    let nvtxs = graph.nvtxs as usize;
    let nparts = ctrl.nparts as usize;
    let bndtype = if omode == OMODE_REFINE { BNDTYPE_REFINE } else { BNDTYPE_BALANCE };

    if nvtxs == 0 {
        return;
    }

    // Setup weight intervals
    let ubfactor = if omode == OMODE_BALANCE {
        ctrl.ubfactors[0]
    } else {
        let cur_imbal = compute_load_imbalance_kway(graph, nparts, &ctrl.pijbm);
        ctrl.ubfactors[0].max(cur_imbal)
    };

    let mut minpwgts = vec![0 as Idx; nparts];
    let mut maxpwgts = vec![0 as Idx; nparts];
    for i in 0..nparts {
        maxpwgts[i] = (ctrl.tpwgts[i] * graph.tvwgt[0] as Real * ubfactor) as Idx;
        minpwgts[i] = (ctrl.tpwgts[i] * graph.tvwgt[0] as Real * (1.0 / ubfactor)) as Idx;
    }

    let mut perm = vec![0 as Idx; nvtxs];
    let mut vstatus = vec![VPQSTATUS_NOTPRESENT; nvtxs];
    let mut updptr = vec![-1 as Idx; nvtxs];
    let mut updind = vec![0 as Idx; nvtxs];

    let mut queue = PQueue::new(nvtxs);

    for _pass in 0..niter {
        if omode == OMODE_BALANCE {
            let mut balanced = true;
            for i in 0..nparts {
                if graph.pwgts[i] > maxpwgts[i] || graph.pwgts[i] < minpwgts[i] {
                    balanced = false;
                    break;
                }
            }
            if balanced {
                break;
            }
        }

        let oldcut = graph.mincut;
        let mut nbnd = graph.nbnd;
        let mut nupd: usize = 0;

        // Insert boundary vertices in PQ
        let nbnd_usize = nbnd as usize;
        for i in 0..nbnd_usize {
            perm[i] = i as Idx;
        }
        if nbnd_usize > 0 {
            let nshuffles = nbnd_usize / 4;
            if nbnd_usize < 10 {
                ctrl.rng.rand_array_permute(nbnd_usize, &mut perm, 0, true);
            } else {
                ctrl.rng.rand_array_permute_with_nshuffles(nbnd_usize, &mut perm, 0, nshuffles, true);
            }
        }

        for ii in 0..nbnd_usize {
            let i = graph.bndind[perm[ii] as usize] as usize;
            if i >= nvtxs { continue; }
            let ckr = &graph.ckrinfo[i];
            let rgain = if ckr.nnbrs > 0 {
                (ckr.ed as f64) / (ckr.nnbrs as f64).sqrt() - ckr.id as f64
            } else {
                -(ckr.id as f64)
            };
            queue.insert(i as Idx, rgain);
            vstatus[i] = VPQSTATUS_PRESENT;
            updind[nupd] = i as Idx;
            updptr[i] = nupd as Idx;
            nupd += 1;
        }

        let mut nmoved = 0;
        loop {
            let i = if let Some((v, _)) = queue.get_top() {
                v as usize
            } else {
                break;
            };
            vstatus[i] = VPQSTATUS_EXTRACTED;

            let from = graph.where_[i] as usize;
            let vwgt = graph.vwgt[i];

            let ckr = &graph.ckrinfo[i];
            let inbr = ckr.inbr;
            if inbr < 0 || ckr.nnbrs == 0 {
                continue;
            }

            let mut best_k: i32 = -1;

            if omode == OMODE_REFINE {
                for kk in (0..ckr.nnbrs as usize).rev() {
                    let nbr = &ctrl.cnbrpool[inbr as usize + kk];
                    let to = nbr.pid as usize;
                    if to >= nparts { continue; }
                    if (nbr.ed > ckr.id
                        && (graph.pwgts[from] - vwgt >= minpwgts[from]
                            || (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real))
                        && (graph.pwgts[to] + vwgt <= maxpwgts[to]
                            || (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real)))
                        || (nbr.ed == ckr.id
                            && (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real))
                    {
                        best_k = kk as i32;
                        break;
                    }
                }
                if best_k < 0 {
                    continue;
                }

                let bk = best_k as usize;
                for jj in (0..bk).rev() {
                    let nbr_j = &ctrl.cnbrpool[inbr as usize + jj];
                    let to = nbr_j.pid as usize;
                    if to >= nparts { continue; }
                    let nbr_k = &ctrl.cnbrpool[inbr as usize + best_k as usize];
                    if (nbr_j.ed > nbr_k.ed
                        && (graph.pwgts[from] - vwgt >= minpwgts[from]
                            || (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real))
                        && (graph.pwgts[to] + vwgt <= maxpwgts[to]
                            || (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real)))
                        || (nbr_j.ed == nbr_k.ed
                            && (ctrl.tpwgts[nbr_k.pid as usize] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * graph.pwgts[nbr_k.pid as usize] as Real))
                    {
                        best_k = jj as i32;
                    }
                }
            } else {
                // OMODE_BALANCE
                for kk in (0..ckr.nnbrs as usize).rev() {
                    let nbr = &ctrl.cnbrpool[inbr as usize + kk];
                    let to = nbr.pid as usize;
                    if to >= nparts { continue; }
                    if from >= nparts
                        || (ctrl.tpwgts[from] * graph.pwgts[to] as Real) < (ctrl.tpwgts[to] * (graph.pwgts[from] - vwgt) as Real)
                    {
                        best_k = kk as i32;
                        break;
                    }
                }
                if best_k < 0 {
                    continue;
                }

                let bk = best_k as usize;
                for jj in (0..bk).rev() {
                    let nbr_j = &ctrl.cnbrpool[inbr as usize + jj];
                    let to = nbr_j.pid as usize;
                    if to >= nparts { continue; }
                    let nbr_k = &ctrl.cnbrpool[inbr as usize + best_k as usize];
                    if (ctrl.tpwgts[nbr_k.pid as usize] * graph.pwgts[to] as Real)
                        < (ctrl.tpwgts[to] * graph.pwgts[nbr_k.pid as usize] as Real)
                    {
                        best_k = jj as i32;
                    }
                }
            }

            let k = best_k as usize;
            let ed_k = ctrl.cnbrpool[inbr as usize + k].ed;
            let to = ctrl.cnbrpool[inbr as usize + k].pid as usize;

            let gain = ed_k - graph.ckrinfo[i].id;
            graph.mincut -= gain;
            nmoved += 1;

            graph.pwgts[to] += vwgt;
            graph.pwgts[from] -= vwgt;

            // UpdateMovedVertexInfoAndBND
            update_moved_vertex(graph, ctrl, i, from, k, to, &mut nbnd, bndtype);

            // Update adjacent vertices
            for jj in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let ii = graph.adjncy[jj] as usize;
                let me = graph.where_[ii] as usize;
                let ewgt = graph.adjwgt[jj];
                let old_nnbrs = graph.ckrinfo[ii].nnbrs;

                update_adjacent_vertex(ctrl, graph, ii, me, from, to, ewgt, &mut nbnd, bndtype);

                let new_nnbrs = graph.ckrinfo[ii].nnbrs;
                if me == to || me == from || old_nnbrs != new_nnbrs {
                    let ckr_ii = &graph.ckrinfo[ii];
                    let rgain = if ckr_ii.nnbrs > 0 {
                        (ckr_ii.ed as f64) / (ckr_ii.nnbrs as f64).sqrt() - ckr_ii.id as f64
                    } else {
                        -(ckr_ii.id as f64)
                    };

                    if bndtype == BNDTYPE_REFINE {
                        if vstatus[ii] == VPQSTATUS_PRESENT {
                            if ckr_ii.ed - ckr_ii.id >= 0 {
                                queue.update(ii as Idx, rgain);
                            } else {
                                queue.delete(ii as Idx);
                                vstatus[ii] = VPQSTATUS_NOTPRESENT;
                                if updptr[ii] >= 0 {
                                    let pos = updptr[ii] as usize;
                                    nupd -= 1;
                                    updind[pos] = updind[nupd];
                                    if (updind[nupd] as usize) < nvtxs {
                                        updptr[updind[nupd] as usize] = pos as Idx;
                                    }
                                    updptr[ii] = -1;
                                }
                            }
                        } else if vstatus[ii] == VPQSTATUS_NOTPRESENT && ckr_ii.ed - ckr_ii.id >= 0 {
                            queue.insert(ii as Idx, rgain);
                            vstatus[ii] = VPQSTATUS_PRESENT;
                            updind[nupd] = ii as Idx;
                            updptr[ii] = nupd as Idx;
                            nupd += 1;
                        }
                    } else {
                        if vstatus[ii] == VPQSTATUS_PRESENT {
                            if ckr_ii.ed > 0 {
                                queue.update(ii as Idx, rgain);
                            } else {
                                queue.delete(ii as Idx);
                                vstatus[ii] = VPQSTATUS_NOTPRESENT;
                                if updptr[ii] >= 0 {
                                    let pos = updptr[ii] as usize;
                                    nupd -= 1;
                                    updind[pos] = updind[nupd];
                                    if (updind[nupd] as usize) < nvtxs {
                                        updptr[updind[nupd] as usize] = pos as Idx;
                                    }
                                    updptr[ii] = -1;
                                }
                            }
                        } else if vstatus[ii] == VPQSTATUS_NOTPRESENT && ckr_ii.ed > 0 {
                            queue.insert(ii as Idx, rgain);
                            vstatus[ii] = VPQSTATUS_PRESENT;
                            updind[nupd] = ii as Idx;
                            updptr[ii] = nupd as Idx;
                            nupd += 1;
                        }
                    }
                }
            }
        }

        graph.nbnd = nbnd;

        // Reset vstatus
        for idx in 0..nupd {
            let v = updind[idx] as usize;
            if v < nvtxs {
                vstatus[v] = VPQSTATUS_NOTPRESENT;
                updptr[v] = -1;
            }
        }

        if nmoved == 0 || (omode == OMODE_REFINE && graph.mincut == oldcut) {
            break;
        }
    }
}

fn update_moved_vertex(
    graph: &mut GraphData,
    ctrl: &mut Ctrl,
    i: usize,
    from: usize,
    k: usize,
    to: usize,
    nbnd: &mut Idx,
    bndtype: Idx,
) {
    let inbr = graph.ckrinfo[i].inbr as usize;
    let ed_k = ctrl.cnbrpool[inbr + k].ed;

    graph.where_[i] = to as Idx;
    graph.ckrinfo[i].ed += graph.ckrinfo[i].id - ed_k;

    let tmp = graph.ckrinfo[i].id;
    graph.ckrinfo[i].id = ed_k;
    ctrl.cnbrpool[inbr + k].ed = tmp;

    if ctrl.cnbrpool[inbr + k].ed == 0 {
        let nnbrs = graph.ckrinfo[i].nnbrs as usize;
        if nnbrs > 0 {
            graph.ckrinfo[i].nnbrs -= 1;
            let last = graph.ckrinfo[i].nnbrs as usize;
            ctrl.cnbrpool[inbr + k] = ctrl.cnbrpool[inbr + last];
        }
    } else {
        ctrl.cnbrpool[inbr + k].pid = from as Idx;
    }

    if bndtype == BNDTYPE_REFINE {
        if graph.bndptr[i] != -1 && graph.ckrinfo[i].ed - graph.ckrinfo[i].id < 0 {
            graph.bnd_delete(i);
            *nbnd -= 1;
        }
        if graph.bndptr[i] == -1 && graph.ckrinfo[i].ed - graph.ckrinfo[i].id >= 0 {
            graph.bnd_insert(i);
            *nbnd += 1;
        }
    } else {
        if graph.bndptr[i] != -1 && graph.ckrinfo[i].ed <= 0 {
            graph.bnd_delete(i);
            *nbnd -= 1;
        }
        if graph.bndptr[i] == -1 && graph.ckrinfo[i].ed > 0 {
            graph.bnd_insert(i);
            *nbnd += 1;
        }
    }
}

fn update_adjacent_vertex(
    ctrl: &mut Ctrl,
    graph: &mut GraphData,
    vid: usize,
    me: usize,
    from: usize,
    to: usize,
    ewgt: Idx,
    nbnd: &mut Idx,
    bndtype: Idx,
) {
    let adjlen = (graph.xadj[vid + 1] - graph.xadj[vid]) as usize;

    if graph.ckrinfo[vid].inbr == -1 {
        graph.ckrinfo[vid].inbr = ctrl.cnbrpool_get_next(adjlen);
        graph.ckrinfo[vid].nnbrs = 0;
    }

    let inbr = graph.ckrinfo[vid].inbr as usize;

    if me == from {
        graph.ckrinfo[vid].ed += ewgt;
        graph.ckrinfo[vid].id -= ewgt;
        if bndtype == BNDTYPE_REFINE {
            if graph.ckrinfo[vid].ed - graph.ckrinfo[vid].id >= 0 && graph.bndptr[vid] == -1 {
                graph.bnd_insert(vid);
                *nbnd += 1;
            }
        } else {
            if graph.ckrinfo[vid].ed > 0 && graph.bndptr[vid] == -1 {
                graph.bnd_insert(vid);
                *nbnd += 1;
            }
        }
    } else if me == to {
        graph.ckrinfo[vid].id += ewgt;
        graph.ckrinfo[vid].ed -= ewgt;
        if bndtype == BNDTYPE_REFINE {
            if graph.ckrinfo[vid].ed - graph.ckrinfo[vid].id < 0 && graph.bndptr[vid] != -1 {
                graph.bnd_delete(vid);
                *nbnd -= 1;
            }
        } else {
            if graph.ckrinfo[vid].ed <= 0 && graph.bndptr[vid] != -1 {
                graph.bnd_delete(vid);
                *nbnd -= 1;
            }
        }
    }

    // Remove contribution from 'from'
    if me != from {
        let nnbrs = graph.ckrinfo[vid].nnbrs as usize;
        for kk in 0..nnbrs {
            if ctrl.cnbrpool[inbr + kk].pid == from as Idx {
                if ctrl.cnbrpool[inbr + kk].ed == ewgt {
                    graph.ckrinfo[vid].nnbrs -= 1;
                    let last = graph.ckrinfo[vid].nnbrs as usize;
                    ctrl.cnbrpool[inbr + kk] = ctrl.cnbrpool[inbr + last];
                } else {
                    ctrl.cnbrpool[inbr + kk].ed -= ewgt;
                }
                break;
            }
        }
    }

    // Add contribution to 'to'
    if me != to {
        let nnbrs = graph.ckrinfo[vid].nnbrs as usize;
        let mut found = false;
        for kk in 0..nnbrs {
            if ctrl.cnbrpool[inbr + kk].pid == to as Idx {
                ctrl.cnbrpool[inbr + kk].ed += ewgt;
                found = true;
                break;
            }
        }
        if !found {
            let kk = graph.ckrinfo[vid].nnbrs as usize;
            ctrl.cnbrpool[inbr + kk].pid = to as Idx;
            ctrl.cnbrpool[inbr + kk].ed = ewgt;
            graph.ckrinfo[vid].nnbrs += 1;
        }
    }
}

fn compute_load_imbalance_kway(graph: &GraphData, nparts: usize, pijbm: &[Real]) -> Real {
    let mut max_imbal: Real = 0.0;
    for i in 0..nparts {
        if i < graph.pwgts.len() && i < pijbm.len() {
            let imbal = graph.pwgts[i] as Real * pijbm[i];
            if imbal > max_imbal {
                max_imbal = imbal;
            }
        }
    }
    max_imbal
}
