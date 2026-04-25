use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::types::{Idx, Real};
use crate::util::pqueue::PQueue;

const BOUNDARY_REFINE: Idx = 1;
const BOUNDARY_BALANCE: Idx = 2;
pub const MODE_REFINE: Idx = 1;
pub const MODE_BALANCE: Idx = 2;

const PQ_PRESENT: i8 = 1;
const PQ_EXTRACTED: i8 = 2;
const PQ_NOT_PRESENT: i8 = 3;

#[expect(dead_code)]
/// Greedy k-way optimization entry point.
pub fn greedy_kway_optimize(
    ctrl: &mut Control,
    graph: &mut GraphData,
    niter: Idx,
    _ffactor: Real,
    omode: Idx,
) {
    greedy_kway_cut_optimize(ctrl, graph, niter, omode);
}

/// Greedy k-way cut optimization matching METIS Greedy_KWayOptimize dispatcher.
pub fn greedy_kway_cut_optimize(ctrl: &mut Control, graph: &mut GraphData, niter: Idx, omode: Idx) {
    if graph.num_constraints > 1 {
        greedy_mc_kway_cut_optimize(ctrl, graph, niter, omode);
        return;
    }

    let num_vertices = graph.num_vertices as usize;
    let nparts = ctrl.num_parts as usize;
    let bndtype = if omode == MODE_REFINE {
        BOUNDARY_REFINE
    } else {
        BOUNDARY_BALANCE
    };

    if num_vertices == 0 {
        return;
    }

    // Setup weight intervals
    let ubfactor = if omode == MODE_BALANCE {
        ctrl.imbalance_tols[0]
    } else {
        let cur_imbal =
            compute_load_imbalance_kway(graph, nparts, &ctrl.partition_ij_balance_multipliers);
        ctrl.imbalance_tols[0].max(cur_imbal)
    };

    let mut minpart_weights = vec![0 as Idx; nparts];
    let mut maxpart_weights = vec![0 as Idx; nparts];
    for i in 0..nparts {
        maxpart_weights[i] =
            (ctrl.target_part_weights[i] * graph.total_vertex_weight[0] as Real * ubfactor) as Idx;
        minpart_weights[i] = (ctrl.target_part_weights[i]
            * graph.total_vertex_weight[0] as Real
            * (1.0 / ubfactor)) as Idx;
    }

    let mut perm = vec![0 as Idx; num_vertices];
    let mut vertex_status = vec![PQ_NOT_PRESENT; num_vertices];
    let mut update_index = vec![-1 as Idx; num_vertices];
    let mut update_list = vec![0 as Idx; num_vertices];

    let mut queue = PQueue::new(num_vertices);

    for _pass in 0..niter {
        if omode == MODE_BALANCE {
            let mut balanced = true;
            for i in 0..nparts {
                if graph.part_weights[i] > maxpart_weights[i]
                    || graph.part_weights[i] < minpart_weights[i]
                {
                    balanced = false;
                    break;
                }
            }
            if balanced {
                break;
            }
        }

        let oldcut = graph.edge_cut;
        let mut num_boundary = graph.num_boundary;
        let mut num_updates: usize = 0;

        // Insert boundary vertices in PQ
        let num_boundary_usize = num_boundary as usize;
        #[allow(clippy::needless_range_loop)]
        for i in 0..num_boundary_usize {
            perm[i] = i as Idx;
        }
        if num_boundary_usize > 0 {
            let nshuffles = num_boundary_usize / 4;
            if num_boundary_usize < 10 {
                ctrl.rng
                    .rand_array_permute(num_boundary_usize, &mut perm, 0, true);
            } else {
                ctrl.rng.rand_array_permute_with_nshuffles(
                    num_boundary_usize,
                    &mut perm,
                    0,
                    nshuffles,
                    true,
                );
            }
        }

        for &p in &perm[..num_boundary_usize] {
            let i = graph.boundary_list[p as usize] as usize;
            if i >= num_vertices {
                continue;
            }
            let ckr = &graph.kway_refinement_info[i];
            let scaled_gain =
                compute_kway_gain(ckr.external_degree, ckr.num_neighbors, ckr.internal_degree);
            queue.insert(i as Idx, scaled_gain);
            vertex_status[i] = PQ_PRESENT;
            update_list[num_updates] = i as Idx;
            update_index[i] = num_updates as Idx;
            num_updates += 1;
        }

        let mut nmoved = 0;
        while let Some((v, _)) = queue.get_top() {
            let i = v as usize;
            vertex_status[i] = PQ_EXTRACTED;

            let from = graph.partition[i] as usize;
            let vertex_weights = graph.vertex_weights[i];

            let ckr = &graph.kway_refinement_info[i];
            let inbr = ckr.neighbor_offset;
            if inbr < 0 || ckr.num_neighbors == 0 {
                continue;
            }

            let mut best_k: i32 = -1;

            if omode == MODE_REFINE {
                for kk in (0..ckr.num_neighbors as usize).rev() {
                    let nbr = &ctrl.neighbor_pool[inbr as usize + kk];
                    let to = nbr.part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    if (nbr.external_degree > ckr.internal_degree
                        && (graph.part_weights[from] - vertex_weights >= minpart_weights[from]
                            || (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * (graph.part_weights[from] - vertex_weights) as Real))
                        && (graph.part_weights[to] + vertex_weights <= maxpart_weights[to]
                            || (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * (graph.part_weights[from] - vertex_weights) as Real)))
                        || (nbr.external_degree == ckr.internal_degree
                            && (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * (graph.part_weights[from] - vertex_weights) as Real))
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
                    let nbr_j = &ctrl.neighbor_pool[inbr as usize + jj];
                    let to = nbr_j.part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let nbr_k = &ctrl.neighbor_pool[inbr as usize + best_k as usize];
                    if (nbr_j.external_degree > nbr_k.external_degree
                        && (graph.part_weights[from] - vertex_weights >= minpart_weights[from]
                            || (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * (graph.part_weights[from] - vertex_weights) as Real))
                        && (graph.part_weights[to] + vertex_weights <= maxpart_weights[to]
                            || (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * (graph.part_weights[from] - vertex_weights) as Real)))
                        || (nbr_j.external_degree == nbr_k.external_degree
                            && (ctrl.target_part_weights[nbr_k.part_id as usize]
                                * graph.part_weights[to] as Real)
                                < (ctrl.target_part_weights[to]
                                    * graph.part_weights[nbr_k.part_id as usize] as Real))
                    {
                        best_k = jj as i32;
                    }
                }
            } else {
                // MODE_BALANCE
                for kk in (0..ckr.num_neighbors as usize).rev() {
                    let nbr = &ctrl.neighbor_pool[inbr as usize + kk];
                    let to = nbr.part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    if from >= nparts
                        || (ctrl.target_part_weights[from] * graph.part_weights[to] as Real)
                            < (ctrl.target_part_weights[to]
                                * (graph.part_weights[from] - vertex_weights) as Real)
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
                    let nbr_j = &ctrl.neighbor_pool[inbr as usize + jj];
                    let to = nbr_j.part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let nbr_k = &ctrl.neighbor_pool[inbr as usize + best_k as usize];
                    if (ctrl.target_part_weights[nbr_k.part_id as usize]
                        * graph.part_weights[to] as Real)
                        < (ctrl.target_part_weights[to]
                            * graph.part_weights[nbr_k.part_id as usize] as Real)
                    {
                        best_k = jj as i32;
                    }
                }
            }

            let k = best_k as usize;
            let ed_k = ctrl.neighbor_pool[inbr as usize + k].external_degree;
            let to = ctrl.neighbor_pool[inbr as usize + k].part_id as usize;

            let gain = ed_k - graph.kway_refinement_info[i].internal_degree;
            graph.edge_cut -= gain;
            nmoved += 1;

            graph.part_weights[to] += vertex_weights;
            graph.part_weights[from] -= vertex_weights;

            // UpdateMovedVertexInfoAndBND
            update_moved_vertex(graph, ctrl, i, from, k, to, &mut num_boundary, bndtype);

            // Update adjacent vertices
            for jj in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let ii = graph.adjacency[jj] as usize;
                let me = graph.partition[ii] as usize;
                let ewgt = graph.edge_weights[jj];
                let old_nnbrs = graph.kway_refinement_info[ii].num_neighbors;

                update_adjacent_vertex(
                    ctrl,
                    graph,
                    ii,
                    me,
                    from,
                    to,
                    ewgt,
                    &mut num_boundary,
                    bndtype,
                );

                let new_nnbrs = graph.kway_refinement_info[ii].num_neighbors;
                if me == to || me == from || old_nnbrs != new_nnbrs {
                    let ckr_ii = &graph.kway_refinement_info[ii];
                    let scaled_gain = compute_kway_gain(
                        ckr_ii.external_degree,
                        ckr_ii.num_neighbors,
                        ckr_ii.internal_degree,
                    );

                    if bndtype == BOUNDARY_REFINE {
                        if vertex_status[ii] == PQ_PRESENT {
                            if ckr_ii.external_degree - ckr_ii.internal_degree >= 0 {
                                queue.update(ii as Idx, scaled_gain);
                            } else {
                                queue.delete(ii as Idx);
                                vertex_status[ii] = PQ_NOT_PRESENT;
                                if update_index[ii] >= 0 {
                                    let pos = update_index[ii] as usize;
                                    num_updates -= 1;
                                    update_list[pos] = update_list[num_updates];
                                    if (update_list[num_updates] as usize) < num_vertices {
                                        update_index[update_list[num_updates] as usize] =
                                            pos as Idx;
                                    }
                                    update_index[ii] = -1;
                                }
                            }
                        } else if vertex_status[ii] == PQ_NOT_PRESENT
                            && ckr_ii.external_degree - ckr_ii.internal_degree >= 0
                        {
                            queue.insert(ii as Idx, scaled_gain);
                            vertex_status[ii] = PQ_PRESENT;
                            update_list[num_updates] = ii as Idx;
                            update_index[ii] = num_updates as Idx;
                            num_updates += 1;
                        }
                    } else if vertex_status[ii] == PQ_PRESENT {
                        if ckr_ii.external_degree > 0 {
                            queue.update(ii as Idx, scaled_gain);
                        } else {
                            queue.delete(ii as Idx);
                            vertex_status[ii] = PQ_NOT_PRESENT;
                            if update_index[ii] >= 0 {
                                let pos = update_index[ii] as usize;
                                num_updates -= 1;
                                update_list[pos] = update_list[num_updates];
                                if (update_list[num_updates] as usize) < num_vertices {
                                    update_index[update_list[num_updates] as usize] = pos as Idx;
                                }
                                update_index[ii] = -1;
                            }
                        }
                    } else if vertex_status[ii] == PQ_NOT_PRESENT && ckr_ii.external_degree > 0 {
                        queue.insert(ii as Idx, scaled_gain);
                        vertex_status[ii] = PQ_PRESENT;
                        update_list[num_updates] = ii as Idx;
                        update_index[ii] = num_updates as Idx;
                        num_updates += 1;
                    }
                }
            }
        }

        graph.num_boundary = num_boundary;

        // Reset vertex_status
        for &u in &update_list[..num_updates] {
            let v = u as usize;
            if v < num_vertices {
                vertex_status[v] = PQ_NOT_PRESENT;
                update_index[v] = -1;
            }
        }

        if nmoved == 0 || (omode == MODE_REFINE && graph.edge_cut == oldcut) {
            break;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn update_moved_vertex(
    graph: &mut GraphData,
    ctrl: &mut Control,
    i: usize,
    from: usize,
    k: usize,
    to: usize,
    num_boundary: &mut Idx,
    bndtype: Idx,
) {
    let inbr = graph.kway_refinement_info[i].neighbor_offset as usize;
    let ed_k = ctrl.neighbor_pool[inbr + k].external_degree;

    graph.partition[i] = to as Idx;
    graph.kway_refinement_info[i].external_degree +=
        graph.kway_refinement_info[i].internal_degree - ed_k;

    let tmp = graph.kway_refinement_info[i].internal_degree;
    graph.kway_refinement_info[i].internal_degree = ed_k;
    ctrl.neighbor_pool[inbr + k].external_degree = tmp;

    if ctrl.neighbor_pool[inbr + k].external_degree == 0 {
        let nnbrs = graph.kway_refinement_info[i].num_neighbors as usize;
        if nnbrs > 0 {
            graph.kway_refinement_info[i].num_neighbors -= 1;
            let last = graph.kway_refinement_info[i].num_neighbors as usize;
            ctrl.neighbor_pool[inbr + k] = ctrl.neighbor_pool[inbr + last];
        }
    } else {
        ctrl.neighbor_pool[inbr + k].part_id = from as Idx;
    }

    if bndtype == BOUNDARY_REFINE {
        if graph.boundary_map[i] != -1
            && graph.kway_refinement_info[i].external_degree
                - graph.kway_refinement_info[i].internal_degree
                < 0
        {
            graph.remove_from_boundary(i);
            *num_boundary -= 1;
        }
        if graph.boundary_map[i] == -1
            && graph.kway_refinement_info[i].external_degree
                - graph.kway_refinement_info[i].internal_degree
                >= 0
        {
            graph.add_to_boundary(i);
            *num_boundary += 1;
        }
    } else {
        if graph.boundary_map[i] != -1 && graph.kway_refinement_info[i].external_degree <= 0 {
            graph.remove_from_boundary(i);
            *num_boundary -= 1;
        }
        if graph.boundary_map[i] == -1 && graph.kway_refinement_info[i].external_degree > 0 {
            graph.add_to_boundary(i);
            *num_boundary += 1;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn update_adjacent_vertex(
    ctrl: &mut Control,
    graph: &mut GraphData,
    vid: usize,
    me: usize,
    from: usize,
    to: usize,
    ewgt: Idx,
    num_boundary: &mut Idx,
    bndtype: Idx,
) {
    let adjlen = (graph.xadj[vid + 1] - graph.xadj[vid]) as usize;

    if graph.kway_refinement_info[vid].neighbor_offset == -1 {
        graph.kway_refinement_info[vid].neighbor_offset = ctrl.alloc_neighbor_info(adjlen);
        graph.kway_refinement_info[vid].num_neighbors = 0;
    }

    let inbr = graph.kway_refinement_info[vid].neighbor_offset as usize;

    if me == from {
        graph.kway_refinement_info[vid].external_degree += ewgt;
        graph.kway_refinement_info[vid].internal_degree -= ewgt;
        if bndtype == BOUNDARY_REFINE {
            if graph.kway_refinement_info[vid].external_degree
                - graph.kway_refinement_info[vid].internal_degree
                >= 0
                && graph.boundary_map[vid] == -1
            {
                graph.add_to_boundary(vid);
                *num_boundary += 1;
            }
        } else if graph.kway_refinement_info[vid].external_degree > 0
            && graph.boundary_map[vid] == -1
        {
            graph.add_to_boundary(vid);
            *num_boundary += 1;
        }
    } else if me == to {
        graph.kway_refinement_info[vid].internal_degree += ewgt;
        graph.kway_refinement_info[vid].external_degree -= ewgt;
        if bndtype == BOUNDARY_REFINE {
            if graph.kway_refinement_info[vid].external_degree
                - graph.kway_refinement_info[vid].internal_degree
                < 0
                && graph.boundary_map[vid] != -1
            {
                graph.remove_from_boundary(vid);
                *num_boundary -= 1;
            }
        } else if graph.kway_refinement_info[vid].external_degree <= 0
            && graph.boundary_map[vid] != -1
        {
            graph.remove_from_boundary(vid);
            *num_boundary -= 1;
        }
    }

    // Remove contribution from 'from'
    if me != from {
        let nnbrs = graph.kway_refinement_info[vid].num_neighbors as usize;
        for kk in 0..nnbrs {
            if ctrl.neighbor_pool[inbr + kk].part_id == from as Idx {
                if ctrl.neighbor_pool[inbr + kk].external_degree == ewgt {
                    graph.kway_refinement_info[vid].num_neighbors -= 1;
                    let last = graph.kway_refinement_info[vid].num_neighbors as usize;
                    ctrl.neighbor_pool[inbr + kk] = ctrl.neighbor_pool[inbr + last];
                } else {
                    ctrl.neighbor_pool[inbr + kk].external_degree -= ewgt;
                }
                break;
            }
        }
    }

    // Add contribution to 'to'
    if me != to {
        let nnbrs = graph.kway_refinement_info[vid].num_neighbors as usize;
        let mut found = false;
        for kk in 0..nnbrs {
            if ctrl.neighbor_pool[inbr + kk].part_id == to as Idx {
                ctrl.neighbor_pool[inbr + kk].external_degree += ewgt;
                found = true;
                break;
            }
        }
        if !found {
            let kk = graph.kway_refinement_info[vid].num_neighbors as usize;
            ctrl.neighbor_pool[inbr + kk].part_id = to as Idx;
            ctrl.neighbor_pool[inbr + kk].external_degree = ewgt;
            graph.kway_refinement_info[vid].num_neighbors += 1;
        }
    }
}

/// Compute k-way PQ gain matching C METIS precision.
/// C METIS computes: rgain = (nnbrs > 0 ? 1.0*ed/sqrt(nnbrs) : 0.0) - id
/// The computation is done in double but stored in real_t (f32) before insertion
/// into the RPQ. We replicate this by truncating to f32 precision.
fn compute_kway_gain(ed: Idx, nnbrs: Idx, id: Idx) -> f64 {
    let gain = if nnbrs > 0 {
        (ed as f64) / (nnbrs as f64).sqrt() - id as f64
    } else {
        -(id as f64)
    };
    // C METIS stores gain as real_t (f32) in the RPQ. Truncate to match.
    gain as f32 as f64
}

fn compute_load_imbalance_kway(
    graph: &GraphData,
    nparts: usize,
    partition_ij_balance_multipliers: &[Real],
) -> Real {
    let mut max_imbal: Real = 0.0;
    for i in 0..nparts {
        if i < graph.part_weights.len() && i < partition_ij_balance_multipliers.len() {
            let imbal = graph.part_weights[i] as Real * partition_ij_balance_multipliers[i];
            if imbal > max_imbal {
                max_imbal = imbal;
            }
        }
    }
    max_imbal
}

// ======= Multi-constraint k-way helpers =======

/// Returns true if a*x[i] + y[i] <= z[i] for all i in 0..n.
fn ivecaxpylez(n: usize, a: Idx, x: &[Idx], y: &[Idx], z: &[Idx]) -> bool {
    for i in (0..n).rev() {
        if a * x[i] + y[i] > z[i] {
            return false;
        }
    }
    true
}

/// Returns true if a*x[i] + y[i] >= z[i] for all i in 0..n.
fn ivecaxpygez(n: usize, a: Idx, x: &[Idx], y: &[Idx], z: &[Idx]) -> bool {
    for i in (0..n).rev() {
        if a * x[i] + y[i] < z[i] {
            return false;
        }
    }
    true
}

/// BetterBalanceKWay: returns true if option 2 (a2,pt2,bm2) has better balance than option 1 (a1,pt1,bm1).
/// Matches C METIS BetterBalanceKWay from mcutil.c.
#[allow(clippy::too_many_arguments)]
fn better_balance_kway(
    ncon: usize,
    vwgt: &[Idx],
    ubvec: &[Real],
    a1: Idx,
    pt1: &[Idx],
    bm1: &[Real],
    a2: Idx,
    pt2: &[Idx],
    bm2: &[Real],
) -> bool {
    let mut nrm1: Real = 0.0;
    let mut nrm2: Real = 0.0;
    let mut max1: Real = 0.0;
    let mut max2: Real = 0.0;

    for i in 0..ncon {
        let tmp1 = bm1[i] * (pt1[i] + a1 * vwgt[i]) as Real - ubvec[i];
        nrm1 += tmp1 * tmp1;
        if tmp1 > max1 {
            max1 = tmp1;
        }

        let tmp2 = bm2[i] * (pt2[i] + a2 * vwgt[i]) as Real - ubvec[i];
        nrm2 += tmp2 * tmp2;
        if tmp2 > max2 {
            max2 = tmp2;
        }
    }

    if max2 < max1 {
        return true;
    }
    if max2 == max1 && nrm2 < nrm1 {
        return true;
    }
    false
}

/// ComputeLoadImbalanceVec: per-constraint max load imbalance across all partitions.
fn compute_load_imbalance_vec_kway(graph: &GraphData, nparts: usize, pijbm: &[Real]) -> Vec<Real> {
    let ncon = graph.num_constraints as usize;
    let mut lbvec = vec![0.0 as Real; ncon];
    for i in 0..ncon {
        lbvec[i] = graph.part_weights[i] as Real * pijbm[i];
        for j in 1..nparts {
            let cur = graph.part_weights[j * ncon + i] as Real * pijbm[j * ncon + i];
            if cur > lbvec[i] {
                lbvec[i] = cur;
            }
        }
    }
    lbvec
}

/// ComputeLoadImbalanceDiff for multi-constraint k-way.
fn compute_load_imbalance_diff_mc(
    graph: &GraphData,
    nparts: usize,
    pijbm: &[Real],
    ubfactors: &[Real],
) -> Real {
    let ncon = graph.num_constraints as usize;
    let mut max_diff = Real::NEG_INFINITY;
    for i in 0..nparts {
        for j in 0..ncon {
            let idx = i * ncon + j;
            if idx < graph.part_weights.len() && idx < pijbm.len() && j < ubfactors.len() {
                let diff = graph.part_weights[idx] as Real * pijbm[idx] - ubfactors[j];
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }
    max_diff
}

/// Greedy multi-constraint k-way cut optimization matching METIS Greedy_McKWayCutOptimize.
#[allow(unused_assignments)]
fn greedy_mc_kway_cut_optimize(ctrl: &mut Control, graph: &mut GraphData, niter: Idx, omode: Idx) {
    let num_vertices = graph.num_vertices as usize;
    let nparts = ctrl.num_parts as usize;
    let ncon = graph.num_constraints as usize;
    let bndtype = if omode == MODE_REFINE {
        BOUNDARY_REFINE
    } else {
        BOUNDARY_BALANCE
    };

    if num_vertices == 0 {
        return;
    }

    // Clone pijbm to avoid borrow conflicts with ctrl mutation
    let pijbm = ctrl.partition_ij_balance_multipliers.clone();

    // Compute per-constraint ubfactors
    let mut ubfactors = compute_load_imbalance_vec_kway(graph, nparts, &pijbm);
    if omode == MODE_BALANCE {
        ubfactors = ctrl.imbalance_tols[..ncon].to_vec();
    } else {
        #[allow(clippy::needless_range_loop)]
        for i in 0..ncon {
            if ctrl.imbalance_tols[i] > ubfactors[i] {
                ubfactors[i] = ctrl.imbalance_tols[i];
            }
        }
    }

    // Setup weight intervals [nparts * ncon]
    let mut minpwgts = vec![0 as Idx; nparts * ncon];
    let mut maxpwgts = vec![0 as Idx; nparts * ncon];
    for i in 0..nparts {
        #[allow(clippy::needless_range_loop)]
        for j in 0..ncon {
            let idx = i * ncon + j;
            maxpwgts[idx] = (ctrl.target_part_weights[idx]
                * graph.total_vertex_weight[j] as Real
                * ubfactors[j]) as Idx;
            minpwgts[idx] =
                (ctrl.target_part_weights[idx] * graph.total_vertex_weight[j] as Real * 0.2) as Idx;
        }
    }

    let mut perm = vec![0 as Idx; num_vertices];
    let mut vertex_status = vec![PQ_NOT_PRESENT; num_vertices];
    let mut update_index = vec![-1 as Idx; num_vertices];
    let mut update_list = vec![0 as Idx; num_vertices];

    let mut queue = PQueue::new(num_vertices);

    for _pass in 0..niter {
        // Balance mode: exit as soon as balance is reached (IsBalanced with ffactor=0)
        if omode == MODE_BALANCE
            && compute_load_imbalance_diff_mc(graph, nparts, &pijbm, &ctrl.imbalance_tols) <= 0.0
        {
            break;
        }

        let oldcut = graph.edge_cut;
        let mut num_boundary = graph.num_boundary;
        let mut num_updates: usize = 0;

        // Insert boundary vertices in PQ
        let num_boundary_usize = num_boundary as usize;
        #[allow(clippy::needless_range_loop)]
        for i in 0..num_boundary_usize {
            perm[i] = i as Idx;
        }
        if num_boundary_usize > 0 {
            let nshuffles = num_boundary_usize / 4;
            if num_boundary_usize < 10 {
                ctrl.rng
                    .rand_array_permute(num_boundary_usize, &mut perm, 0, true);
            } else {
                ctrl.rng.rand_array_permute_with_nshuffles(
                    num_boundary_usize,
                    &mut perm,
                    0,
                    nshuffles,
                    true,
                );
            }
        }

        for &p in &perm[..num_boundary_usize] {
            let i = graph.boundary_list[p as usize] as usize;
            if i >= num_vertices {
                continue;
            }
            let ckr = &graph.kway_refinement_info[i];
            let scaled_gain =
                compute_kway_gain(ckr.external_degree, ckr.num_neighbors, ckr.internal_degree);
            queue.insert(i as Idx, scaled_gain);
            vertex_status[i] = PQ_PRESENT;
            update_list[num_updates] = i as Idx;
            update_index[i] = num_updates as Idx;
            num_updates += 1;
        }

        let mut nmoved = 0;
        let mut iii: usize = 0; // vertex extraction counter (matches C for-loop iii)
        while let Some((v, _)) = queue.get_top() {
            let i = v as usize;
            vertex_status[i] = PQ_EXTRACTED;

            let from = graph.partition[i] as usize;

            let id_i = graph.kway_refinement_info[i].internal_degree;
            let nnbrs_i = graph.kway_refinement_info[i].num_neighbors;
            let inbr = graph.kway_refinement_info[i].neighbor_offset;
            if inbr < 0 || nnbrs_i == 0 {
                iii += 1;
                continue;
            }

            // Copy vertex weights for borrow safety
            let vw: Vec<Idx> = graph.vertex_weights[i * ncon..(i + 1) * ncon].to_vec();
            let from_start = from * ncon;

            // Pre-move underbalance check
            if omode == MODE_REFINE {
                if id_i > 0
                    && !ivecaxpygez(
                        ncon,
                        -1,
                        &vw,
                        &graph.part_weights[from_start..from_start + ncon],
                        &minpwgts[from_start..from_start + ncon],
                    )
                {
                    iii += 1;
                    continue;
                }
            } else if !ivecaxpygez(
                ncon,
                -1,
                &vw,
                &graph.part_weights[from_start..from_start + ncon],
                &minpwgts[from_start..from_start + ncon],
            ) {
                iii += 1;
                continue;
            }

            let mut best_k: i32 = -1;
            let mut final_to: usize = 0;

            if omode == MODE_REFINE {
                // Find first valid candidate from the end
                for kk in (0..nnbrs_i as usize).rev() {
                    let to = ctrl.neighbor_pool[inbr as usize + kk].part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let gain = ctrl.neighbor_pool[inbr as usize + kk].external_degree - id_i;
                    let to_start = to * ncon;
                    if gain >= 0
                        && ivecaxpylez(
                            ncon,
                            1,
                            &vw,
                            &graph.part_weights[to_start..to_start + ncon],
                            &maxpwgts[to_start..to_start + ncon],
                        )
                    {
                        best_k = kk as i32;
                        break;
                    }
                }
                if best_k < 0 {
                    iii += 1;
                    continue;
                }

                // Compare with earlier candidates for better destination
                let mut cto = ctrl.neighbor_pool[inbr as usize + best_k as usize].part_id as usize;
                for jj in (0..best_k as usize).rev() {
                    let to = ctrl.neighbor_pool[inbr as usize + jj].part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let to_start = to * ncon;
                    let ed_j = ctrl.neighbor_pool[inbr as usize + jj].external_degree;
                    let ed_k = ctrl.neighbor_pool[inbr as usize + best_k as usize].external_degree;
                    if (ed_j > ed_k
                        && ivecaxpylez(
                            ncon,
                            1,
                            &vw,
                            &graph.part_weights[to_start..to_start + ncon],
                            &maxpwgts[to_start..to_start + ncon],
                        ))
                        || (ed_j == ed_k
                            && better_balance_kway(
                                ncon,
                                &vw,
                                &ubfactors,
                                1,
                                &graph.part_weights[cto * ncon..(cto + 1) * ncon],
                                &pijbm[cto * ncon..(cto + 1) * ncon],
                                1,
                                &graph.part_weights[to_start..to_start + ncon],
                                &pijbm[to_start..to_start + ncon],
                            ))
                    {
                        best_k = jj as i32;
                        cto = to;
                    }
                }
                final_to = cto;

                // Final acceptance check
                let gain =
                    ctrl.neighbor_pool[inbr as usize + best_k as usize].external_degree - id_i;
                let to_start = final_to * ncon;
                if !(gain > 0
                    || (gain == 0
                        && (
                            better_balance_kway(
                                ncon,
                                &vw,
                                &ubfactors,
                                -1,
                                &graph.part_weights[from_start..from_start + ncon],
                                &pijbm[from_start..from_start + ncon],
                                1,
                                &graph.part_weights[to_start..to_start + ncon],
                                &pijbm[to_start..to_start + ncon],
                            ) || iii.is_multiple_of(2)
                            // safetos always 2 when no minconn
                        )))
                {
                    iii += 1;
                    continue;
                }
            } else {
                // MODE_BALANCE
                for kk in (0..nnbrs_i as usize).rev() {
                    let to = ctrl.neighbor_pool[inbr as usize + kk].part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let to_start = to * ncon;
                    if ivecaxpylez(
                        ncon,
                        1,
                        &vw,
                        &graph.part_weights[to_start..to_start + ncon],
                        &maxpwgts[to_start..to_start + ncon],
                    ) || better_balance_kway(
                        ncon,
                        &vw,
                        &ubfactors,
                        -1,
                        &graph.part_weights[from_start..from_start + ncon],
                        &pijbm[from_start..from_start + ncon],
                        1,
                        &graph.part_weights[to_start..to_start + ncon],
                        &pijbm[to_start..to_start + ncon],
                    ) {
                        best_k = kk as i32;
                        break;
                    }
                }
                if best_k < 0 {
                    iii += 1;
                    continue;
                }

                let mut cto = ctrl.neighbor_pool[inbr as usize + best_k as usize].part_id as usize;
                for jj in (0..best_k as usize).rev() {
                    let to = ctrl.neighbor_pool[inbr as usize + jj].part_id as usize;
                    if to >= nparts {
                        continue;
                    }
                    let to_start = to * ncon;
                    if better_balance_kway(
                        ncon,
                        &vw,
                        &ubfactors,
                        1,
                        &graph.part_weights[cto * ncon..(cto + 1) * ncon],
                        &pijbm[cto * ncon..(cto + 1) * ncon],
                        1,
                        &graph.part_weights[to_start..to_start + ncon],
                        &pijbm[to_start..to_start + ncon],
                    ) {
                        best_k = jj as i32;
                        cto = to;
                    }
                }
                final_to = cto;

                // Final acceptance: if gain < 0, require BetterBalance
                let gain =
                    ctrl.neighbor_pool[inbr as usize + best_k as usize].external_degree - id_i;
                let to_start = final_to * ncon;
                if gain < 0
                    && !better_balance_kway(
                        ncon,
                        &vw,
                        &ubfactors,
                        -1,
                        &graph.part_weights[from_start..from_start + ncon],
                        &pijbm[from_start..from_start + ncon],
                        1,
                        &graph.part_weights[to_start..to_start + ncon],
                        &pijbm[to_start..to_start + ncon],
                    )
                {
                    iii += 1;
                    continue;
                }
            }

            // Perform the move
            let k = best_k as usize;
            let to = final_to;
            let to_start = to * ncon;
            let ed_k = ctrl.neighbor_pool[inbr as usize + k].external_degree;

            let gain = ed_k - graph.kway_refinement_info[i].internal_degree;
            graph.edge_cut -= gain;
            nmoved += 1;

            // Update part weights (multi-constraint: ncon values)
            #[allow(clippy::needless_range_loop)]
            for j in 0..ncon {
                graph.part_weights[to_start + j] += vw[j];
                graph.part_weights[from_start + j] -= vw[j];
            }

            // UpdateMovedVertexInfoAndBND
            update_moved_vertex(graph, ctrl, i, from, k, to, &mut num_boundary, bndtype);

            // Update adjacent vertices
            for jj in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let ii = graph.adjacency[jj] as usize;
                let me = graph.partition[ii] as usize;
                let ewgt = graph.edge_weights[jj];
                let old_nnbrs = graph.kway_refinement_info[ii].num_neighbors;

                update_adjacent_vertex(
                    ctrl,
                    graph,
                    ii,
                    me,
                    from,
                    to,
                    ewgt,
                    &mut num_boundary,
                    bndtype,
                );

                let new_nnbrs = graph.kway_refinement_info[ii].num_neighbors;
                if me == to || me == from || old_nnbrs != new_nnbrs {
                    let ckr_ii = &graph.kway_refinement_info[ii];
                    let scaled_gain = compute_kway_gain(
                        ckr_ii.external_degree,
                        ckr_ii.num_neighbors,
                        ckr_ii.internal_degree,
                    );

                    if bndtype == BOUNDARY_REFINE {
                        if vertex_status[ii] == PQ_PRESENT {
                            if ckr_ii.external_degree - ckr_ii.internal_degree >= 0 {
                                queue.update(ii as Idx, scaled_gain);
                            } else {
                                queue.delete(ii as Idx);
                                vertex_status[ii] = PQ_NOT_PRESENT;
                                if update_index[ii] >= 0 {
                                    let pos = update_index[ii] as usize;
                                    num_updates -= 1;
                                    update_list[pos] = update_list[num_updates];
                                    if (update_list[num_updates] as usize) < num_vertices {
                                        update_index[update_list[num_updates] as usize] =
                                            pos as Idx;
                                    }
                                    update_index[ii] = -1;
                                }
                            }
                        } else if vertex_status[ii] == PQ_NOT_PRESENT
                            && ckr_ii.external_degree - ckr_ii.internal_degree >= 0
                        {
                            queue.insert(ii as Idx, scaled_gain);
                            vertex_status[ii] = PQ_PRESENT;
                            update_list[num_updates] = ii as Idx;
                            update_index[ii] = num_updates as Idx;
                            num_updates += 1;
                        }
                    } else if vertex_status[ii] == PQ_PRESENT {
                        if ckr_ii.external_degree > 0 {
                            queue.update(ii as Idx, scaled_gain);
                        } else {
                            queue.delete(ii as Idx);
                            vertex_status[ii] = PQ_NOT_PRESENT;
                            if update_index[ii] >= 0 {
                                let pos = update_index[ii] as usize;
                                num_updates -= 1;
                                update_list[pos] = update_list[num_updates];
                                if (update_list[num_updates] as usize) < num_vertices {
                                    update_index[update_list[num_updates] as usize] = pos as Idx;
                                }
                                update_index[ii] = -1;
                            }
                        }
                    } else if vertex_status[ii] == PQ_NOT_PRESENT && ckr_ii.external_degree > 0 {
                        queue.insert(ii as Idx, scaled_gain);
                        vertex_status[ii] = PQ_PRESENT;
                        update_list[num_updates] = ii as Idx;
                        update_index[ii] = num_updates as Idx;
                        num_updates += 1;
                    }
                }
            }

            iii += 1;
        }

        graph.num_boundary = num_boundary;

        // Reset vertex_status
        for &u in &update_list[..num_updates] {
            let v = u as usize;
            if v < num_vertices {
                vertex_status[v] = PQ_NOT_PRESENT;
                update_index[v] = -1;
            }
        }

        if nmoved == 0 || (omode == MODE_REFINE && graph.edge_cut == oldcut) {
            break;
        }
    }
}
