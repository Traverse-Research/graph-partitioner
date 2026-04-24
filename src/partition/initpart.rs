use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;

/// Initialize a 2-way partition on the coarsest graph. Tries niparts attempts.
///
/// Matches C METIS Init2WayPartition from initpart.c:
/// - iptype==RANDOM && ncon==1: RandomBisection
/// - iptype==RANDOM && ncon>1: McRandomBisection
/// - iptype==GROW && ncon==1: GrowBisection (fallback to Random if no edges)
/// - iptype==GROW && ncon>1: McGrowBisection (fallback to McRandom if no edges)
pub fn init_2way_partition(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    graph.alloc_2way();

    let ncon = graph.num_constraints as usize;

    match ctrl.init_part_type {
        1 => {
            // METIS_IPTYPE_RANDOM
            if ncon == 1 {
                random_bisection_loop(ctrl, graph, target_part_weights, niparts);
            } else {
                mc_random_bisection(ctrl, graph, target_part_weights, niparts);
            }
        }
        _ => {
            // METIS_IPTYPE_GROW (or default)
            if graph.num_edges == 0 {
                if ncon == 1 {
                    random_bisection_loop(ctrl, graph, target_part_weights, niparts);
                } else {
                    mc_random_bisection(ctrl, graph, target_part_weights, niparts);
                }
            } else {
                if ncon == 1 {
                    grow_bisection_loop(ctrl, graph, target_part_weights, niparts);
                } else {
                    mc_grow_bisection(ctrl, graph, target_part_weights, niparts);
                }
            }
        }
    }
}

/// RandomBisection loop: for ncon==1, try niparts random bisections.
fn random_bisection_loop(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; graph.num_vertices as usize];

    for inbfs in 0..niparts {
        random_bisection(ctrl, graph, target_part_weights);

        compute_2way_partition_params(ctrl, graph);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, 4);

        if inbfs == 0 || bestcut > graph.edge_cut {
            bestcut = graph.edge_cut;
            bestwhere.copy_from_slice(&graph.partition);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.edge_cut = bestcut;
    graph.partition.copy_from_slice(&bestwhere);
}

/// GrowBisection loop: for ncon==1, try niparts grow bisections.
fn grow_bisection_loop(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; graph.num_vertices as usize];

    for inbfs in 0..niparts {
        grow_bisection(ctrl, graph, target_part_weights);

        compute_2way_partition_params(ctrl, graph);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);

        if inbfs == 0 || bestcut > graph.edge_cut {
            bestcut = graph.edge_cut;
            bestwhere.copy_from_slice(&graph.partition);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.edge_cut = bestcut;
    graph.partition.copy_from_slice(&bestwhere);
}

/// McRandomBisection: multi-constraint random bisection.
///
/// Matches C METIS McRandomBisection from initpart.c:
/// - Loops 2*niparts times
/// - Each iteration: permute vertices, assign by round-robin on dominant constraint
/// - Compute2WayPartitionParams
/// - FM_2WayRefine → Balance2Way → FM → Balance → FM (3 FM + 2 Balance)
/// - Best cut: bestcut >= graph->mincut (note >=, not >)
fn mc_random_bisection(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; num_vertices];
    let mut perm = vec![0 as Idx; num_vertices];
    let mut counts = vec![0 as Idx; ncon];

    for inbfs in 0..(2 * niparts) {
        // Permute vertices
        ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut perm, 0, num_vertices / 2, true);

        // Reset counts
        for c in counts.iter_mut() { *c = 0; }

        // Assign vertices by round-robin on dominant constraint
        for ii in 0..num_vertices {
            let i = perm[ii] as usize;
            // iargmax(ncon, vwgt+i*ncon, 1) - find constraint with maximum weight
            let mut qnum = 0;
            for j in 1..ncon {
                if graph.vertex_weights[i * ncon + j] > graph.vertex_weights[i * ncon + qnum] {
                    qnum = j;
                }
            }
            graph.partition[i] = counts[qnum] % 2;
            counts[qnum] += 1;
        }

        compute_2way_partition_params(ctrl, graph);
        // FM → Balance → FM → Balance → FM
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);

        // Note: C METIS uses >= for mc, not > like single-constraint
        if inbfs == 0 || bestcut >= graph.edge_cut {
            bestcut = graph.edge_cut;
            bestwhere.copy_from_slice(&graph.partition);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.edge_cut = bestcut;
    graph.partition.copy_from_slice(&bestwhere);
}

/// McGrowBisection: multi-constraint grow bisection.
///
/// Matches C METIS McGrowBisection from initpart.c:
/// - Loops 2*niparts times
/// - Each iteration: set all to partition 1, randomly pick one vertex to partition 0
/// - Compute2WayPartitionParams
/// - Balance2Way → FM → Balance → FM (2 FM + 2 Balance)
/// - Best cut: bestcut >= graph->mincut
fn mc_grow_bisection(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    let num_vertices = graph.num_vertices as usize;

    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; num_vertices];

    for inbfs in 0..(2 * niparts) {
        // Set all to partition 1
        for i in 0..num_vertices {
            graph.partition[i] = 1;
        }
        // Pick one random vertex for partition 0
        let seed_v = ctrl.rng.rand_in_range(graph.num_vertices) as usize;
        graph.partition[seed_v] = 0;

        compute_2way_partition_params(ctrl, graph);

        // Balance → FM → Balance → FM
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_refine(ctrl, graph, target_part_weights, ctrl.num_iter);

        if inbfs == 0 || bestcut >= graph.edge_cut {
            bestcut = graph.edge_cut;
            bestwhere.copy_from_slice(&graph.partition);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.edge_cut = bestcut;
    graph.partition.copy_from_slice(&bestwhere);
}

/// GrowBisection: BFS from random seed, growing partition 0.
fn grow_bisection(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Target weights for partition 1
    let max_target_weight = (ctrl.imbalance_tols[0] * graph.total_vertex_weight[0] as Real * target_part_weights[ncon]) as Idx;
    let min_target_weight = ((1.0 / ctrl.imbalance_tols[0]) * graph.total_vertex_weight[0] as Real * target_part_weights[ncon]) as Idx;

    // Initialize all to partition 1
    for i in 0..num_vertices {
        graph.partition[i] = 1;
    }

    let mut touched = vec![0u8; num_vertices];
    let mut queue = vec![0usize; num_vertices];
    let mut part_weights = [0 as Idx; 2];
    part_weights[1] = graph.total_vertex_weight[0];

    // Pick random seed
    let start = ctrl.rng.rand_in_range(graph.num_vertices) as usize;
    queue[0] = start;
    touched[start] = 1;
    let mut first = 0usize;
    let mut last = 1usize;
    let mut nleft = num_vertices - 1;
    let mut drain = false;

    loop {
        if first == last {
            if nleft == 0 || drain {
                break;
            }
            // Find random untouched vertex
            let k = ctrl.rng.rand_in_range(nleft as Idx) as usize;
            let mut cnt = 0;
            for i in 0..num_vertices {
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
        if part_weights[0] > 0 && part_weights[1] - graph.vertex_weights[i] < min_target_weight {
            drain = true;
            continue;
        }

        // Move to partition 0
        graph.partition[i] = 0;
        part_weights[0] += graph.vertex_weights[i];
        part_weights[1] -= graph.vertex_weights[i];

        if part_weights[1] <= max_target_weight {
            break;
        }

        drain = false;

        // Expand BFS
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let j = graph.adjacency[k] as usize;
            if touched[j] == 0 {
                touched[j] = 1;
                queue[last] = j;
                last += 1;
                nleft -= 1;
            }
        }
    }

    // Degenerate cases
    if part_weights[1] == 0 {
        let v = ctrl.rng.rand_in_range(graph.num_vertices) as usize;
        graph.partition[v] = 1;
    }
    if part_weights[0] == 0 {
        let v = ctrl.rng.rand_in_range(graph.num_vertices) as usize;
        graph.partition[v] = 0;
    }
}

/// RandomBisection: randomly assign vertices to partition 0.
fn random_bisection(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real]) {
    let num_vertices = graph.num_vertices as usize;

    let max_zero_weight = (ctrl.imbalance_tols[0] * graph.total_vertex_weight[0] as Real * target_part_weights[0]) as Idx;

    for i in 0..num_vertices {
        graph.partition[i] = 1;
    }

    if num_vertices < 10 {
        let mut perm = vec![0 as Idx; num_vertices];
        ctrl.rng.rand_array_permute(num_vertices, &mut perm, 0, true);

        let mut pwgt0 = 0 as Idx;
        for &pi in &perm {
            let i = pi as usize;
            if pwgt0 + graph.vertex_weights[i] <= max_zero_weight {
                graph.partition[i] = 0;
                pwgt0 += graph.vertex_weights[i];
            }
        }
    } else {
        let mut perm = vec![0 as Idx; num_vertices];
        let nshuffles = num_vertices / 2;
        ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut perm, 0, nshuffles, true);

        let mut pwgt0 = 0 as Idx;
        for &pi in &perm {
            let i = pi as usize;
            if pwgt0 + graph.vertex_weights[i] <= max_zero_weight {
                graph.partition[i] = 0;
                pwgt0 += graph.vertex_weights[i];
            }
        }
    }
}

/// Compute 2-way partition parameters from scratch.
pub fn compute_2way_partition_params(_ctrl: &Control, graph: &mut GraphData) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Reset
    graph.part_weights = vec![0; 2 * ncon];
    graph.boundary_map = vec![-1; num_vertices];
    graph.boundary_list = vec![0; num_vertices];
    graph.num_boundary = 0;
    graph.internal_degree = vec![0; num_vertices];
    graph.external_degree = vec![0; num_vertices];

    // Compute part_weights
    if ncon == 1 {
        for i in 0..num_vertices {
            graph.part_weights[graph.partition[i] as usize] += graph.vertex_weights[i];
        }
    } else {
        for i in 0..num_vertices {
            let p = graph.partition[i] as usize;
            for j in 0..ncon {
                graph.part_weights[p * ncon + j] += graph.vertex_weights[i * ncon + j];
            }
        }
    }

    // Compute id, ed, boundary, edge_cut
    let mut edge_cut: Idx = 0;

    for i in 0..num_vertices {
        let me = graph.partition[i];
        let mut tid: Idx = 0;
        let mut ted: Idx = 0;

        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;

        for k in istart..iend {
            if graph.partition[graph.adjacency[k] as usize] == me {
                tid += graph.edge_weights[k];
            } else {
                ted += graph.edge_weights[k];
            }
        }

        graph.internal_degree[i] = tid;
        graph.external_degree[i] = ted;

        if ted > 0 || istart == iend {
            graph.add_to_boundary(i);
            edge_cut += ted;
        }
    }

    graph.edge_cut = edge_cut / 2;
}
