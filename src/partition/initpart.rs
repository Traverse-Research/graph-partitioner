use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::graph::GraphData;

/// Initialize a 2-way partition on the coarsest graph. Tries niparts attempts.
pub fn init_2way_partition(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real], niparts: Idx) {
    graph.alloc_2way();
    let mut bestcut = Idx::MAX;
    let mut bestwhere = vec![0 as Idx; graph.num_vertices as usize];

    for _inbfs in 0..niparts {
        if ctrl.init_part_type == 1 {
            random_bisection(ctrl, graph, target_part_weights);
        } else {
            grow_bisection(ctrl, graph, target_part_weights);
        }

        compute_2way_partition_params(ctrl, graph);
        super::balance::balance_2way(ctrl, graph, target_part_weights);
        super::fm::fm_2way_cut_refine(ctrl, graph, target_part_weights, ctrl.num_iter);

        if graph.edge_cut < bestcut {
            bestcut = graph.edge_cut;
            bestwhere.copy_from_slice(&graph.partition);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.partition.copy_from_slice(&bestwhere);
    compute_2way_partition_params(ctrl, graph);
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
