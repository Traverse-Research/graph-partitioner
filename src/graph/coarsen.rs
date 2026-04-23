use crate::types::Idx;
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::graph::contract;

const COARSEN_FRACTION: f64 = 0.85;
const UNMATCHED: Idx = -1;
const UNMATCHED_THRESHOLD_2HOP: f64 = 0.10;

/// Coarsen the graph iteratively. Returns the coarser levels as a Vec.
///
/// `graph` is the finest level. Its `coarse_map` is set by the matching step
/// to map fine vertices to the first coarser level (`levels[0]`).
/// Each `levels[i].coarse_map` maps to `levels[i+1]`.
/// The coarsest level is `levels[levels.len() - 1]`.
pub fn coarsen_graph(ctrl: &mut Control, graph: &mut GraphData) -> Vec<GraphData> {
    let ncon = graph.num_constraints as usize;

    // Check if all edge weights are equal
    let eqewgts = if graph.num_edges > 0 {
        let first = graph.edge_weights[0];
        graph.edge_weights.iter().all(|&w| w == first)
    } else {
        true
    };

    // Compute max_vertex_weight: 1.5 * total_vertex_weight[i] / CoarsenTo
    for i in 0..ncon {
        ctrl.max_vertex_weight[i] = ((1.5 * graph.total_vertex_weight[i] as f64) / ctrl.coarsen_to.max(1) as f64) as Idx;
    }

    let mut levels: Vec<GraphData> = Vec::new();
    let mut cur_eqewgts = eqewgts;

    // First iteration: match and coarsen the finest graph
    {
        let num_vertices = graph.num_vertices as usize;
        graph.coarse_map = vec![0; num_vertices];

        if ctrl.coarsen_type == 1 && !cur_eqewgts && graph.num_edges > 0 {
            matchingshem(ctrl, graph);
        } else {
            matchingrm(ctrl, graph);
        }

        let cnum_vertices = graph.coarse_map.iter().copied().max().unwrap_or(-1) + 1;
        let coarse = contract::create_coarse_graph(ctrl, graph, cnum_vertices);

        let continue_coarsening =
            coarse.num_vertices > ctrl.coarsen_to
            && (coarse.num_vertices as f64) < COARSEN_FRACTION * graph.num_vertices as f64
            && coarse.num_edges > coarse.num_vertices / 2;

        levels.push(coarse);

        if !continue_coarsening {
            return levels;
        }
        cur_eqewgts = false;
    }

    // Subsequent iterations: coarsen the last level in the arena
    loop {
        let g = levels.last_mut().unwrap();
        let num_vertices = g.num_vertices as usize;
        g.coarse_map = vec![0; num_vertices];

        if ctrl.coarsen_type == 1 && !cur_eqewgts && g.num_edges > 0 {
            matchingshem(ctrl, g);
        } else {
            matchingrm(ctrl, g);
        }

        let cnum_vertices = g.coarse_map.iter().copied().max().unwrap_or(-1) + 1;
        let coarse = contract::create_coarse_graph(ctrl, g, cnum_vertices);

        let continue_coarsening =
            coarse.num_vertices > ctrl.coarsen_to
            && (coarse.num_vertices as f64) < COARSEN_FRACTION * g.num_vertices as f64
            && coarse.num_edges > coarse.num_vertices / 2;

        levels.push(coarse);

        if !continue_coarsening {
            break;
        }
    }

    levels
}

/// Sorted Heavy Edge Matching (SHEM).
fn matchingshem(ctrl: &mut Control, graph: &mut GraphData) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    let mut matching = vec![UNMATCHED; num_vertices];
    let mut perm = vec![0 as Idx; num_vertices];
    let mut tperm = vec![0 as Idx; num_vertices];

    // Create degree-biased random permutation
    // Step 1: random permutation of tperm
    let nshuffles = num_vertices / 8;
    ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut tperm, 0, nshuffles, true);

    // Step 2: degree-based bucket sort
    // C METIS: avgdegree = 4.0*(xadj[num_vertices]/num_vertices) — integer division first
    let avgdegree = if num_vertices > 0 && graph.num_edges > 0 {
        (4.0 * (graph.num_edges / graph.num_vertices) as f64) as Idx
    } else {
        1
    };
    let avgdegree = avgdegree.max(1);

    let mut degrees = vec![0 as Idx; num_vertices];
    for i in 0..num_vertices {
        let deg = graph.xadj[i + 1] - graph.xadj[i];
        let bnum = ((1.0 + deg as f64).sqrt()) as Idx;
        degrees[i] = bnum.min(avgdegree);
    }

    // Bucket sort tperm by degrees (stable, increasing)
    bucket_sort_perm(num_vertices, avgdegree, &degrees, &tperm, &mut perm);

    let mut cnum_vertices: Idx = 0;
    let mut last_unmatched = 0usize;
    let mut nunmatched = 0;

    for pi in 0..num_vertices {
        let i = perm[pi] as usize;
        if matching[i] != UNMATCHED {
            continue;
        }

        let mut maxidx: i64 = i as i64;

        // Weight guard: skip if already at max_vertex_weight
        let mut at_max = false;
        if ncon == 1 {
            if graph.vertex_weights[i] >= ctrl.max_vertex_weight[0] {
                at_max = true;
            }
        } else {
            at_max = (0..ncon).all(|j| graph.vertex_weights[i * ncon + j] >= ctrl.max_vertex_weight[j]);
        }

        if at_max {
            // Self-match
        } else if graph.xadj[i] == graph.xadj[i + 1] {
            // Island vertex: find another unmatched vertex
            maxidx = find_island_partner(&perm, &matching, &mut last_unmatched, num_vertices, i) as i64;
        } else {
            // SHEM: find heaviest unmatched neighbor
            let mut maxwgt: Idx = -1;
            maxidx = -1; // Will be set to partner or UNMATCHED for 2-hop candidate

            if ncon == 1 {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjacency[k] as usize;
                    if matching[j] == UNMATCHED
                        && maxwgt < graph.edge_weights[k]
                        && graph.vertex_weights[i] + graph.vertex_weights[j] <= ctrl.max_vertex_weight[0]
                    {
                        maxidx = j as i64;
                        maxwgt = graph.edge_weights[k];
                    }
                }
            } else {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjacency[k] as usize;
                    if matching[j] == UNMATCHED
                        && ivec_le_sum(ncon, &graph.vertex_weights, i, j, &ctrl.max_vertex_weight)
                        && maxwgt < graph.edge_weights[k]
                    {
                        maxidx = j as i64;
                        maxwgt = graph.edge_weights[k];
                    }
                }
            }

            if maxidx == -1 {
                // Check if 2-hop candidate
                if ncon == 1 {
                    if 2 * graph.vertex_weights[i] < ctrl.max_vertex_weight[0] {
                        nunmatched += 1;
                        maxidx = UNMATCHED as i64;
                    } else {
                        maxidx = i as i64;
                    }
                } else {
                    maxidx = i as i64;
                }
            }
        }

        if maxidx != UNMATCHED as i64 {
            let maxidx = if maxidx < 0 { i } else { maxidx as usize };
            graph.coarse_map[i] = cnum_vertices;
            graph.coarse_map[maxidx] = cnum_vertices;
            matching[i] = maxidx as Idx;
            matching[maxidx] = i as Idx;
            cnum_vertices += 1;
        }
    }

    // 2-hop matching
    if !ctrl.disable_2hop && nunmatched as f64 > UNMATCHED_THRESHOLD_2HOP * num_vertices as f64 {
        matching2hop(ctrl, graph, &perm, &mut matching, &mut cnum_vertices, &mut nunmatched);
    }

    // Final renumbering: assign sequential coarse vertex IDs
    let mut new_cnum_vertices: Idx = 0;
    for i in 0..num_vertices {
        if matching[i] == UNMATCHED {
            matching[i] = i as Idx;
            graph.coarse_map[i] = new_cnum_vertices;
            new_cnum_vertices += 1;
        } else if i as Idx <= matching[i] {
            graph.coarse_map[i] = new_cnum_vertices;
            graph.coarse_map[matching[i] as usize] = new_cnum_vertices;
            new_cnum_vertices += 1;
        }
    }

    graph.matching = matching;
}

/// Random Matching (RM).
fn matchingrm(ctrl: &mut Control, graph: &mut GraphData) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    let mut matching = vec![UNMATCHED; num_vertices];
    let mut perm = vec![0 as Idx; num_vertices];
    let mut tperm = vec![0 as Idx; num_vertices];

    let nshuffles = num_vertices / 8;
    ctrl.rng.rand_array_permute_with_nshuffles(num_vertices, &mut tperm, 0, nshuffles, true);

    // C METIS: avgdegree = 4.0*(xadj[num_vertices]/num_vertices) — integer division first
    let avgdegree = if num_vertices > 0 && graph.num_edges > 0 {
        (4.0 * (graph.num_edges / graph.num_vertices) as f64) as Idx
    } else {
        1
    };
    let avgdegree = avgdegree.max(1);

    let mut degrees = vec![0 as Idx; num_vertices];
    for i in 0..num_vertices {
        let deg = graph.xadj[i + 1] - graph.xadj[i];
        degrees[i] = ((1.0 + deg as f64).sqrt() as Idx).min(avgdegree);
    }

    bucket_sort_perm(num_vertices, avgdegree, &degrees, &tperm, &mut perm);

    let mut cnum_vertices: Idx = 0;
    let mut last_unmatched = 0usize;
    let mut nunmatched = 0;

    for pi in 0..num_vertices {
        let i = perm[pi] as usize;
        if matching[i] != UNMATCHED {
            continue;
        }

        let mut maxidx: i64 = i as i64;

        let mut at_max = false;
        if ncon == 1 {
            if graph.vertex_weights[i] >= ctrl.max_vertex_weight[0] {
                at_max = true;
            }
        } else {
            at_max = (0..ncon).all(|j| graph.vertex_weights[i * ncon + j] >= ctrl.max_vertex_weight[j]);
        }

        if at_max {
            // Self-match
        } else if graph.xadj[i] == graph.xadj[i + 1] {
            maxidx = find_island_partner(&perm, &matching, &mut last_unmatched, num_vertices, i) as i64;
        } else {
            // RM: find first valid unmatched neighbor
            maxidx = -1;

            if ncon == 1 {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjacency[k] as usize;
                    if matching[j] == UNMATCHED
                        && graph.vertex_weights[i] + graph.vertex_weights[j] <= ctrl.max_vertex_weight[0]
                    {
                        maxidx = j as i64;
                        break;
                    }
                }
            } else {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjacency[k] as usize;
                    if matching[j] == UNMATCHED
                        && ivec_le_sum(ncon, &graph.vertex_weights, i, j, &ctrl.max_vertex_weight)
                    {
                        maxidx = j as i64;
                        break;
                    }
                }
            }

            if maxidx == -1 {
                if ncon == 1 {
                    if 2 * graph.vertex_weights[i] < ctrl.max_vertex_weight[0] {
                        nunmatched += 1;
                        maxidx = UNMATCHED as i64;
                    } else {
                        maxidx = i as i64;
                    }
                } else {
                    maxidx = i as i64;
                }
            }
        }

        if maxidx != UNMATCHED as i64 {
            let maxidx = if maxidx < 0 { i } else { maxidx as usize };
            graph.coarse_map[i] = cnum_vertices;
            graph.coarse_map[maxidx] = cnum_vertices;
            matching[i] = maxidx as Idx;
            matching[maxidx] = i as Idx;
            cnum_vertices += 1;
        }
    }

    if !ctrl.disable_2hop && nunmatched as f64 > UNMATCHED_THRESHOLD_2HOP * num_vertices as f64 {
        matching2hop(ctrl, graph, &perm, &mut matching, &mut cnum_vertices, &mut nunmatched);
    }

    let mut new_cnum_vertices: Idx = 0;
    for i in 0..num_vertices {
        if matching[i] == UNMATCHED {
            matching[i] = i as Idx;
            graph.coarse_map[i] = new_cnum_vertices;
            new_cnum_vertices += 1;
        } else if i as Idx <= matching[i] {
            graph.coarse_map[i] = new_cnum_vertices;
            graph.coarse_map[matching[i] as usize] = new_cnum_vertices;
            new_cnum_vertices += 1;
        }
    }

    graph.matching = matching;
}

/// Find an island (degree-0) partner by scanning forward in perm.
fn find_island_partner(
    perm: &[Idx],
    matching: &[Idx],
    last_unmatched: &mut usize,
    num_vertices: usize,
    i: usize,
) -> usize {
    while *last_unmatched < num_vertices {
        let j = perm[*last_unmatched] as usize;
        if matching[j] == UNMATCHED && j != i {
            return j;
        }
        *last_unmatched += 1;
    }
    i // Self-match if no partner found
}

/// 2-hop matching for pathological cases.
fn matching2hop(
    _ctrl: &mut Control,
    graph: &mut GraphData,
    perm: &[Idx],
    matching: &mut [Idx],
    cnum_vertices: &mut Idx,
    nunmatched: &mut usize,
) {
    let num_vertices = graph.num_vertices as usize;

    // Match_2HopAny: match vertices sharing a common neighbor
    matching2hop_any(graph, perm, matching, cnum_vertices, nunmatched, 2);
    matching2hop_all(graph, perm, matching, cnum_vertices, nunmatched, 64);
    if *nunmatched as f64 > 1.5 * UNMATCHED_THRESHOLD_2HOP * num_vertices as f64 {
        matching2hop_any(graph, perm, matching, cnum_vertices, nunmatched, 3);
    }
    if *nunmatched as f64 > 2.0 * UNMATCHED_THRESHOLD_2HOP * num_vertices as f64 {
        matching2hop_any(graph, perm, matching, cnum_vertices, nunmatched, num_vertices);
    }
}

fn matching2hop_any(
    graph: &mut GraphData,
    perm: &[Idx],
    matching: &mut [Idx],
    cnum_vertices: &mut Idx,
    nunmatched: &mut usize,
    maxdegree: usize,
) {
    let num_vertices = graph.num_vertices as usize;

    // Build inverted index: for each node v, which unmatched vertices are neighbors of v?
    let mut colptr = vec![0 as Idx; num_vertices + 1];

    // Count phase
    for pi in 0..num_vertices {
        let i = perm[pi] as usize;
        if matching[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg >= maxdegree {
            continue;
        }
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            colptr[graph.adjacency[k] as usize + 1] += 1;
        }
    }

    // Prefix sum
    for i in 0..num_vertices {
        colptr[i + 1] += colptr[i];
    }

    let total = colptr[num_vertices] as usize;
    if total == 0 {
        return;
    }
    let mut rowind = vec![0 as Idx; total];
    let mut colptr_copy = colptr.clone();

    // Fill phase
    for pi in 0..num_vertices {
        let i = perm[pi] as usize;
        if matching[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg >= maxdegree {
            continue;
        }
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let v = graph.adjacency[k] as usize;
            rowind[colptr_copy[v] as usize] = i as Idx;
            colptr_copy[v] += 1;
        }
    }

    // Match pairs sharing a common neighbor
    for v in 0..num_vertices {
        let start = colptr[v] as usize;
        let end = colptr[v + 1] as usize;
        if end - start < 2 {
            continue;
        }

        let mut j = start;
        let mut jj = end - 1;

        while j < jj {
            // Find unmatched from front
            while j < jj && matching[rowind[j] as usize] != UNMATCHED {
                j += 1;
            }
            // Find unmatched from back
            while jj > j && matching[rowind[jj] as usize] != UNMATCHED {
                jj -= 1;
            }

            if j < jj {
                let vi = rowind[j] as usize;
                let vj = rowind[jj] as usize;
                graph.coarse_map[vi] = *cnum_vertices;
                graph.coarse_map[vj] = *cnum_vertices;
                matching[vi] = vj as Idx;
                matching[vj] = vi as Idx;
                *cnum_vertices += 1;
                *nunmatched -= 2;
                j += 1;
                jj -= 1;
            }
        }
    }
}

fn matching2hop_all(
    graph: &mut GraphData,
    perm: &[Idx],
    matching: &mut [Idx],
    cnum_vertices: &mut Idx,
    nunmatched: &mut usize,
    maxdegree: usize,
) {
    let num_vertices = graph.num_vertices as usize;

    // Collect candidates with degree > 1 and < maxdegree
    struct KeyVal { key: Idx, val: Idx }
    let mut cands: Vec<KeyVal> = Vec::new();
    let mask = if maxdegree > 0 { Idx::MAX / maxdegree as Idx } else { 1 };

    for pi in 0..num_vertices {
        let i = perm[pi] as usize;
        if matching[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg <= 1 || deg >= maxdegree {
            continue;
        }
        let mut k: Idx = 0;
        for jj in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            k += graph.adjacency[jj] % mask;
        }
        cands.push(KeyVal {
            key: (k % mask) * maxdegree as Idx + deg as Idx,
            val: i as Idx,
        });
    }

    cands.sort_by_key(|c| c.key);

    let mut mark = vec![-1 as Idx; num_vertices];

    let ncand = cands.len();
    let mut pi = 0;
    while pi < ncand {
        let i = cands[pi].val as usize;
        if matching[i] != UNMATCHED {
            pi += 1;
            continue;
        }

        // Mark neighbors
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            mark[graph.adjacency[k] as usize] = i as Idx;
        }

        let ideg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;

        let mut pk = pi + 1;
        while pk < ncand && cands[pk].key == cands[pi].key {
            let kk = cands[pk].val as usize;
            if matching[kk] != UNMATCHED {
                pk += 1;
                continue;
            }

            let kdeg = (graph.xadj[kk + 1] - graph.xadj[kk]) as usize;
            if kdeg != ideg {
                pk += 1;
                continue;
            }

            // Verify identical adjacency
            let mut identical = true;
            for jj in graph.xadj[kk] as usize..graph.xadj[kk + 1] as usize {
                if mark[graph.adjacency[jj] as usize] != i as Idx {
                    identical = false;
                    break;
                }
            }

            if identical {
                graph.coarse_map[i] = *cnum_vertices;
                graph.coarse_map[kk] = *cnum_vertices;
                matching[i] = kk as Idx;
                matching[kk] = i as Idx;
                *cnum_vertices += 1;
                *nunmatched -= 2;
                break;
            }
            pk += 1;
        }

        // Clean mark
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            mark[graph.adjacency[k] as usize] = -1;
        }
        pi += 1;
    }
}

/// Check that vertex_weights[i*ncon..] + vertex_weights[j*ncon..] <= max_vertex_weight[..] for all constraints.
fn ivec_le_sum(ncon: usize, vertex_weights: &[Idx], i: usize, j: usize, max_vertex_weight: &[Idx]) -> bool {
    for k in 0..ncon {
        if vertex_weights[i * ncon + k] + vertex_weights[j * ncon + k] > max_vertex_weight[k] {
            return false;
        }
    }
    true
}

/// Bucket sort: sort input by increasing keys, producing perm.
fn bucket_sort_perm(n: usize, max_key: Idx, keys: &[Idx], input: &[Idx], output: &mut [Idx]) {
    let nbuckets = max_key as usize + 1;
    let mut counts = vec![0usize; nbuckets];

    for i in 0..n {
        counts[keys[input[i] as usize] as usize] += 1;
    }

    let mut sum = 0;
    for i in 0..nbuckets {
        let c = counts[i];
        counts[i] = sum;
        sum += c;
    }

    for i in 0..n {
        let k = keys[input[i] as usize] as usize;
        output[counts[k]] = input[i];
        counts[k] += 1;
    }
}
