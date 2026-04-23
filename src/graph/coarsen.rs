use crate::types::Idx;
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::graph::contract;

const COARSEN_FRACTION: f64 = 0.85;
const UNMATCHED: Idx = -1;
const UNMATCHEDFOR2HOP: f64 = 0.10;

/// Coarsen the graph iteratively. Returns pointer to coarsest graph.
pub fn coarsen_graph(ctrl: &mut Ctrl, graph: &mut GraphData) {
    let ncon = graph.ncon as usize;

    // Check if all edge weights are equal
    let eqewgts = if graph.nedges > 0 {
        let first = graph.adjwgt[0];
        graph.adjwgt.iter().all(|&w| w == first)
    } else {
        true
    };

    // Compute maxvwgt: 1.5 * tvwgt[i] / CoarsenTo
    for i in 0..ncon {
        ctrl.maxvwgt[i] = ((1.5 * graph.tvwgt[i] as f64) / ctrl.coarsen_to.max(1) as f64) as Idx;
    }

    let mut cur = graph as *mut GraphData;
    let mut cur_eqewgts = eqewgts;

    loop {
        let g = unsafe { &mut *cur };

        let nvtxs = g.nvtxs as usize;
        g.cmap = vec![0; nvtxs];

        // Choose matching
        if ctrl.ctype == 1 && !cur_eqewgts && g.nedges > 0 {
            match_shem(ctrl, g);
        } else {
            match_rm(ctrl, g);
        }

        // Build coarse graph
        let cnvtxs = g.cmap.iter().copied().max().unwrap_or(-1) + 1;
        let coarse = contract::create_coarse_graph(ctrl, g, cnvtxs);
        g.coarser = Some(coarse);
        g.coarser.as_mut().unwrap().finer = cur;

        cur_eqewgts = false;

        let coarser = g.coarser.as_ref().unwrap();
        let continue_coarsening =
            coarser.nvtxs > ctrl.coarsen_to
            && (coarser.nvtxs as f64) < COARSEN_FRACTION * g.nvtxs as f64
            && coarser.nedges > coarser.nvtxs / 2;

        if !continue_coarsening {
            break;
        }

        cur = &mut **g.coarser.as_mut().unwrap() as *mut GraphData;
    }
}

/// Sorted Heavy Edge Matching (SHEM).
fn match_shem(ctrl: &mut Ctrl, graph: &mut GraphData) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    let mut match_ = vec![UNMATCHED; nvtxs];
    let mut perm = vec![0 as Idx; nvtxs];
    let mut tperm = vec![0 as Idx; nvtxs];

    // Create degree-biased random permutation
    // Step 1: random permutation of tperm
    let nshuffles = nvtxs / 8;
    ctrl.rng.rand_array_permute_with_nshuffles(nvtxs, &mut tperm, 0, nshuffles, true);

    // Step 2: degree-based bucket sort
    // C METIS: avgdegree = 4.0*(xadj[nvtxs]/nvtxs) — integer division first
    let avgdegree = if nvtxs > 0 && graph.nedges > 0 {
        (4.0 * (graph.nedges / graph.nvtxs) as f64) as Idx
    } else {
        1
    };
    let avgdegree = avgdegree.max(1);

    let mut degrees = vec![0 as Idx; nvtxs];
    for i in 0..nvtxs {
        let deg = graph.xadj[i + 1] - graph.xadj[i];
        let bnum = ((1.0 + deg as f64).sqrt()) as Idx;
        degrees[i] = bnum.min(avgdegree);
    }

    // Bucket sort tperm by degrees (stable, increasing)
    bucket_sort_perm(nvtxs, avgdegree, &degrees, &tperm, &mut perm);

    let mut cnvtxs: Idx = 0;
    let mut last_unmatched = 0usize;
    let mut nunmatched = 0;

    for pi in 0..nvtxs {
        let i = perm[pi] as usize;
        if match_[i] != UNMATCHED {
            continue;
        }

        let mut maxidx: i64 = i as i64;

        // Weight guard: skip if already at maxvwgt
        let mut at_max = false;
        if ncon == 1 {
            if graph.vwgt[i] >= ctrl.maxvwgt[0] {
                at_max = true;
            }
        } else {
            at_max = (0..ncon).all(|j| graph.vwgt[i * ncon + j] >= ctrl.maxvwgt[j]);
        }

        if at_max {
            // Self-match
        } else if graph.xadj[i] == graph.xadj[i + 1] {
            // Island vertex: find another unmatched vertex
            maxidx = find_island_partner(&perm, &match_, &mut last_unmatched, nvtxs, i) as i64;
        } else {
            // SHEM: find heaviest unmatched neighbor
            let mut maxwgt: Idx = -1;
            maxidx = -1; // Will be set to partner or UNMATCHED for 2-hop candidate

            if ncon == 1 {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjncy[k] as usize;
                    if match_[j] == UNMATCHED
                        && maxwgt < graph.adjwgt[k]
                        && graph.vwgt[i] + graph.vwgt[j] <= ctrl.maxvwgt[0]
                    {
                        maxidx = j as i64;
                        maxwgt = graph.adjwgt[k];
                    }
                }
            } else {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjncy[k] as usize;
                    if match_[j] == UNMATCHED
                        && ivec_le_sum(ncon, &graph.vwgt, i, j, &ctrl.maxvwgt)
                        && maxwgt < graph.adjwgt[k]
                    {
                        maxidx = j as i64;
                        maxwgt = graph.adjwgt[k];
                    }
                }
            }

            if maxidx == -1 {
                // Check if 2-hop candidate
                if ncon == 1 {
                    if 2 * graph.vwgt[i] < ctrl.maxvwgt[0] {
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
            graph.cmap[i] = cnvtxs;
            graph.cmap[maxidx] = cnvtxs;
            match_[i] = maxidx as Idx;
            match_[maxidx] = i as Idx;
            cnvtxs += 1;
        }
    }

    // 2-hop matching
    if !ctrl.no2hop && nunmatched as f64 > UNMATCHEDFOR2HOP * nvtxs as f64 {
        match_2hop(ctrl, graph, &perm, &mut match_, &mut cnvtxs, &mut nunmatched);
    }

    // Final renumbering: assign sequential coarse vertex IDs
    let mut new_cnvtxs: Idx = 0;
    for i in 0..nvtxs {
        if match_[i] == UNMATCHED {
            match_[i] = i as Idx;
            graph.cmap[i] = new_cnvtxs;
            new_cnvtxs += 1;
        } else if i as Idx <= match_[i] {
            graph.cmap[i] = new_cnvtxs;
            graph.cmap[match_[i] as usize] = new_cnvtxs;
            new_cnvtxs += 1;
        }
    }

    graph.match_ = match_;
}

/// Random Matching (RM).
fn match_rm(ctrl: &mut Ctrl, graph: &mut GraphData) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    let mut match_ = vec![UNMATCHED; nvtxs];
    let mut perm = vec![0 as Idx; nvtxs];
    let mut tperm = vec![0 as Idx; nvtxs];

    let nshuffles = nvtxs / 8;
    ctrl.rng.rand_array_permute_with_nshuffles(nvtxs, &mut tperm, 0, nshuffles, true);

    // C METIS: avgdegree = 4.0*(xadj[nvtxs]/nvtxs) — integer division first
    let avgdegree = if nvtxs > 0 && graph.nedges > 0 {
        (4.0 * (graph.nedges / graph.nvtxs) as f64) as Idx
    } else {
        1
    };
    let avgdegree = avgdegree.max(1);

    let mut degrees = vec![0 as Idx; nvtxs];
    for i in 0..nvtxs {
        let deg = graph.xadj[i + 1] - graph.xadj[i];
        degrees[i] = ((1.0 + deg as f64).sqrt() as Idx).min(avgdegree);
    }

    bucket_sort_perm(nvtxs, avgdegree, &degrees, &tperm, &mut perm);

    let mut cnvtxs: Idx = 0;
    let mut last_unmatched = 0usize;
    let mut nunmatched = 0;

    for pi in 0..nvtxs {
        let i = perm[pi] as usize;
        if match_[i] != UNMATCHED {
            continue;
        }

        let mut maxidx: i64 = i as i64;

        let mut at_max = false;
        if ncon == 1 {
            if graph.vwgt[i] >= ctrl.maxvwgt[0] {
                at_max = true;
            }
        } else {
            at_max = (0..ncon).all(|j| graph.vwgt[i * ncon + j] >= ctrl.maxvwgt[j]);
        }

        if at_max {
            // Self-match
        } else if graph.xadj[i] == graph.xadj[i + 1] {
            maxidx = find_island_partner(&perm, &match_, &mut last_unmatched, nvtxs, i) as i64;
        } else {
            // RM: find first valid unmatched neighbor
            maxidx = -1;

            if ncon == 1 {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjncy[k] as usize;
                    if match_[j] == UNMATCHED
                        && graph.vwgt[i] + graph.vwgt[j] <= ctrl.maxvwgt[0]
                    {
                        maxidx = j as i64;
                        break;
                    }
                }
            } else {
                for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                    let j = graph.adjncy[k] as usize;
                    if match_[j] == UNMATCHED
                        && ivec_le_sum(ncon, &graph.vwgt, i, j, &ctrl.maxvwgt)
                    {
                        maxidx = j as i64;
                        break;
                    }
                }
            }

            if maxidx == -1 {
                if ncon == 1 {
                    if 2 * graph.vwgt[i] < ctrl.maxvwgt[0] {
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
            graph.cmap[i] = cnvtxs;
            graph.cmap[maxidx] = cnvtxs;
            match_[i] = maxidx as Idx;
            match_[maxidx] = i as Idx;
            cnvtxs += 1;
        }
    }

    if !ctrl.no2hop && nunmatched as f64 > UNMATCHEDFOR2HOP * nvtxs as f64 {
        match_2hop(ctrl, graph, &perm, &mut match_, &mut cnvtxs, &mut nunmatched);
    }

    let mut new_cnvtxs: Idx = 0;
    for i in 0..nvtxs {
        if match_[i] == UNMATCHED {
            match_[i] = i as Idx;
            graph.cmap[i] = new_cnvtxs;
            new_cnvtxs += 1;
        } else if i as Idx <= match_[i] {
            graph.cmap[i] = new_cnvtxs;
            graph.cmap[match_[i] as usize] = new_cnvtxs;
            new_cnvtxs += 1;
        }
    }

    graph.match_ = match_;
}

/// Find an island (degree-0) partner by scanning forward in perm.
fn find_island_partner(
    perm: &[Idx],
    match_: &[Idx],
    last_unmatched: &mut usize,
    nvtxs: usize,
    i: usize,
) -> usize {
    while *last_unmatched < nvtxs {
        let j = perm[*last_unmatched] as usize;
        if match_[j] == UNMATCHED && j != i {
            return j;
        }
        *last_unmatched += 1;
    }
    i // Self-match if no partner found
}

/// 2-hop matching for pathological cases.
fn match_2hop(
    _ctrl: &mut Ctrl,
    graph: &mut GraphData,
    perm: &[Idx],
    match_: &mut [Idx],
    cnvtxs: &mut Idx,
    nunmatched: &mut usize,
) {
    let nvtxs = graph.nvtxs as usize;

    // Match_2HopAny: match vertices sharing a common neighbor
    match_2hop_any(graph, perm, match_, cnvtxs, nunmatched, 2);
    match_2hop_all(graph, perm, match_, cnvtxs, nunmatched, 64);
    if *nunmatched as f64 > 1.5 * UNMATCHEDFOR2HOP * nvtxs as f64 {
        match_2hop_any(graph, perm, match_, cnvtxs, nunmatched, 3);
    }
    if *nunmatched as f64 > 2.0 * UNMATCHEDFOR2HOP * nvtxs as f64 {
        match_2hop_any(graph, perm, match_, cnvtxs, nunmatched, nvtxs);
    }
}

fn match_2hop_any(
    graph: &mut GraphData,
    perm: &[Idx],
    match_: &mut [Idx],
    cnvtxs: &mut Idx,
    nunmatched: &mut usize,
    maxdegree: usize,
) {
    let nvtxs = graph.nvtxs as usize;

    // Build inverted index: for each node v, which unmatched vertices are neighbors of v?
    let mut colptr = vec![0 as Idx; nvtxs + 1];

    // Count phase
    for pi in 0..nvtxs {
        let i = perm[pi] as usize;
        if match_[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg >= maxdegree {
            continue;
        }
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            colptr[graph.adjncy[k] as usize + 1] += 1;
        }
    }

    // Prefix sum
    for i in 0..nvtxs {
        colptr[i + 1] += colptr[i];
    }

    let total = colptr[nvtxs] as usize;
    if total == 0 {
        return;
    }
    let mut rowind = vec![0 as Idx; total];
    let mut colptr_copy = colptr.clone();

    // Fill phase
    for pi in 0..nvtxs {
        let i = perm[pi] as usize;
        if match_[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg >= maxdegree {
            continue;
        }
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let v = graph.adjncy[k] as usize;
            rowind[colptr_copy[v] as usize] = i as Idx;
            colptr_copy[v] += 1;
        }
    }

    // Match pairs sharing a common neighbor
    for v in 0..nvtxs {
        let start = colptr[v] as usize;
        let end = colptr[v + 1] as usize;
        if end - start < 2 {
            continue;
        }

        let mut j = start;
        let mut jj = end - 1;

        while j < jj {
            // Find unmatched from front
            while j < jj && match_[rowind[j] as usize] != UNMATCHED {
                j += 1;
            }
            // Find unmatched from back
            while jj > j && match_[rowind[jj] as usize] != UNMATCHED {
                jj -= 1;
            }

            if j < jj {
                let vi = rowind[j] as usize;
                let vj = rowind[jj] as usize;
                graph.cmap[vi] = *cnvtxs;
                graph.cmap[vj] = *cnvtxs;
                match_[vi] = vj as Idx;
                match_[vj] = vi as Idx;
                *cnvtxs += 1;
                *nunmatched -= 2;
                j += 1;
                jj -= 1;
            }
        }
    }
}

fn match_2hop_all(
    graph: &mut GraphData,
    perm: &[Idx],
    match_: &mut [Idx],
    cnvtxs: &mut Idx,
    nunmatched: &mut usize,
    maxdegree: usize,
) {
    let nvtxs = graph.nvtxs as usize;

    // Collect candidates with degree > 1 and < maxdegree
    struct KeyVal { key: Idx, val: Idx }
    let mut cands: Vec<KeyVal> = Vec::new();
    let mask = if maxdegree > 0 { Idx::MAX / maxdegree as Idx } else { 1 };

    for pi in 0..nvtxs {
        let i = perm[pi] as usize;
        if match_[i] != UNMATCHED {
            continue;
        }
        let deg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;
        if deg <= 1 || deg >= maxdegree {
            continue;
        }
        let mut k: Idx = 0;
        for jj in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            k += graph.adjncy[jj] % mask;
        }
        cands.push(KeyVal {
            key: (k % mask) * maxdegree as Idx + deg as Idx,
            val: i as Idx,
        });
    }

    cands.sort_by_key(|c| c.key);

    let mut mark = vec![-1 as Idx; nvtxs];

    let ncand = cands.len();
    let mut pi = 0;
    while pi < ncand {
        let i = cands[pi].val as usize;
        if match_[i] != UNMATCHED {
            pi += 1;
            continue;
        }

        // Mark neighbors
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            mark[graph.adjncy[k] as usize] = i as Idx;
        }

        let ideg = (graph.xadj[i + 1] - graph.xadj[i]) as usize;

        let mut pk = pi + 1;
        while pk < ncand && cands[pk].key == cands[pi].key {
            let kk = cands[pk].val as usize;
            if match_[kk] != UNMATCHED {
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
                if mark[graph.adjncy[jj] as usize] != i as Idx {
                    identical = false;
                    break;
                }
            }

            if identical {
                graph.cmap[i] = *cnvtxs;
                graph.cmap[kk] = *cnvtxs;
                match_[i] = kk as Idx;
                match_[kk] = i as Idx;
                *cnvtxs += 1;
                *nunmatched -= 2;
                break;
            }
            pk += 1;
        }

        // Clean mark
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            mark[graph.adjncy[k] as usize] = -1;
        }
        pi += 1;
    }
}

/// Check that vwgt[i*ncon..] + vwgt[j*ncon..] <= maxvwgt[..] for all constraints.
fn ivec_le_sum(ncon: usize, vwgt: &[Idx], i: usize, j: usize, maxvwgt: &[Idx]) -> bool {
    for k in 0..ncon {
        if vwgt[i * ncon + k] + vwgt[j * ncon + k] > maxvwgt[k] {
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
