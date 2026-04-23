use crate::types::{Idx, Real, NOPTIONS};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::graph::coarsen::coarsen_graph;
use crate::partition::initpart;
use crate::partition::refine2way;
use crate::partition::split;

/// SMALLNIPARTS / LARGENIPARTS from C METIS
const SMALLNIPARTS: Idx = 5;
const LARGENIPARTS: Idx = 7;

/// Initialize a k-way partition by recursive bisection on the coarsest graph.
///
/// This replicates InitKWayPartitioning from kmetis.c:
/// - Creates a new Ctrl with PMETIS settings (optype=1, ufactor=1, rtype=FM)
/// - Computes ubvec for recursive bisection: ubvec[i] = pow(ctrl.ubfactors[i], 1.0/log(nparts))
/// - Sets options[NCUTS] = ctrl.nIparts and options[NITER] = ctrl.niter
/// - Computes coarsen_to = max(nvtxs/(40*log2(nparts)), 20*nparts)
/// - Runs METIS_PartGraphRecursive (MlevelRecursiveBisection) on the graph
/// - Result goes into graph.where_
pub fn init_kway_partition(ctrl: &mut Ctrl, graph: &mut GraphData, nparts: Idx) {
    let nvtxs = graph.nvtxs as usize;
    graph.where_ = vec![0; nvtxs];

    if nparts == 1 {
        return;
    }

    let ncon = graph.ncon as usize;

    // Create a temporary PMETIS ctrl, replicating SetupCtrl for METIS_OP_PMETIS.
    // In C METIS, InitKWayPartitioning calls METIS_PartGraphRecursive on the
    // coarsest graph, which creates a fresh ctrl with PMETIS settings.
    let mut options = [-1 as Idx; NOPTIONS];

    // options[NCUTS] = ctrl->nIparts (the niparts from the kway ctrl)
    // In METIS_PartGraphKway, nIparts defaults to -1, then SetupCtrl sets it
    // based on the graph. For InitKWayPartitioning, ncuts = ctrl->nIparts.
    // If niparts is -1 (default), METIS computes it as 4 or 5 depending on
    // whether the graph is large. We use the parent ctrl's niparts if set,
    // otherwise default to 1 (the PMETIS ncuts default).
    let ncuts = if ctrl.niparts > 0 { ctrl.niparts } else { 1 };
    options[8] = ncuts; // NCuts index = 8

    // options[NITER] = ctrl->niter
    options[7] = ctrl.niter; // NIter index = 7

    // Seed: C METIS does NOT propagate the seed to the inner PMETIS ctrl.
    // options[SEED] remains -1, so InitRandom(-1) uses default seed 4321.
    // Do NOT set options[9].

    // Create PMETIS ctrl (is_kway = false => optype=1, ufactor=1, rtype=FM)
    let mut pctrl = Ctrl::new(&options, graph.ncon, nparts, false);

    // Compute ubvec for recursive bisection:
    // ubvec[i] = pow(ctrl->ubfactors[i], 1.0/log(nparts))
    let log_nparts = (nparts as f64).ln();
    if log_nparts > 0.0 {
        let inv_log = 1.0 / log_nparts;
        for i in 0..ncon {
            pctrl.ubfactors[i] = (ctrl.ubfactors[i] as f64).powf(inv_log) as Real;
            // C METIS SetupCtrl always adds 0.0000499 AFTER any ubvec override
            pctrl.ubfactors[i] += 0.0000499;
        }
    }

    // C METIS SetupCtrl for METIS_OP_PMETIS hardcodes CoarsenTo:
    //   ncon == 1: CoarsenTo = 20
    //   ncon >  1: CoarsenTo = 100
    // METIS_PartGraphRecursive does NOT override this value.
    pctrl.coarsen_to = if ncon == 1 { 20 } else { 100 };

    // Compute tpwgts: uniform 1/nparts for each partition and constraint
    pctrl.tpwgts = vec![1.0 / nparts as Real; (nparts * graph.ncon) as usize];

    // In C METIS, METIS_PartGraphRecursive calls SetupCtrl which resets the
    // global RNG to the original seed. pctrl already has the correct seed from
    // Ctrl::new (via options[9] = ctrl.seed). No override needed.

    // Always initialize labels to identity for recursive bisection.
    // contract.rs sets labels to vec![0; n], so we must overwrite unconditionally.
    graph.label = (0..graph.nvtxs).collect();

    // Run recursive bisection: MlevelRecursiveBisection
    let mut part = vec![0 as Idx; nvtxs];
    let tpwgts = pctrl.tpwgts.clone();
    mlevel_recursive_bisection(&mut pctrl, graph, nparts, &mut part, &tpwgts, 0);

    // In C METIS, the global RNG state after recursive bisection persists.
    // Propagate pctrl's RNG state back to the parent ctrl so RefineKWay
    // sees the same RNG state as C METIS.
    ctrl.rng = pctrl.rng;

    // Copy result into graph.where_
    graph.where_.copy_from_slice(&part);
}

/// Recursively bisect graph into nparts partitions.
///
/// Replicates MlevelRecursiveBisection from pmetis.c:
/// 1. Compute tpwgts2[i] = sum of first half of tpwgts for each constraint
/// 2. Call MultilevelBisect(ctrl, graph, tpwgts2)
/// 3. Map partition back using label[i]: part[label[i]] = where[i] + fpart
/// 4. If nparts > 2: SplitGraphPart, then recurse on each subgraph
/// 5. Scale tpwgts for recursive calls
fn mlevel_recursive_bisection(
    ctrl: &mut Ctrl,
    graph: &mut GraphData,
    nparts: Idx,
    part: &mut [Idx],
    tpwgts: &[Real],
    fpart: Idx,
) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Base case: single partition
    if nparts == 1 {
        // Map using label: part[label[i]] = fpart
        for i in 0..nvtxs {
            let label_i = graph.label[i] as usize;
            if label_i < part.len() {
                part[label_i] = fpart;
            }
        }
        return;
    }

    // Step 1: Compute tpwgts2 for the bisection.
    // tpwgts2[0*ncon+j] = sum of tpwgts[i*ncon+j] for i in 0..nparts/2
    // tpwgts2[1*ncon+j] = sum of tpwgts[i*ncon+j] for i in nparts/2..nparts
    let nparts_left = nparts / 2;
    let nparts_right = nparts - nparts_left;

    let mut tpwgts2 = vec![0.0 as Real; 2 * ncon];
    for j in 0..ncon {
        for i in 0..nparts_left as usize {
            tpwgts2[j] += tpwgts[i * ncon + j];
        }
        for i in nparts_left as usize..nparts as usize {
            tpwgts2[ncon + j] += tpwgts[i * ncon + j];
        }
    }

    // Step 2: MultilevelBisect
    multilevel_bisect(ctrl, graph, &tpwgts2);
    // Step 3: Map partition back using label.
    // part[label[i]] = where[i] + fpart
    for i in 0..nvtxs {
        let label_i = graph.label[i] as usize;
        if label_i < part.len() {
            part[label_i] = graph.where_[i] + fpart;
        }
    }

    // Step 4: If nparts > 2, split and recurse
    if nparts > 2 {
        let (mut lgraph, mut rgraph) = split::split_graph_part(graph);

        // Step 5: Scale tpwgts for recursive calls.
        // Left subgraph gets tpwgts[0..nparts_left], normalized by sum(tpwgts2[0..ncon])
        // Right subgraph gets tpwgts[nparts_left..nparts], normalized by sum(tpwgts2[ncon..2*ncon])
        let mut left_tpwgts = vec![0.0 as Real; (nparts_left * graph.ncon) as usize];
        let mut right_tpwgts = vec![0.0 as Real; (nparts_right * graph.ncon) as usize];

        for j in 0..ncon {
            let wgt_left = if tpwgts2[j] > 0.0 { 1.0 / tpwgts2[j] } else { 0.0 };
            let wgt_right = if tpwgts2[ncon + j] > 0.0 { 1.0 / tpwgts2[ncon + j] } else { 0.0 };

            for i in 0..nparts_left as usize {
                left_tpwgts[i * ncon + j] = tpwgts[i * ncon + j] * wgt_left;
            }
            for i in 0..nparts_right as usize {
                right_tpwgts[i * ncon + j] =
                    tpwgts[(nparts_left as usize + i) * ncon + j] * wgt_right;
            }
        }

        // Recurse on left subgraph
        mlevel_recursive_bisection(
            ctrl,
            &mut lgraph,
            nparts_left,
            part,
            &left_tpwgts,
            fpart,
        );

        // Recurse on right subgraph
        mlevel_recursive_bisection(
            ctrl,
            &mut rgraph,
            nparts_right,
            part,
            &right_tpwgts,
            fpart + nparts_left,
        );
    }
}

/// Multilevel bisection: coarsen -> init partition -> refine.
///
/// Replicates MultilevelBisect from pmetis.c:
/// 1. Setup2WayBalMultipliers(ctrl, graph, tpwgts)
/// 2. For ncuts iterations:
///    a. cgraph = CoarsenGraph(ctrl, graph)
///    b. niparts = (cgraph->nvtxs <= ctrl->CoarsenTo) ? SMALLNIPARTS : LARGENIPARTS
///    c. Init2WayPartition(ctrl, cgraph, tpwgts, niparts)
///    d. Refine2Way(ctrl, graph, cgraph, tpwgts) - uncoarsens and refines
///    e. Keep best based on balance/cut
/// Matches C METIS MultilevelBisect from pmetis.c:
/// - Uses ComputeLoadImbalanceDiff for balance comparison
/// - Selection: i==0 || (curbal<=0.0005 && bestobj>curobj) || (bestbal>0.0005 && curbal<bestbal)
/// - Early break on bestobj==0
fn multilevel_bisect(ctrl: &mut Ctrl, graph: &mut GraphData, tpwgts: &[Real]) {
    let ncuts = ctrl.ncuts;

    // Step 1: Setup2WayBalMultipliers
    ctrl.setup_2way_bal_multipliers(&graph.invtvwgt, tpwgts);

    let mut best_where: Vec<Idx> = Vec::new();
    let mut bestobj: Idx = 0;
    let mut bestbal: Real = 0.0;

    // Track last iteration's result for C METIS-compatible restoration logic.
    // C METIS works on the same graph object in-place. After the loop, graph->where
    // contains the LAST iteration's partition. It only restores bestwhere if
    // bestobj != curobj (last iteration's obj). We must replicate this.
    let mut last_curobj: Idx = 0;
    let mut last_where: Vec<Idx> = Vec::new();

    // Step 2: For ncuts iterations
    for i in 0..ncuts {
        // Step 2a: Make a working copy and coarsen
        let mut work = clone_graph_for_bisect(graph);
        coarsen_graph(ctrl, &mut work);

        let orggraph_ptr: *mut GraphData = &mut work as *mut GraphData;
        let cgraph_ptr: *mut GraphData = get_coarsest_ptr(orggraph_ptr);
        let cgraph = unsafe { &mut *cgraph_ptr };

        // Step 2b: Compute niparts based on coarsest graph size
        let niparts = if cgraph.nvtxs <= ctrl.coarsen_to {
            SMALLNIPARTS
        } else {
            LARGENIPARTS
        };

        // Step 2c: Init2WayPartition on the coarsest graph
        initpart::init_2way_partition(ctrl, cgraph, tpwgts, niparts);

        // Step 2d: Refine2Way - uncoarsen and refine through the chain.
        refine2way::refine_2way(ctrl, orggraph_ptr, cgraph, tpwgts);

        // After refinement, work (orggraph) has the final partition result.
        let curobj = work.mincut;
        let curbal = compute_load_imbalance_diff_2way(&work, &ctrl.pijbm, &ctrl.ubfactors);

        // Step 2e: Keep best (matching C METIS criterion exactly)
        // C METIS only saves bestwhere on non-last iterations (i < ncuts-1).
        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            bestobj = curobj;
            bestbal = curbal;
            if i < ncuts - 1 {
                best_where = work.where_.clone();
            }
        }

        // Always track last iteration's result
        last_curobj = curobj;
        last_where = work.where_;

        if bestobj == 0 {
            break;
        }
    }

    // Post-loop restoration: match C METIS exactly.
    // C works on same graph, so graph->where has last iteration's result.
    // It only restores bestwhere if bestobj != curobj (last iter's obj).
    if bestobj != last_curobj {
        graph.where_ = best_where;
    } else {
        graph.where_ = last_where;
    }
    graph.mincut = bestobj;
}

/// ComputeLoadImbalanceDiff for 2-way partition.
/// Returns max over all partitions and constraints of:
///   pwgts[i*ncon+j] * pijbm[i*ncon+j] - ubfactors[j]
/// Negative means balanced, positive means imbalanced.
fn compute_load_imbalance_diff_2way(graph: &GraphData, pijbm: &[Real], ubfactors: &[Real]) -> Real {
    let ncon = graph.ncon as usize;
    let mut max_diff: Real = -1.0 * ncon as Real;
    for i in 0..2 {
        for j in 0..ncon {
            let idx = i * ncon + j;
            if idx < graph.pwgts.len() && idx < pijbm.len() && j < ubfactors.len() {
                let diff = graph.pwgts[idx] as Real * pijbm[idx] - ubfactors[j];
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }
    }
    max_diff
}

/// Clone graph data for bisection (without coarsening chain).
fn clone_graph_for_bisect(g: &GraphData) -> GraphData {
    let mut ng = GraphData::new();
    ng.nvtxs = g.nvtxs;
    ng.nedges = g.nedges;
    ng.ncon = g.ncon;
    ng.xadj = g.xadj.clone();
    ng.adjncy = g.adjncy.clone();
    ng.vwgt = g.vwgt.clone();
    ng.vsize = g.vsize.clone();
    ng.adjwgt = g.adjwgt.clone();
    ng.tvwgt = g.tvwgt.clone();
    ng.invtvwgt = g.invtvwgt.clone();
    ng.label = (0..g.nvtxs).collect();
    ng
}

/// Traverse coarsening chain to find coarsest graph (raw pointer version).
/// This avoids borrow checker issues when we need both the finest and coarsest
/// graphs from the same chain simultaneously.
fn get_coarsest_ptr(g: *mut GraphData) -> *mut GraphData {
    let mut cur = g;
    unsafe {
        while (*cur).coarser.is_some() {
            cur = &mut **(*cur).coarser.as_mut().unwrap() as *mut GraphData;
        }
    }
    cur
}
