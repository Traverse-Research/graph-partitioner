use crate::types::{Idx, Real};
use crate::ctrl::Control;
use crate::option::Options;
use crate::graph::GraphData;
use crate::graph::coarsen::coarsen_graph;
use crate::partition::initpart;
use crate::partition::refine2way;
use crate::partition::split;

/// SMALL_NUM_INIT_PARTS / LARGE_NUM_INIT_PARTS from C METIS
const SMALL_NUM_INIT_PARTS: Idx = 5;
const LARGE_NUM_INIT_PARTS: Idx = 7;

/// Initialize a k-way partition by recursive bisection on the coarsest graph.
///
/// This replicates InitKWayPartitioning from kmetis.c:
/// - Creates a new Control with PMETIS settings (optype=1, ufactor=1, rtype=FM)
/// - Computes ubvec for recursive bisection: ubvec[i] = pow(ctrl.imbalance_tols[i], 1.0/log(nparts))
/// - Sets options[NCUTS] = ctrl.nIparts and options[NITER] = ctrl.num_iter
/// - Computes coarsen_to = max(num_vertices/(40*log2(nparts)), 20*nparts)
/// - Runs METIS_PartGraphRecursive (MlevelRecursiveBisection) on the graph
/// - Result goes into graph.partition
pub fn init_kway_partition(ctrl: &mut Control, graph: &mut GraphData, nparts: Idx) {
    let num_vertices = graph.num_vertices as usize;
    graph.partition = vec![0; num_vertices];

    if nparts == 1 {
        return;
    }

    let ncon = graph.num_constraints as usize;

    // Create a temporary PMETIS ctrl, replicating SetupCtrl for METIS_OP_PMETIS.
    // In C METIS, InitKWayPartitioning calls METIS_PartGraphRecursive on the
    // coarsest graph, which creates a fresh ctrl with PMETIS settings.
    let ncuts = if ctrl.num_init_parts > 0 { ctrl.num_init_parts } else { 1 };

    // Seed: C METIS does NOT propagate the seed to the inner PMETIS ctrl.
    // seed remains None (default -1), so InitRandom(-1) uses default seed 4321.
    let pmetis_options = Options {
        num_cuts: Some(ncuts),
        num_iter: Some(ctrl.num_iter),
        ..Options::default()
    };

    // Create PMETIS ctrl (is_kway = false => optype=1, ufactor=1, rtype=FM)
    let mut pctrl = Control::new(&pmetis_options, graph.num_constraints, nparts, false);

    // C METIS SetupCtrl for METIS_OP_PMETIS uses different defaults based on ncon:
    //   ncon == 1: iptype=GROW(0),  ufactor=1,  CoarsenTo=20
    //   ncon >  1: iptype=RANDOM(1), ufactor=10, CoarsenTo=100
    if ncon > 1 {
        pctrl.init_part_type = 1; // METIS_IPTYPE_RANDOM
        // ufactor=10 (MCPMETIS_DEFAULT_UFACTOR) - but imbalance_tols are overridden below
    }
    pctrl.coarsen_to = if ncon == 1 { 20 } else { 100 };

    // Compute ubvec for recursive bisection:
    // ubvec[i] = pow(ctrl->ubfactors[i], 1.0/log(nparts))
    // In C METIS, ctrl->ubfactors[i] already includes +0.0000499 (added in SetupCtrl).
    // Then the inner PMETIS SetupCtrl adds +0.0000499 again to the result.
    let log_nparts = (nparts as f64).ln();
    if log_nparts > 0.0 {
        let inv_log = 1.0 / log_nparts;
        for i in 0..ncon {
            // C: ubvec[i] = pow(ctrl->ubfactors[i], 1.0/log(nparts))
            // ctrl->ubfactors[i] includes +0.0000499 from outer SetupCtrl.
            // Our ctrl.imbalance_tols[i] also includes +0.0000499.
            let base_ubfactor = ctrl.imbalance_tols[i] as f64;
            // C: inner SetupCtrl adds +0.0000499 to ubvec result.
            pctrl.imbalance_tols[i] = base_ubfactor.powf(inv_log) as Real + 0.0000499;
        }
    }

    // Compute target_part_weights: uniform 1/nparts for each partition and constraint
    pctrl.target_part_weights = vec![1.0 / nparts as Real; (nparts * graph.num_constraints) as usize];

    // In C METIS, METIS_PartGraphRecursive calls SetupCtrl which resets the
    // global RNG to the original seed. pctrl already has the correct seed from
    // Control::new (via options[9] = ctrl.seed). No override needed.

    // Always initialize labels to identity for recursive bisection.
    // contract.rs sets labels to vec![0; n], so we must overwrite unconditionally.
    graph.label = (0..graph.num_vertices).collect();

    // Run recursive bisection: MlevelRecursiveBisection
    let mut part = vec![0 as Idx; num_vertices];
    let target_part_weights = pctrl.target_part_weights.clone();
    mlevel_recursive_bisection(&mut pctrl, graph, nparts, &mut part, &target_part_weights, 0);

    // In C METIS, the global RNG state after recursive bisection persists.
    // Propagate pctrl's RNG state back to the parent ctrl so RefineKWay
    // sees the same RNG state as C METIS.
    ctrl.rng = pctrl.rng;

    // Copy result into graph.partition
    graph.partition.copy_from_slice(&part);
}

/// Recursively bisect graph into nparts partitions.
///
/// Replicates MlevelRecursiveBisection from pmetis.c:
/// 1. Compute target_part_weights2[i] = sum of first half of target_part_weights for each constraint
/// 2. Call MultilevelBisect(ctrl, graph, target_part_weights2)
/// 3. Map partition back using label[i]: part[label[i]] = where[i] + fpart
/// 4. If nparts > 2: SplitGraphPart, then recurse on each subgraph
/// 5. Scale target_part_weights for recursive calls
fn mlevel_recursive_bisection(
    ctrl: &mut Control,
    graph: &mut GraphData,
    nparts: Idx,
    part: &mut [Idx],
    target_part_weights: &[Real],
    fpart: Idx,
) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Base case: single partition
    if nparts == 1 {
        // Map using label: part[label[i]] = fpart
        for i in 0..num_vertices {
            let label_i = graph.label[i] as usize;
            if label_i < part.len() {
                part[label_i] = fpart;
            }
        }
        return;
    }

    // Step 1: Compute target_part_weights2 for the bisection.
    // target_part_weights2[0*ncon+j] = sum of target_part_weights[i*ncon+j] for i in 0..nparts/2
    // target_part_weights2[1*ncon+j] = sum of target_part_weights[i*ncon+j] for i in nparts/2..nparts
    let nparts_left = nparts / 2;
    let nparts_right = nparts - nparts_left;

    let mut target_part_weights2 = vec![0.0 as Real; 2 * ncon];
    for j in 0..ncon {
        for i in 0..nparts_left as usize {
            target_part_weights2[j] += target_part_weights[i * ncon + j];
        }
        for i in nparts_left as usize..nparts as usize {
            target_part_weights2[ncon + j] += target_part_weights[i * ncon + j];
        }
    }

    // Step 2: MultilevelBisect
    multilevel_bisect(ctrl, graph, &target_part_weights2);

    // Step 3: Map partition back using label.
    // part[label[i]] = where[i] + fpart
    for i in 0..num_vertices {
        let label_i = graph.label[i] as usize;
        if label_i < part.len() {
            part[label_i] = graph.partition[i] + fpart;
        }
    }

    // Step 4: If nparts > 2, split and recurse
    if nparts > 2 {
        let (mut lgraph, mut rgraph) = split::split_graph_part(graph);

        // Step 5: Scale target_part_weights for recursive calls.
        // Left subgraph gets target_part_weights[0..nparts_left], normalized by sum(target_part_weights2[0..ncon])
        // Right subgraph gets target_part_weights[nparts_left..nparts], normalized by sum(target_part_weights2[ncon..2*ncon])
        let mut left_target_part_weights = vec![0.0 as Real; (nparts_left * graph.num_constraints) as usize];
        let mut right_target_part_weights = vec![0.0 as Real; (nparts_right * graph.num_constraints) as usize];

        for j in 0..ncon {
            let wgt_left = if target_part_weights2[j] > 0.0 { 1.0 / target_part_weights2[j] } else { 0.0 };
            let wgt_right = if target_part_weights2[ncon + j] > 0.0 { 1.0 / target_part_weights2[ncon + j] } else { 0.0 };

            for i in 0..nparts_left as usize {
                left_target_part_weights[i * ncon + j] = target_part_weights[i * ncon + j] * wgt_left;
            }
            for i in 0..nparts_right as usize {
                right_target_part_weights[i * ncon + j] =
                    target_part_weights[(nparts_left as usize + i) * ncon + j] * wgt_right;
            }
        }

        // Recurse on left subgraph
        mlevel_recursive_bisection(
            ctrl,
            &mut lgraph,
            nparts_left,
            part,
            &left_target_part_weights,
            fpart,
        );

        // Recurse on right subgraph
        mlevel_recursive_bisection(
            ctrl,
            &mut rgraph,
            nparts_right,
            part,
            &right_target_part_weights,
            fpart + nparts_left,
        );
    }
}

/// Multilevel bisection: coarsen -> init partition -> refine.
///
/// Replicates MultilevelBisect from pmetis.c:
/// 1. Setup2WayBalMultipliers(ctrl, graph, target_part_weights)
/// 2. For ncuts iterations:
///    a. cgraph = CoarsenGraph(ctrl, graph)
///    b. niparts = (cgraph->num_vertices <= ctrl->CoarsenTo) ? SMALL_NUM_INIT_PARTS : LARGE_NUM_INIT_PARTS
///    c. Init2WayPartition(ctrl, cgraph, target_part_weights, niparts)
///    d. Refine2Way(ctrl, graph, cgraph, target_part_weights) - uncoarsens and refines
///    e. Keep best based on balance/cut
/// Matches C METIS MultilevelBisect from pmetis.c:
/// - Uses ComputeLoadImbalanceDiff for balance comparison
/// - Selection: i==0 || (curbal<=0.0005 && bestobj>curobj) || (bestbal>0.0005 && curbal<bestbal)
/// - Early break on bestobj==0
fn multilevel_bisect(ctrl: &mut Control, graph: &mut GraphData, target_part_weights: &[Real]) {
    let ncuts = ctrl.num_cuts;

    // Step 1: Setup2WayBalMultipliers
    ctrl.setup_2way_balance_multipliers(&graph.inv_total_vertex_weight, target_part_weights);

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
        let mut levels = coarsen_graph(ctrl, &mut work);

        // Step 2b: Compute niparts based on coarsest graph size
        let niparts = if levels.last().unwrap().num_vertices <= ctrl.coarsen_to {
            SMALL_NUM_INIT_PARTS
        } else {
            LARGE_NUM_INIT_PARTS
        };

        // Step 2c: Init2WayPartition on the coarsest graph
        initpart::init_2way_partition(ctrl, levels.last_mut().unwrap(), target_part_weights, niparts);

        // Step 2d: Refine2Way - uncoarsen and refine through the chain.
        refine2way::refine_2way(ctrl, &mut work, &mut levels, target_part_weights);

        // After refinement, work (orggraph) has the final partition result.
        let curobj = work.edge_cut;
        let curbal = compute_load_imbalance_diff_2way(&work, &ctrl.partition_ij_balance_multipliers, &ctrl.imbalance_tols);

        // Step 2e: Keep best (matching C METIS criterion exactly)
        // C METIS only saves bestwhere on non-last iterations (i < ncuts-1).
        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            bestobj = curobj;
            bestbal = curbal;
            if i < ncuts - 1 {
                best_where = work.partition.clone();
            }
        }

        // Always track last iteration's result
        last_curobj = curobj;
        last_where = work.partition;

        if bestobj == 0 {
            break;
        }
    }

    // Post-loop restoration: match C METIS exactly.
    // C works on same graph, so graph->where has last iteration's result.
    // It only restores bestwhere if bestobj != curobj (last iter's obj).
    if bestobj != last_curobj {
        graph.partition = best_where;
    } else {
        graph.partition = last_where;
    }
    graph.edge_cut = bestobj;
}

/// ComputeLoadImbalanceDiff for 2-way partition.
/// Returns max over all partitions and constraints of:
///   part_weights[i*ncon+j] * partition_ij_balance_multipliers[i*ncon+j] - imbalance_tols[j]
/// Negative means balanced, positive means imbalanced.
fn compute_load_imbalance_diff_2way(graph: &GraphData, partition_ij_balance_multipliers: &[Real], imbalance_tols: &[Real]) -> Real {
    let ncon = graph.num_constraints as usize;
    let mut max_diff: Real = -1.0 * ncon as Real;
    for i in 0..2 {
        for j in 0..ncon {
            let idx = i * ncon + j;
            if idx < graph.part_weights.len() && idx < partition_ij_balance_multipliers.len() && j < imbalance_tols.len() {
                let diff = graph.part_weights[idx] as Real * partition_ij_balance_multipliers[idx] - imbalance_tols[j];
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
    ng.num_vertices = g.num_vertices;
    ng.num_edges = g.num_edges;
    ng.num_constraints = g.num_constraints;
    ng.xadj = g.xadj.clone();
    ng.adjacency = g.adjacency.clone();
    ng.vertex_weights = g.vertex_weights.clone();
    ng.vertex_sizes = g.vertex_sizes.clone();
    ng.edge_weights = g.edge_weights.clone();
    ng.total_vertex_weight = g.total_vertex_weight.clone();
    ng.inv_total_vertex_weight = g.inv_total_vertex_weight.clone();
    ng.label = (0..g.num_vertices).collect();
    ng
}

