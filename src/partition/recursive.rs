use crate::ctrl::Control;
use crate::graph::setup::setup_graph;
use crate::partition::pmetis::mlevel_recursive_bisection;
use crate::types::{Idx, Result};
use crate::Graph;

/// Entry point for METIS_PartGraphRecursive logic.
///
/// Matches C METIS METIS_PartGraphRecursive from pmetis.c:
/// 1. SetupCtrl with METIS_OP_PMETIS settings (is_kway=false)
/// 2. SetupGraph
/// 3. Set labels to identity
/// 4. MlevelRecursiveBisection
/// 5. Compute edge cut from partition
pub fn part_recursive(graph: Graph, part: &mut [Idx]) -> Result<Idx> {
    let ncon = graph.num_constraints;
    let nparts = graph.num_parts;

    if nparts == 1 {
        for p in part.iter_mut() {
            *p = 0;
        }
        return Ok(0);
    }

    // Set up the control structure (is_kway = false => PMETIS settings)
    let mut ctrl = Control::new(&graph.options, ncon, nparts, false);

    // Handle ncon > 1 overrides matching C METIS SetupCtrl for METIS_OP_PMETIS:
    // - iptype forced to RANDOM(1) for multi-constraint
    // - CoarsenTo = 100 instead of 20
    if ncon > 1 {
        if graph.options.init_part_type.is_none() {
            ctrl.init_part_type = 1; // METIS_IPTYPE_RANDOM
        }
    }
    ctrl.coarsen_to = if ncon == 1 { 20 } else { 100 };

    // Override target_part_weights if user provided
    if let Some(tp) = graph.target_part_weights {
        ctrl.set_target_part_weights(tp);
    }

    // Override ubvec if user provided
    if let Some(ub) = graph.ubvec {
        ctrl.set_ubvec(ub);
    }

    // Create internal graph
    let mut gdata = setup_graph(
        ncon,
        graph.xadj,
        graph.adjacency,
        graph.vertex_weights,
        graph.vertex_sizes,
        graph.edge_weights,
    );

    // Labels = identity (required for recursive bisection)
    gdata.label = (0..gdata.num_vertices).collect();

    // Run recursive bisection
    let target_part_weights = ctrl.target_part_weights.clone();
    mlevel_recursive_bisection(&mut ctrl, &mut gdata, nparts, part, &target_part_weights, 0);

    // Compute edge cut from the final partition
    let edge_cut = compute_edge_cut(&gdata, part);

    Ok(edge_cut)
}

/// Compute edge cut: sum of edge weights for edges crossing partition boundaries.
/// Each edge (u,v) appears twice in CSR, so we count each direction and divide by 2.
fn compute_edge_cut(graph: &crate::graph::GraphData, part: &[Idx]) -> Idx {
    let n = graph.num_vertices as usize;
    let mut cut: Idx = 0;
    for u in 0..n {
        let start = graph.xadj[u] as usize;
        let end = graph.xadj[u + 1] as usize;
        for j in start..end {
            let v = graph.adjacency[j] as usize;
            if part[u] != part[v] {
                cut += graph.edge_weights[j];
            }
        }
    }
    cut / 2
}
