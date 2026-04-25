use crate::ctrl::Control;
use crate::graph::setup::setup_graph;
use crate::partition;
use crate::types::{Idx, Result};
use crate::Graph;

/// Entry point for METIS_PartGraphKway logic.
///
/// Matches C METIS METIS_PartGraphKway from kmetis.c:
/// 1. SetupCtrl with METIS_OP_KMETIS settings
/// 2. SetupGraph
/// 3. SetupKWayBalMultipliers
/// 4. Compute CoarsenTo and nIparts
/// 5. MlevelKWayPartitioning
pub fn part_kway(graph: Graph, part: &mut [Idx]) -> Result<Idx> {
    let ncon = graph.num_constraints;
    let nparts = graph.num_parts;

    if nparts == 1 {
        for p in part.iter_mut() {
            *p = 0;
        }
        return Ok(0);
    }

    // Set up the control structure (is_kway = true)
    let mut ctrl = Control::new(&graph.options, ncon, nparts, true);

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

    // SetupKWayBalMultipliers
    let inv_total_vertex_weight = gdata.inv_total_vertex_weight.clone();
    ctrl.setup_kway_balance_multipliers(&inv_total_vertex_weight);

    // Compute CoarsenTo = max(num_vertices/(40*log2(nparts)), 30*nparts)
    let log2_nparts = (nparts as f64).log2().max(1.0) as Idx;
    let coarsen_to = ((gdata.num_vertices / (40 * log2_nparts)).max(30 * nparts)).max(1);
    ctrl.coarsen_to = coarsen_to;

    // Compute nIparts = (CoarsenTo == 30*nparts ? 4 : 5)
    ctrl.num_init_parts = if coarsen_to == 30 * nparts { 4 } else { 5 };

    // Run multilevel k-way partitioning
    let cut = partition::mlevel_kway_partitioning(&mut ctrl, &mut gdata, part);

    Ok(cut)
}
