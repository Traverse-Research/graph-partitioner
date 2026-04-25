use crate::graph::GraphData;
use crate::types::{Idx, Real};

/// Create an internal graph from CSR arrays.
pub fn setup_graph(
    ncon: Idx,
    xadj: &[Idx],
    adjacency: &[Idx],
    vertex_weights: Option<&[Idx]>,
    vertex_sizes: Option<&[Idx]>,
    edge_weights: Option<&[Idx]>,
) -> GraphData {
    let num_vertices = (xadj.len() - 1) as Idx;
    let num_edges = adjacency.len() as Idx;

    let mut g = GraphData::new();
    g.num_vertices = num_vertices;
    g.num_edges = num_edges;
    g.num_constraints = ncon;

    g.xadj = xadj.to_vec();
    g.adjacency = adjacency.to_vec();

    g.vertex_weights = match vertex_weights {
        Some(w) => w.to_vec(),
        None => vec![1; (num_vertices * ncon) as usize],
    };

    g.vertex_sizes = match vertex_sizes {
        Some(s) => s.to_vec(),
        None => vec![1; num_vertices as usize],
    };

    g.edge_weights = match edge_weights {
        Some(w) => w.to_vec(),
        None => vec![1; num_edges as usize],
    };

    g.label = (0..num_vertices).collect();

    setup_graph_total_vertex_weight(&mut g);

    g
}

/// Compute total vertex weights and their inverses.
pub fn setup_graph_total_vertex_weight(g: &mut GraphData) {
    let ncon = g.num_constraints as usize;
    let num_vertices = g.num_vertices as usize;

    g.total_vertex_weight = vec![0; ncon];
    g.inv_total_vertex_weight = vec![0.0; ncon];

    for j in 0..ncon {
        for i in 0..num_vertices {
            g.total_vertex_weight[j] += g.vertex_weights[i * ncon + j];
        }
        g.inv_total_vertex_weight[j] = 1.0 / g.total_vertex_weight[j].max(1) as Real;
    }
}
