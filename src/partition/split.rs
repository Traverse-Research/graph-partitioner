use crate::types::Idx;
use crate::graph::GraphData;
use crate::graph::setup::setup_graph_total_vertex_weight;

/// Split a bisected graph into two subgraphs based on partition[i] == 0 or 1.
pub fn split_graph_part(graph: &GraphData) -> (GraphData, GraphData) {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;

    // Build vertex maps
    let mut map = vec![0 as Idx; num_vertices]; // fine vertex -> subgraph vertex
    let mut lcount = 0;
    let mut rcount = 0;

    for i in 0..num_vertices {
        if graph.partition[i] == 0 {
            map[i] = lcount;
            lcount += 1;
        } else {
            map[i] = rcount;
            rcount += 1;
        }
    }

    // Build left subgraph
    let mut lxadj = vec![0 as Idx; lcount as usize + 1];
    let mut ladjacency = Vec::new();
    let mut ledge_weights = Vec::new();
    let mut lvertex_weights = vec![0 as Idx; lcount as usize * ncon];
    let mut lvertex_sizes = vec![0 as Idx; lcount as usize];
    let mut llabel = vec![0 as Idx; lcount as usize];

    let mut rxadj = vec![0 as Idx; rcount as usize + 1];
    let mut radjacency = Vec::new();
    let mut redge_weights = Vec::new();
    let mut rvertex_weights = vec![0 as Idx; rcount as usize * ncon];
    let mut rvertex_sizes = vec![0 as Idx; rcount as usize];
    let mut rlabel = vec![0 as Idx; rcount as usize];

    let mut li = 0usize;
    let mut ri = 0usize;

    for i in 0..num_vertices {
        if graph.partition[i] == 0 {
            for j in 0..ncon {
                lvertex_weights[li * ncon + j] = graph.vertex_weights[i * ncon + j];
            }
            lvertex_sizes[li] = graph.vertex_sizes[i];
            llabel[li] = graph.label[i];

            for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let nbr = graph.adjacency[k] as usize;
                if graph.partition[nbr] == 0 {
                    ladjacency.push(map[nbr]);
                    ledge_weights.push(graph.edge_weights[k]);
                }
            }
            li += 1;
            lxadj[li] = ladjacency.len() as Idx;
        } else {
            for j in 0..ncon {
                rvertex_weights[ri * ncon + j] = graph.vertex_weights[i * ncon + j];
            }
            rvertex_sizes[ri] = graph.vertex_sizes[i];
            rlabel[ri] = graph.label[i];

            for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let nbr = graph.adjacency[k] as usize;
                if graph.partition[nbr] == 1 {
                    radjacency.push(map[nbr]);
                    redge_weights.push(graph.edge_weights[k]);
                }
            }
            ri += 1;
            rxadj[ri] = radjacency.len() as Idx;
        }
    }

    let mut lgraph = GraphData::new();
    lgraph.num_vertices = lcount;
    lgraph.num_edges = ladjacency.len() as Idx;
    lgraph.num_constraints = graph.num_constraints;
    lgraph.xadj = lxadj;
    lgraph.adjacency = ladjacency;
    lgraph.vertex_weights = lvertex_weights;
    lgraph.vertex_sizes = lvertex_sizes;
    lgraph.edge_weights = ledge_weights;
    lgraph.label = llabel;
    lgraph.partition = vec![0; lcount as usize];
    setup_graph_total_vertex_weight(&mut lgraph);

    let mut rgraph = GraphData::new();
    rgraph.num_vertices = rcount;
    rgraph.num_edges = radjacency.len() as Idx;
    rgraph.num_constraints = graph.num_constraints;
    rgraph.xadj = rxadj;
    rgraph.adjacency = radjacency;
    rgraph.vertex_weights = rvertex_weights;
    rgraph.vertex_sizes = rvertex_sizes;
    rgraph.edge_weights = redge_weights;
    rgraph.label = rlabel;
    rgraph.partition = vec![0; rcount as usize];
    setup_graph_total_vertex_weight(&mut rgraph);

    (lgraph, rgraph)
}
