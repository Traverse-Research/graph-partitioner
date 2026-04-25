use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::types::Idx;

#[expect(dead_code)]
/// Compute the subdomain graph (partition adjacency graph).
pub fn compute_subdomain_graph(
    graph: &GraphData,
    partition: &[Idx],
    nparts: Idx,
) -> (Vec<Vec<Idx>>, Vec<Vec<Idx>>) {
    let num_vertices = graph.num_vertices as usize;
    let np = nparts as usize;

    let mut sadj: Vec<Vec<Idx>> = vec![Vec::new(); np];
    let mut swgt: Vec<Vec<Idx>> = vec![Vec::new(); np];

    let mut marker = vec![-1 as Idx; np];

    for i in 0..num_vertices {
        let me = partition[i] as usize;
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let other = partition[graph.adjacency[k] as usize] as usize;
            if other != me {
                if marker[other] == -1 {
                    marker[other] = sadj[me].len() as Idx;
                    sadj[me].push(other as Idx);
                    swgt[me].push(graph.edge_weights[k]);
                } else {
                    swgt[me][marker[other] as usize] += graph.edge_weights[k];
                }
            }
        }
        // Clean up marker
        for &s in &sadj[me] {
            marker[s as usize] = -1;
        }
    }

    (sadj, swgt)
}

#[expect(dead_code)]
/// Try to minimize the maximum degree in the subdomain graph.
pub fn eliminate_subdomain_edges(_ctrl: &mut Control, _graph: &mut GraphData, _nparts: Idx) {
    // This is a complex optimization that tries to reduce the
    // connectivity between partitions. For now, this is a no-op
    // as it's an optimization that doesn't affect correctness.
}
