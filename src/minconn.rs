use crate::types::Idx;
use crate::ctrl::Ctrl;
use crate::graph::GraphData;

#[allow(dead_code)]
/// Compute the subdomain graph (partition adjacency graph).
pub fn compute_subdomain_graph(
    graph: &GraphData,
    where_: &[Idx],
    nparts: Idx,
) -> (Vec<Vec<Idx>>, Vec<Vec<Idx>>) {
    let nvtxs = graph.nvtxs as usize;
    let np = nparts as usize;

    let mut sadj: Vec<Vec<Idx>> = vec![Vec::new(); np];
    let mut swgt: Vec<Vec<Idx>> = vec![Vec::new(); np];

    let mut marker = vec![-1 as Idx; np];

    for i in 0..nvtxs {
        let me = where_[i] as usize;
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let other = where_[graph.adjncy[k] as usize] as usize;
            if other != me {
                if marker[other] == -1 {
                    marker[other] = sadj[me].len() as Idx;
                    sadj[me].push(other as Idx);
                    swgt[me].push(graph.adjwgt[k]);
                } else {
                    swgt[me][marker[other] as usize] += graph.adjwgt[k];
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

#[allow(dead_code)]
/// Try to minimize the maximum degree in the subdomain graph.
pub fn eliminate_subdomain_edges(_ctrl: &mut Ctrl, _graph: &mut GraphData, _nparts: Idx) {
    // This is a complex optimization that tries to reduce the
    // connectivity between partitions. For now, this is a no-op
    // as it's an optimization that doesn't affect correctness.
}
