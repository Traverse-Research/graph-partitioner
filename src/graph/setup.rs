use crate::types::{Idx, Real};
use crate::graph::GraphData;

/// Create an internal graph from CSR arrays.
pub fn setup_graph(
    ncon: Idx,
    xadj: &[Idx],
    adjncy: &[Idx],
    vwgt: Option<&[Idx]>,
    vsize: Option<&[Idx]>,
    adjwgt: Option<&[Idx]>,
) -> GraphData {
    let nvtxs = (xadj.len() - 1) as Idx;
    let nedges = adjncy.len() as Idx;

    let mut g = GraphData::new();
    g.nvtxs = nvtxs;
    g.nedges = nedges;
    g.ncon = ncon;

    g.xadj = xadj.to_vec();
    g.adjncy = adjncy.to_vec();

    g.vwgt = match vwgt {
        Some(w) => w.to_vec(),
        None => vec![1; (nvtxs * ncon) as usize],
    };

    g.vsize = match vsize {
        Some(s) => s.to_vec(),
        None => vec![1; nvtxs as usize],
    };

    g.adjwgt = match adjwgt {
        Some(w) => w.to_vec(),
        None => vec![1; nedges as usize],
    };

    g.label = (0..nvtxs).collect();

    setup_graph_tvwgt(&mut g);

    g
}

/// Compute total vertex weights and their inverses.
pub fn setup_graph_tvwgt(g: &mut GraphData) {
    let ncon = g.ncon as usize;
    let nvtxs = g.nvtxs as usize;

    g.tvwgt = vec![0; ncon];
    g.invtvwgt = vec![0.0; ncon];

    for j in 0..ncon {
        for i in 0..nvtxs {
            g.tvwgt[j] += g.vwgt[i * ncon + j];
        }
        g.invtvwgt[j] = 1.0 / g.tvwgt[j].max(1) as Real;
    }
}
