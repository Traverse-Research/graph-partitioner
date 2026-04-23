use crate::types::Idx;
use crate::graph::GraphData;
use crate::graph::setup::setup_graph_tvwgt;

/// Split a bisected graph into two subgraphs based on where_[i] == 0 or 1.
pub fn split_graph_part(graph: &GraphData) -> (GraphData, GraphData) {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    // Build vertex maps
    let mut map = vec![0 as Idx; nvtxs]; // fine vertex -> subgraph vertex
    let mut lcount = 0;
    let mut rcount = 0;

    for i in 0..nvtxs {
        if graph.where_[i] == 0 {
            map[i] = lcount;
            lcount += 1;
        } else {
            map[i] = rcount;
            rcount += 1;
        }
    }

    // Build left subgraph
    let mut lxadj = vec![0 as Idx; lcount as usize + 1];
    let mut ladjncy = Vec::new();
    let mut ladjwgt = Vec::new();
    let mut lvwgt = vec![0 as Idx; lcount as usize * ncon];
    let mut lvsize = vec![0 as Idx; lcount as usize];
    let mut llabel = vec![0 as Idx; lcount as usize];

    let mut rxadj = vec![0 as Idx; rcount as usize + 1];
    let mut radjncy = Vec::new();
    let mut radjwgt = Vec::new();
    let mut rvwgt = vec![0 as Idx; rcount as usize * ncon];
    let mut rvsize = vec![0 as Idx; rcount as usize];
    let mut rlabel = vec![0 as Idx; rcount as usize];

    let mut li = 0usize;
    let mut ri = 0usize;

    for i in 0..nvtxs {
        if graph.where_[i] == 0 {
            for j in 0..ncon {
                lvwgt[li * ncon + j] = graph.vwgt[i * ncon + j];
            }
            lvsize[li] = graph.vsize[i];
            llabel[li] = graph.label[i];

            for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let nbr = graph.adjncy[k] as usize;
                if graph.where_[nbr] == 0 {
                    ladjncy.push(map[nbr]);
                    ladjwgt.push(graph.adjwgt[k]);
                }
            }
            li += 1;
            lxadj[li] = ladjncy.len() as Idx;
        } else {
            for j in 0..ncon {
                rvwgt[ri * ncon + j] = graph.vwgt[i * ncon + j];
            }
            rvsize[ri] = graph.vsize[i];
            rlabel[ri] = graph.label[i];

            for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
                let nbr = graph.adjncy[k] as usize;
                if graph.where_[nbr] == 1 {
                    radjncy.push(map[nbr]);
                    radjwgt.push(graph.adjwgt[k]);
                }
            }
            ri += 1;
            rxadj[ri] = radjncy.len() as Idx;
        }
    }

    let mut lgraph = GraphData::new();
    lgraph.nvtxs = lcount;
    lgraph.nedges = ladjncy.len() as Idx;
    lgraph.ncon = graph.ncon;
    lgraph.xadj = lxadj;
    lgraph.adjncy = ladjncy;
    lgraph.vwgt = lvwgt;
    lgraph.vsize = lvsize;
    lgraph.adjwgt = ladjwgt;
    lgraph.label = llabel;
    lgraph.where_ = vec![0; lcount as usize];
    setup_graph_tvwgt(&mut lgraph);

    let mut rgraph = GraphData::new();
    rgraph.nvtxs = rcount;
    rgraph.nedges = radjncy.len() as Idx;
    rgraph.ncon = graph.ncon;
    rgraph.xadj = rxadj;
    rgraph.adjncy = radjncy;
    rgraph.vwgt = rvwgt;
    rgraph.vsize = rvsize;
    rgraph.adjwgt = radjwgt;
    rgraph.label = rlabel;
    rgraph.where_ = vec![0; rcount as usize];
    setup_graph_tvwgt(&mut rgraph);

    (lgraph, rgraph)
}
