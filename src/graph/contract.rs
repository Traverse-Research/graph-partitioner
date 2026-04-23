use crate::types::Idx;
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::graph::setup::setup_graph_tvwgt;

const HTLENGTH: usize = (1 << 13) - 1; // 8191

/// Build a coarsened graph from the matching stored in graph.cmap and graph.match_.
pub fn create_coarse_graph(_ctrl: &mut Ctrl, graph: &GraphData, cnvtxs: Idx) -> Box<GraphData> {
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    let cnvtxs_usize = cnvtxs as usize;
    let mask = HTLENGTH;

    let mut cg = Box::new(GraphData::new());
    cg.nvtxs = cnvtxs;
    cg.ncon = graph.ncon;

    let mut cxadj = vec![0 as Idx; cnvtxs_usize + 1];
    let mut cvwgt = vec![0 as Idx; cnvtxs_usize * ncon];
    let mut cvsize = vec![0 as Idx; cnvtxs_usize];

    // Pre-allocate adjacency (will be trimmed)
    let max_edges = graph.nedges as usize + 1;
    let mut cadjncy = vec![0 as Idx; max_edges];
    let mut cadjwgt = vec![0 as Idx; max_edges];

    // Hash table for edge merging
    let mut htable = vec![-1 as i32; mask + 1];

    let mut cnedges: usize = 0;

    for v in 0..nvtxs {
        let u = graph.match_[v] as usize;
        if u < v {
            continue; // Already processed as part of the pair
        }

        let cv = graph.cmap[v] as usize;

        // Accumulate vertex weights
        for j in 0..ncon {
            cvwgt[cv * ncon + j] += graph.vwgt[v * ncon + j];
        }
        cvsize[cv] += graph.vsize[v];
        if v != u {
            for j in 0..ncon {
                cvwgt[cv * ncon + j] += graph.vwgt[u * ncon + j];
            }
            cvsize[cv] += graph.vsize[u];
        }

        // Edge merging using hash table
        let edge_start = cnedges;

        // Seed with self-loop sentinel
        let hk = cv & mask;
        htable[hk] = cnedges as i32;
        cadjncy[cnedges] = cv as Idx;
        cadjwgt[cnedges] = 0;
        cnedges += 1;

        // Process edges of v
        for jj in graph.xadj[v] as usize..graph.xadj[v + 1] as usize {
            let k = graph.cmap[graph.adjncy[jj] as usize] as usize;
            let mut kk = k & mask;

            // Linear probe to find k or empty slot
            loop {
                if htable[kk] == -1 {
                    // New neighbor
                    htable[kk] = cnedges as i32;
                    cadjncy[cnedges] = k as Idx;
                    cadjwgt[cnedges] = graph.adjwgt[jj];
                    cnedges += 1;
                    break;
                } else if cadjncy[htable[kk] as usize] == k as Idx {
                    // Existing neighbor, accumulate weight
                    cadjwgt[htable[kk] as usize] += graph.adjwgt[jj];
                    break;
                }
                kk = (kk + 1) & mask;
            }
        }

        // Process edges of u (if matched pair)
        if v != u {
            for jj in graph.xadj[u] as usize..graph.xadj[u + 1] as usize {
                let k = graph.cmap[graph.adjncy[jj] as usize] as usize;
                let mut kk = k & mask;

                loop {
                    if htable[kk] == -1 {
                        htable[kk] = cnedges as i32;
                        cadjncy[cnedges] = k as Idx;
                        cadjwgt[cnedges] = graph.adjwgt[jj];
                        cnedges += 1;
                        break;
                    } else if cadjncy[htable[kk] as usize] == k as Idx {
                        cadjwgt[htable[kk] as usize] += graph.adjwgt[jj];
                        break;
                    }
                    kk = (kk + 1) & mask;
                }
            }
        }

        // Clean up hash table (LIFO order for linear-probe correctness)
        for j in (edge_start..cnedges).rev() {
            let k = cadjncy[j] as usize;
            let mut kk = k & mask;
            loop {
                if htable[kk] != -1 && cadjncy[htable[kk] as usize] as usize == k {
                    htable[kk] = -1;
                    break;
                }
                kk = (kk + 1) & mask;
            }
        }

        // Remove self-loop (at edge_start position)
        cnedges -= 1;
        if edge_start < cnedges {
            cadjncy[edge_start] = cadjncy[cnedges];
            cadjwgt[edge_start] = cadjwgt[cnedges];
        }
        cxadj[cv + 1] = cnedges as Idx;
    }

    cadjncy.truncate(cnedges);
    cadjwgt.truncate(cnedges);

    cg.xadj = cxadj;
    cg.adjncy = cadjncy;
    cg.adjwgt = cadjwgt;
    cg.vwgt = cvwgt;
    cg.vsize = cvsize;
    cg.nedges = cnedges as Idx;
    cg.label = vec![0; cnvtxs_usize];

    setup_graph_tvwgt(&mut cg);

    cg
}
