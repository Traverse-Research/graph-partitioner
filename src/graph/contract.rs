use crate::types::Idx;
use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::graph::setup::setup_graph_total_vertex_weight;

const HTLENGTH: usize = (1 << 13) - 1; // 8191

/// Build a coarsened graph from the matching stored in graph.coarse_map and graph.matching.
pub fn create_coarse_graph(_ctrl: &mut Control, graph: &GraphData, cnum_vertices: Idx) -> GraphData {
    let num_vertices = graph.num_vertices as usize;
    let ncon = graph.num_constraints as usize;
    let cnum_vertices_usize = cnum_vertices as usize;
    let mask = HTLENGTH;

    let mut cg = GraphData::new();
    cg.num_vertices = cnum_vertices;
    cg.num_constraints = graph.num_constraints;

    let mut cxadj = vec![0 as Idx; cnum_vertices_usize + 1];
    let mut cvertex_weights = vec![0 as Idx; cnum_vertices_usize * ncon];
    let mut cvertex_sizes = vec![0 as Idx; cnum_vertices_usize];

    // Pre-allocate adjacency (will be trimmed)
    let max_edges = graph.num_edges as usize + 1;
    let mut cadjacency = vec![0 as Idx; max_edges];
    let mut cedge_weights = vec![0 as Idx; max_edges];

    // Hash table for edge merging
    let mut htable = vec![-1 as i32; mask + 1];

    let mut cnum_edges: usize = 0;

    for v in 0..num_vertices {
        let u = graph.matching[v] as usize;
        if u < v {
            continue; // Already processed as part of the pair
        }

        let cv = graph.coarse_map[v] as usize;

        // Accumulate vertex weights
        for j in 0..ncon {
            cvertex_weights[cv * ncon + j] += graph.vertex_weights[v * ncon + j];
        }
        cvertex_sizes[cv] += graph.vertex_sizes[v];
        if v != u {
            for j in 0..ncon {
                cvertex_weights[cv * ncon + j] += graph.vertex_weights[u * ncon + j];
            }
            cvertex_sizes[cv] += graph.vertex_sizes[u];
        }

        // Edge merging using hash table
        let edge_start = cnum_edges;

        // Seed with self-loop sentinel
        let hk = cv & mask;
        htable[hk] = cnum_edges as i32;
        cadjacency[cnum_edges] = cv as Idx;
        cedge_weights[cnum_edges] = 0;
        cnum_edges += 1;

        // Process edges of v
        for jj in graph.xadj[v] as usize..graph.xadj[v + 1] as usize {
            let k = graph.coarse_map[graph.adjacency[jj] as usize] as usize;
            let mut kk = k & mask;

            // Linear probe to find k or empty slot
            loop {
                if htable[kk] == -1 {
                    // New neighbor
                    htable[kk] = cnum_edges as i32;
                    cadjacency[cnum_edges] = k as Idx;
                    cedge_weights[cnum_edges] = graph.edge_weights[jj];
                    cnum_edges += 1;
                    break;
                } else if cadjacency[htable[kk] as usize] == k as Idx {
                    // Existing neighbor, accumulate weight
                    cedge_weights[htable[kk] as usize] += graph.edge_weights[jj];
                    break;
                }
                kk = (kk + 1) & mask;
            }
        }

        // Process edges of u (if matched pair)
        if v != u {
            for jj in graph.xadj[u] as usize..graph.xadj[u + 1] as usize {
                let k = graph.coarse_map[graph.adjacency[jj] as usize] as usize;
                let mut kk = k & mask;

                loop {
                    if htable[kk] == -1 {
                        htable[kk] = cnum_edges as i32;
                        cadjacency[cnum_edges] = k as Idx;
                        cedge_weights[cnum_edges] = graph.edge_weights[jj];
                        cnum_edges += 1;
                        break;
                    } else if cadjacency[htable[kk] as usize] == k as Idx {
                        cedge_weights[htable[kk] as usize] += graph.edge_weights[jj];
                        break;
                    }
                    kk = (kk + 1) & mask;
                }
            }
        }

        // Clean up hash table (LIFO order for linear-probe correctness)
        for j in (edge_start..cnum_edges).rev() {
            let k = cadjacency[j] as usize;
            let mut kk = k & mask;
            loop {
                if htable[kk] != -1 && cadjacency[htable[kk] as usize] as usize == k {
                    htable[kk] = -1;
                    break;
                }
                kk = (kk + 1) & mask;
            }
        }

        // Remove self-loop (at edge_start position)
        cnum_edges -= 1;
        if edge_start < cnum_edges {
            cadjacency[edge_start] = cadjacency[cnum_edges];
            cedge_weights[edge_start] = cedge_weights[cnum_edges];
        }
        cxadj[cv + 1] = cnum_edges as Idx;
    }

    cadjacency.truncate(cnum_edges);
    cedge_weights.truncate(cnum_edges);

    cg.xadj = cxadj;
    cg.adjacency = cadjacency;
    cg.edge_weights = cedge_weights;
    cg.vertex_weights = cvertex_weights;
    cg.vertex_sizes = cvertex_sizes;
    cg.num_edges = cnum_edges as Idx;
    cg.label = vec![0; cnum_vertices_usize];

    setup_graph_total_vertex_weight(&mut cg);

    cg
}
