use crate::ctrl::Control;
use crate::graph::GraphData;
use crate::types::Idx;

/// Find connected components within each partition.
pub fn find_partition_induced_components(
    graph: &GraphData,
    partition: &[Idx],
    nparts: Idx,
) -> (Vec<Vec<Vec<usize>>>,) {
    let num_vertices = graph.num_vertices as usize;
    let np = nparts as usize;

    let mut components: Vec<Vec<Vec<usize>>> = vec![Vec::new(); np];
    let mut visited = vec![false; num_vertices];

    for p in 0..np {
        for i in 0..num_vertices {
            if partition[i] as usize == p && !visited[i] {
                // BFS to find component
                let mut component = Vec::new();
                let mut queue = vec![i];
                visited[i] = true;

                while let Some(v) = queue.pop() {
                    component.push(v);
                    for k in graph.xadj[v] as usize..graph.xadj[v + 1] as usize {
                        let u = graph.adjacency[k] as usize;
                        if !visited[u] && partition[u] as usize == p {
                            visited[u] = true;
                            queue.push(u);
                        }
                    }
                }

                components[p].push(component);
            }
        }
    }

    (components,)
}

#[expect(dead_code)]
/// Eliminate disconnected components by reassigning them to neighboring partitions.
pub fn eliminate_components(_ctrl: &mut Control, graph: &mut GraphData, nparts: Idx) {
    let (components,) = find_partition_induced_components(graph, &graph.partition.clone(), nparts);

    for p in 0..nparts as usize {
        if components[p].len() <= 1 {
            continue;
        }

        // Keep the largest component, move others
        let mut largest_idx = 0;
        let mut largest_size = 0;
        for (ci, comp) in components[p].iter().enumerate() {
            if comp.len() > largest_size {
                largest_size = comp.len();
                largest_idx = ci;
            }
        }

        for (ci, comp) in components[p].iter().enumerate() {
            if ci == largest_idx {
                continue;
            }

            // Find the best neighboring partition for this component
            for &v in comp {
                let mut best_part = p as Idx;
                let mut best_conn = 0;

                for k in graph.xadj[v] as usize..graph.xadj[v + 1] as usize {
                    let u = graph.adjacency[k] as usize;
                    let up = graph.partition[u];
                    if up as usize != p {
                        // Count connection
                        if graph.edge_weights[k] > best_conn {
                            best_conn = graph.edge_weights[k];
                            best_part = up;
                        }
                    }
                }

                if best_part as usize != p {
                    graph.partition[v] = best_part;
                }
            }
        }
    }
}
