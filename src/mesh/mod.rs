pub mod meshpart;

use crate::types::Idx;

/// Build the dual graph from a mesh, matching C METIS CreateGraphDual exactly.
/// Elements become vertices. Two elements share an edge in the dual if they
/// share at least `min_common_nodes` nodes.
///
/// Adjacency order matches C METIS: neighbors are listed in the order they are
/// first discovered by iterating through the element's nodes (in element_indices order)
/// and for each node, iterating through elements sharing that node (in element
/// index order from the node-element CSR).
pub fn create_graph_dual(
    ne: Idx,
    nn: Idx,
    element_offsets: &[Idx],
    element_indices: &[Idx],
    min_common_nodes: Idx,
) -> (Vec<Idx>, Vec<Idx>) {
    let ne_usize = ne as usize;

    // Build node-to-element CSR mapping (nptr/nind)
    let (nptr, nind) = build_node_element_csr(ne, nn, element_offsets, element_indices);

    // Build dual graph using marker-based discovery (matching C METIS FindCommonElements)
    let mut xadj = vec![0 as Idx; ne_usize + 1];
    let mut adjacency_vec = Vec::new();

    // marker[j] tracks the number of shared nodes between current element and element j
    // 0 = not yet encountered
    let mut marker = vec![0i32; ne_usize];
    let mut nbrs: Vec<usize> = Vec::new();

    for i in 0..ne_usize {
        nbrs.clear();

        // For each node of element i, find all elements sharing that node
        for k in element_offsets[i] as usize..element_offsets[i + 1] as usize {
            let node = element_indices[k] as usize;
            for j_idx in nptr[node] as usize..nptr[node + 1] as usize {
                let j = nind[j_idx] as usize;
                if j == i {
                    continue;
                }
                if marker[j] == 0 {
                    // First time seeing this neighbor — record in discovery order
                    nbrs.push(j);
                }
                marker[j] += 1;
            }
        }

        // Filter: keep only neighbors with shared-node count >= min_common_nodes
        // Preserve the discovery order (matching C METIS)
        for &j in &nbrs {
            if marker[j] >= min_common_nodes {
                adjacency_vec.push(j as Idx);
            }
        }

        // Reset markers for next element
        for &j in &nbrs {
            marker[j] = 0;
        }

        xadj[i + 1] = adjacency_vec.len() as Idx;
    }

    (xadj, adjacency_vec)
}

/// Build node-to-element CSR mapping.
/// For each node, stores the list of elements containing it (in element index order).
pub fn build_node_element_csr(
    ne: Idx,
    nn: Idx,
    element_offsets: &[Idx],
    element_indices: &[Idx],
) -> (Vec<Idx>, Vec<Idx>) {
    let ne_usize = ne as usize;
    let nn_usize = nn as usize;

    let mut nptr = vec![0 as Idx; nn_usize + 1];
    for i in 0..ne_usize {
        for k in element_offsets[i] as usize..element_offsets[i + 1] as usize {
            nptr[element_indices[k] as usize + 1] += 1;
        }
    }
    // Prefix sum (MAKECSR equivalent)
    for i in 0..nn_usize {
        nptr[i + 1] += nptr[i];
    }

    let mut nind = vec![0 as Idx; nptr[nn_usize] as usize];
    let mut nptr_copy = nptr.clone();
    for i in 0..ne_usize {
        for k in element_offsets[i] as usize..element_offsets[i + 1] as usize {
            let node = element_indices[k] as usize;
            nind[nptr_copy[node] as usize] = i as Idx;
            nptr_copy[node] += 1;
        }
    }

    (nptr, nind)
}
