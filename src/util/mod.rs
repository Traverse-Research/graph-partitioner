pub mod pqueue;

use crate::types::{Idx, Real};
use crate::graph::GraphData;

#[expect(dead_code)]
/// Return the index of the maximum element.
pub fn iargmax(n: usize, arr: &[Idx]) -> usize {
    let mut max_idx = 0;
    let mut max_val = arr[0];
    for i in 1..n {
        if arr[i] > max_val {
            max_val = arr[i];
            max_idx = i;
        }
    }
    max_idx
}

#[expect(dead_code)]
/// Return the index of the maximum element (Real).
pub fn rargmax(n: usize, arr: &[Real]) -> usize {
    let mut max_idx = 0;
    let mut max_val = arr[0];
    for i in 1..n {
        if arr[i] > max_val {
            max_val = arr[i];
            max_idx = i;
        }
    }
    max_idx
}

#[expect(dead_code)]
/// Compute the edge-cut of a partitioning.
pub fn compute_cut(graph: &GraphData, partition: &[Idx]) -> Idx {
    let num_vertices = graph.num_vertices as usize;
    let mut cut = 0;
    for i in 0..num_vertices {
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            if partition[i] != partition[graph.adjacency[k] as usize] {
                cut += graph.edge_weights[k];
            }
        }
    }
    cut / 2
}

#[expect(dead_code)]
/// Compute the total communication volume of a partitioning.
pub fn compute_volume(graph: &GraphData, partition: &[Idx]) -> Idx {
    let num_vertices = graph.num_vertices as usize;
    let nparts = partition.iter().copied().max().unwrap_or(0) + 1;
    let mut vol = 0;

    let mut marker = vec![false; nparts as usize];

    for i in 0..num_vertices {
        let me = partition[i] as usize;
        let mut cnt = 0;

        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            let other = partition[graph.adjacency[k] as usize] as usize;
            if other != me && !marker[other] {
                marker[other] = true;
                cnt += 1;
            }
        }

        // Clean up marker
        for k in graph.xadj[i] as usize..graph.xadj[i + 1] as usize {
            marker[partition[graph.adjacency[k] as usize] as usize] = false;
        }

        if cnt > 0 {
            vol += graph.vertex_sizes[i] * cnt as Idx;
        }
    }

    vol
}

#[expect(dead_code)]
/// Bucket sort keys in increasing order, returning permutation array.
pub fn bucket_sort_keys_inc(n: usize, max: Idx, keys: &[Idx], tperm: &mut [Idx]) {
    if n == 0 {
        return;
    }

    let max_usize = max as usize + 1;
    let mut counts = vec![0usize; max_usize];

    for i in 0..n {
        counts[keys[i] as usize] += 1;
    }

    // Prefix sum
    let mut sum = 0;
    for i in 0..max_usize {
        let c = counts[i];
        counts[i] = sum;
        sum += c;
    }

    // Place elements
    for i in 0..n {
        let k = keys[i] as usize;
        tperm[counts[k]] = i as Idx;
        counts[k] += 1;
    }
}

#[expect(dead_code)]
/// Compute partition weights.
pub fn compute_partition_weights(graph: &GraphData, partition: &[Idx], nparts: Idx) -> Vec<Idx> {
    let ncon = graph.num_constraints as usize;
    let mut part_weights = vec![0 as Idx; nparts as usize * ncon];
    for i in 0..graph.num_vertices as usize {
        let p = partition[i] as usize;
        for j in 0..ncon {
            part_weights[p * ncon + j] += graph.vertex_weights[i * ncon + j];
        }
    }
    part_weights
}

#[expect(dead_code)]
/// Check if a 2-way partition is balanced.
pub fn is_balanced_2way(
    total_vertex_weight_j: Idx,
    pwgt0_j: Idx,
    pwgt1_j: Idx,
    target_ratio: Real,
    ubfactor: Real,
) -> bool {
    let target0 = (total_vertex_weight_j as Real * target_ratio) as Idx;
    let target1 = total_vertex_weight_j - target0;
    let ub = ubfactor;

    pwgt0_j as Real <= ub * target0 as Real && pwgt1_j as Real <= ub * target1 as Real
}
