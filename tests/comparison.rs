mod fixtures;

use fixtures::*;

/// Helper: run part_kway on both implementations and compare.
fn compare_part_kway(xadj: &[i32], adjacency: &[i32], nparts: i32, seed: i32) {
    let n = xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    let c_cut = metis::Graph::new(1, nparts, xadj, adjacency)
        .unwrap()
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    let r_cut = metis_clone::Graph::new(1, nparts, xadj, adjacency)
        .unwrap()
        .seed(seed)
        .part_kway(&mut r_part)
        .unwrap();

    // First check: valid partitions
    for i in 0..n {
        assert!(
            r_part[i] >= 0 && r_part[i] < nparts,
            "Rust partition[{}] = {} out of range [0, {})",
            i, r_part[i], nparts
        );
    }

    // Check edge cut matches
    assert_eq!(c_cut, r_cut, "Edge cuts differ (seed={}, nparts={}): C={}, Rust={}", seed, nparts, c_cut, r_cut);

    // Check partition arrays match
    assert_eq!(c_part, r_part, "Partitions differ (seed={}, nparts={})", seed, nparts);
}

#[allow(dead_code)]
/// Helper for weighted graphs.
fn compare_part_kway_weighted(
    xadj: &[i32],
    adjacency: &[i32],
    vertex_weights: &[i32],
    edge_weights: &[i32],
    nparts: i32,
    seed: i32,
) {
    let n = xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    let c_cut = metis::Graph::new(1, nparts, xadj, adjacency)
        .unwrap()
        .set_vwgt(vertex_weights)
        .set_adjwgt(edge_weights)
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    let r_cut = metis_clone::Graph::new(1, nparts, xadj, adjacency)
        .unwrap()
        .set_vertex_weights(vertex_weights)
        .set_edge_weights(edge_weights)
        .seed(seed)
        .part_kway(&mut r_part)
        .unwrap();

    for i in 0..n {
        assert!(r_part[i] >= 0 && r_part[i] < nparts);
    }

    assert_eq!(c_cut, r_cut, "Edge cuts differ (weighted, seed={}, nparts={})", seed, nparts);
    assert_eq!(c_part, r_part, "Partitions differ (weighted, seed={}, nparts={})", seed, nparts);
}

/// Helper for mesh part_dual.
fn compare_part_dual(element_offsets: &[i32], element_indices: &[i32], nn: i32, nparts: i32, seed: i32) {
    let ne = element_offsets.len() - 1;
    let mut c_epart = vec![0i32; ne];
    let mut c_npart = vec![0i32; nn as usize];
    let mut r_epart = vec![0i32; ne];
    let mut r_npart = vec![0i32; nn as usize];

    let c_cut = metis::Mesh::new(nparts, element_offsets, element_indices)
        .unwrap()
        .set_option(metis::option::Seed(seed))
        .part_dual(&mut c_epart, &mut c_npart)
        .unwrap();

    let r_cut = metis_clone::Mesh::new(nparts, element_offsets, element_indices)
        .unwrap()
        .seed(seed)
        .part_dual(&mut r_epart, &mut r_npart)
        .unwrap();

    // Valid partitions
    for i in 0..ne {
        assert!(r_epart[i] >= 0 && r_epart[i] < nparts);
    }
    for i in 0..nn as usize {
        assert!(r_npart[i] >= 0 && r_npart[i] < nparts);
    }

    assert_eq!(c_cut, r_cut, "Mesh edge cuts differ (seed={}, nparts={})", seed, nparts);
    assert_eq!(c_epart, r_epart, "Element partitions differ (seed={}, nparts={})", seed, nparts);
    assert_eq!(c_npart, r_npart, "Node partitions differ (seed={}, nparts={})", seed, nparts);
}

// ========= Validity tests (check our output is at least a valid partition) =========

#[test]
fn test_trivial_2v_valid() {
    let (xadj, adjacency) = trivial_2v();
    let n = xadj.len() - 1;
    let mut part = vec![0i32; n];
    let cut = metis_clone::Graph::new(1, 2, &xadj, &adjacency)
        .unwrap()
        .seed(42)
        .part_kway(&mut part)
        .unwrap();
    // Should have 2 different partitions
    assert!(part.iter().any(|&p| p == 0));
    assert!(part.iter().any(|&p| p == 1));
    assert!(cut >= 0);
}

#[test]
fn test_path_5v_valid() {
    let (xadj, adjacency) = path_5v();
    let n = xadj.len() - 1;
    let mut part = vec![0i32; n];
    let cut = metis_clone::Graph::new(1, 2, &xadj, &adjacency)
        .unwrap()
        .seed(42)
        .part_kway(&mut part)
        .unwrap();
    for &p in &part {
        assert!(p >= 0 && p < 2);
    }
    assert!(cut > 0);
}

#[test]
fn test_grid_3x5_valid() {
    let (xadj, adjacency) = grid_3x5();
    let n = xadj.len() - 1;
    for &nparts in &[2, 3, 4] {
        let mut part = vec![0i32; n];
        let cut = metis_clone::Graph::new(1, nparts, &xadj, &adjacency)
            .unwrap()
            .seed(42)
            .part_kway(&mut part)
            .unwrap();
        for &p in &part {
            assert!(p >= 0 && p < nparts, "partition {} out of range for nparts={}", p, nparts);
        }
        assert!(cut >= 0);
    }
}

#[test]
fn test_nparts_1_returns_zeros() {
    let (xadj, adjacency) = grid_3x5();
    let n = xadj.len() - 1;
    let mut part = vec![-1i32; n];
    let cut = metis_clone::Graph::new(1, 1, &xadj, &adjacency)
        .unwrap()
        .part_kway(&mut part)
        .unwrap();
    assert_eq!(cut, 0);
    assert!(part.iter().all(|&p| p == 0));
}

#[test]
fn test_mesh_nparts_1_returns_zeros() {
    let (element_offsets, element_indices, nn) = tri_mesh();
    let ne = element_offsets.len() - 1;
    let mut epart = vec![-1i32; ne];
    let mut npart = vec![-1i32; nn as usize];
    let cut = metis_clone::Mesh::new(1, &element_offsets, &element_indices)
        .unwrap()
        .part_dual(&mut epart, &mut npart)
        .unwrap();
    assert_eq!(cut, 0);
    assert!(epart.iter().all(|&p| p == 0));
    assert!(npart.iter().all(|&p| p == 0));
}

// ========= Diagnostic tests =========

/// Test on the exact 10v subgraph that appears in the nparts=3 pipeline.
/// This isolates whether our bisection code itself works correctly.
#[test]
fn compare_10v_subgraph() {
    // 10v subgraph (original adjacency order, no reordering)
    let xadj = vec![0, 2, 5, 8, 9, 12, 16, 19, 21, 24, 26];
    let adjacency = vec![
        1, 4,           // v0
        0, 2, 5,        // v1
        1, 3, 6,        // v2
        2,              // v3
        0, 5, 7,        // v4
        1, 4, 6, 8,     // v5
        2, 5, 9,        // v6
        4, 8,           // v7
        5, 7, 9,        // v8
        6, 8,           // v9
    ];
    for &seed in &[42, 0, -1, 1128597453] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
    }
}

// ========= Comparison tests =========

#[test]
fn compare_trivial_2v() {
    let (xadj, adjacency) = trivial_2v();
    for &seed in &[42, 0, 12345] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
    }
}

#[test]
fn compare_path_5v() {
    let (xadj, adjacency) = path_5v();
    for &seed in &[42, 0, 12345] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
    }
}

#[test]
fn compare_grid_3x5_2parts() {
    let (xadj, adjacency) = grid_3x5();
    for &seed in &[42, 0, 12345] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
    }
}

#[test]
fn compare_grid_3x5_3parts() {
    let (xadj, adjacency) = grid_3x5();
    for &seed in &[42, 0, 12345] {
        compare_part_kway(&xadj, &adjacency, 3, seed);
    }
}

#[test]
fn compare_grid_3x5_4parts() {
    let (xadj, adjacency) = grid_3x5();
    for &seed in &[42, 0] {
        compare_part_kway(&xadj, &adjacency, 4, seed);
    }
}

#[test]
fn compare_irregular_6v() {
    let (xadj, adjacency) = irregular_6v();
    for &seed in &[42, 0, 12345] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
        compare_part_kway(&xadj, &adjacency, 3, seed);
    }
}

#[test]
fn compare_grid_10x10_2parts() {
    let (xadj, adjacency) = grid_10x10();
    for &seed in &[42, 0] {
        compare_part_kway(&xadj, &adjacency, 2, seed);
    }
}

#[test]
fn compare_grid_10x10_4parts() {
    let (xadj, adjacency) = grid_10x10();
    for &seed in &[42, 0] {
        compare_part_kway(&xadj, &adjacency, 4, seed);
    }
}

#[test]
fn compare_grid_10x10_8parts() {
    let (xadj, adjacency) = grid_10x10();
    for &seed in &[42] {
        compare_part_kway(&xadj, &adjacency, 8, seed);
    }
}

#[test]
fn compare_tri_mesh_dual() {
    let (element_offsets, element_indices, nn) = tri_mesh();
    for &seed in &[42, 0] {
        compare_part_dual(&element_offsets, &element_indices, nn, 2, seed);
    }
}

#[test]
fn compare_quad_mesh_dual() {
    let (element_offsets, element_indices, nn) = quad_mesh();
    for &seed in &[42, 0] {
        compare_part_dual(&element_offsets, &element_indices, nn, 2, seed);
    }
}
