mod fixtures;

use fixtures::*;

// ========= Test helpers =========

/// Mirrors the production partition_mesh calling pattern exactly.
/// Calls part_dual with ObjType::Cut on BOTH C metis and graph_partitioner,
/// asserts bit-identical results (edge cut, element partitions, node partitions).
fn compare_mesh_partition(
    element_offsets: &[i32],
    element_indices: &[i32],
    nparts: i32,
    seed: i32,
) {
    let ne = element_offsets.len() - 1;
    let nn = element_indices.iter().copied().max().unwrap_or(-1) + 1;

    let mut c_epart = vec![0i32; ne];
    let mut c_npart = vec![0i32; nn as usize];
    let mut r_epart = vec![0i32; ne];
    let mut r_npart = vec![0i32; nn as usize];

    let c_cut = metis::Mesh::new(nparts, element_offsets, element_indices)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Seed(seed))
        .part_dual(&mut c_epart, &mut c_npart)
        .unwrap();

    let r_cut = graph_partitioner::Mesh::new(nparts, element_offsets, element_indices)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .seed(seed)
        .part_dual(&mut r_epart, &mut r_npart)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "edge cuts differ (nparts={}, seed={}): C={}, Rust={}",
        nparts, seed, c_cut, r_cut
    );
    assert_eq!(
        c_epart, r_epart,
        "element partitions differ (nparts={}, seed={})",
        nparts, seed
    );
    assert_eq!(
        c_npart, r_npart,
        "node partitions differ (nparts={}, seed={})",
        nparts, seed
    );
}

/// Mirrors the production partition_mesh function that takes raw u32 triangle indices,
/// remaps to compact [0, N), and calls part_dual with ObjType::Cut.
/// Asserts bit-identical results between C metis and graph_partitioner.
fn compare_partition_mesh_remapped(indices: &[u32], num_parts: i32, seed: i32) {
    let num_tris = indices.len() / 3;
    assert_eq!(indices.len() % 3, 0, "indices length must be a multiple of 3");

    // Remap to compact [0, N) — mirrors the production code exactly
    let max_v = indices.iter().copied().max().unwrap_or(0) as usize;
    let mut remap: Vec<i32> = vec![-1; max_v + 1];
    let mut used_v: i32 = 0;
    let mut remapped: Vec<i32> = Vec::with_capacity(indices.len());
    for &i in indices {
        let slot = &mut remap[i as usize];
        if *slot < 0 {
            *slot = used_v;
            used_v += 1;
        }
        remapped.push(*slot);
    }

    let eptr: Vec<i32> = (0..=num_tris).map(|i| (i * 3) as i32).collect();

    let ne = num_tris;
    let nn = used_v;

    let mut c_epart = vec![0i32; ne];
    let mut c_npart = vec![0i32; nn as usize];
    let mut r_epart = vec![0i32; ne];
    let mut r_npart = vec![0i32; nn as usize];

    let c_cut = metis::Mesh::new(num_parts, &eptr, &remapped)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Seed(seed))
        .part_dual(&mut c_epart, &mut c_npart)
        .unwrap();

    let r_cut = graph_partitioner::Mesh::new(num_parts, &eptr, &remapped)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .seed(seed)
        .part_dual(&mut r_epart, &mut r_npart)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "remapped mesh: edge cuts differ (nparts={}, seed={})",
        num_parts, seed
    );
    assert_eq!(
        c_epart, r_epart,
        "remapped mesh: element partitions differ (nparts={}, seed={})",
        num_parts, seed
    );
    assert_eq!(
        c_npart, r_npart,
        "remapped mesh: node partitions differ (nparts={}, seed={})",
        num_parts, seed
    );
}

/// Mirrors the production group_component calling pattern exactly.
/// Weighted CSR graph, part_kway with ObjType::Cut + edge weights on BOTH
/// C metis and graph_partitioner, asserts bit-identical results.
fn compare_weighted_kway(
    xadj: &[i32],
    adjncy: &[i32],
    adjwgt: &[i32],
    nparts: i32,
    seed: i32,
) {
    let n = xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    let c_cut = metis::Graph::new(1, nparts, xadj, adjncy)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_adjwgt(adjwgt)
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    let r_cut = graph_partitioner::Graph::new(1, nparts, xadj, adjncy)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .set_edge_weights(adjwgt)
        .seed(seed)
        .part_kway(&mut r_part)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "weighted kway: edge cuts differ (nparts={}, seed={}): C={}, Rust={}",
        nparts, seed, c_cut, r_cut
    );
    assert_eq!(
        c_part, r_part,
        "weighted kway: partitions differ (nparts={}, seed={})",
        nparts, seed
    );
}

// ========= Grid/strip mesh tests =========

#[test]
fn mesh_strip_20_2parts() {
    let (eoff, eind) = tri_strip_mesh(20);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 2, seed);
    }
}

#[test]
fn mesh_strip_20_4parts() {
    let (eoff, eind) = tri_strip_mesh(20);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_grid_5x5_4parts() {
    let (eoff, eind) = tri_grid_mesh(5, 5);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_grid_5x5_8parts() {
    let (eoff, eind) = tri_grid_mesh(5, 5);
    compare_mesh_partition(&eoff, &eind, 8, 42);
}

#[test]
fn mesh_grid_10x10_16parts() {
    let (eoff, eind) = tri_grid_mesh(10, 10);
    compare_mesh_partition(&eoff, &eind, 16, 42);
}

#[test]
fn mesh_sparse_ids_4parts() {
    let sparse_indices = sparse_tri_mesh(30);
    for &seed in &[42, 0] {
        compare_partition_mesh_remapped(&sparse_indices, 4, seed);
    }
}

#[test]
fn mesh_single_tri() {
    let eoff = vec![0, 3];
    let eind = vec![0, 1, 2];
    // 1 element, 2 parts => nparts > ne, should still work
    compare_mesh_partition(&eoff, &eind, 2, 42);
}

#[test]
fn mesh_two_tris() {
    let eoff = vec![0, 3, 6];
    let eind = vec![0, 1, 2, 1, 3, 2];
    compare_mesh_partition(&eoff, &eind, 2, 42);
}

// ========= Procedural surface mesh tests =========

#[test]
fn mesh_sphere_small_4parts() {
    let (eoff, eind) = uv_sphere_mesh(8, 16);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_sphere_large_8parts() {
    let (eoff, eind) = uv_sphere_mesh(16, 32);
    compare_mesh_partition(&eoff, &eind, 8, 42);
}

#[test]
fn mesh_sphere_large_16parts() {
    let (eoff, eind) = uv_sphere_mesh(16, 32);
    compare_mesh_partition(&eoff, &eind, 16, 42);
}

#[test]
fn mesh_torus_4parts() {
    let (eoff, eind) = torus_mesh(16, 8);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_torus_8parts() {
    let (eoff, eind) = torus_mesh(24, 12);
    compare_mesh_partition(&eoff, &eind, 8, 42);
}

#[test]
fn mesh_disjoint_spheres_2parts() {
    let (eoff, eind) = two_disjoint_spheres_mesh(8);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 2, seed);
    }
}

#[test]
fn mesh_disjoint_spheres_4parts() {
    let (eoff, eind) = two_disjoint_spheres_mesh(8);
    compare_mesh_partition(&eoff, &eind, 4, 42);
}

#[test]
fn mesh_mobius_4parts() {
    let (eoff, eind) = mobius_strip_mesh(16, 2);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_mobius_8parts() {
    let (eoff, eind) = mobius_strip_mesh(32, 3);
    compare_mesh_partition(&eoff, &eind, 8, 42);
}

#[test]
fn mesh_cylinder_4parts() {
    let (eoff, eind) = cylinder_mesh(8, 16);
    for &seed in &[42, 0] {
        compare_mesh_partition(&eoff, &eind, 4, seed);
    }
}

#[test]
fn mesh_cylinder_8parts() {
    let (eoff, eind) = cylinder_mesh(16, 24);
    compare_mesh_partition(&eoff, &eind, 8, 42);
}

// ========= Weighted graph tests =========

#[test]
fn weighted_grid_3x5_2parts() {
    let (xadj, adjncy, adjwgt) = weighted_grid_graph(3, 5);
    for &seed in &[42, 0] {
        compare_weighted_kway(&xadj, &adjncy, &adjwgt, 2, seed);
    }
}

#[test]
fn weighted_grid_3x5_3parts() {
    let (xadj, adjncy, adjwgt) = weighted_grid_graph(3, 5);
    for &seed in &[42, 0] {
        compare_weighted_kway(&xadj, &adjncy, &adjwgt, 3, seed);
    }
}

#[test]
fn weighted_grid_5x5_4parts() {
    let (xadj, adjncy, adjwgt) = weighted_grid_graph(5, 5);
    compare_weighted_kway(&xadj, &adjncy, &adjwgt, 4, 42);
}

#[test]
fn weighted_grid_10x10_8parts() {
    let (xadj, adjncy, adjwgt) = weighted_grid_graph(10, 10);
    compare_weighted_kway(&xadj, &adjncy, &adjwgt, 8, 42);
}

#[test]
fn weighted_irregular_adjwgt() {
    let (xadj, adjncy) = irregular_6v();
    // Custom edge weights: varying weights for each edge
    // irregular_6v adjacency: [1,3, 0,2,3, 1,4, 0,1,4,5, 2,3,5, 3,4]
    let adjwgt = vec![2, 5, 2, 3, 1, 3, 7, 5, 1, 7, 4, 3, 7, 4, 7, 3];
    assert_eq!(adjncy.len(), adjwgt.len());
    for &seed in &[42, 0] {
        compare_weighted_kway(&xadj, &adjncy, &adjwgt, 2, seed);
        compare_weighted_kway(&xadj, &adjncy, &adjwgt, 3, seed);
    }
}

#[test]
fn weighted_single_edge() {
    let xadj = vec![0, 1, 2];
    let adjncy = vec![1, 0];
    let adjwgt = vec![10, 10];
    compare_weighted_kway(&xadj, &adjncy, &adjwgt, 2, 42);
}

// ========= Diagnostic: unweighted graph kway tests =========

/// Compare unweighted graph kway partition (no edge weights, no vertex weights).
fn compare_unweighted_kway(
    xadj: &[i32],
    adjncy: &[i32],
    nparts: i32,
    seed: i32,
) {
    let n = xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    let c_cut = metis::Graph::new(1, nparts, xadj, adjncy)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    let r_cut = graph_partitioner::Graph::new(1, nparts, xadj, adjncy)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .seed(seed)
        .part_kway(&mut r_part)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "unweighted kway: edge cuts differ (nparts={}, seed={}): C={}, Rust={}",
        nparts, seed, c_cut, r_cut
    );
    assert_eq!(
        c_part, r_part,
        "unweighted kway: partitions differ (nparts={}, seed={})",
        nparts, seed
    );
}

#[test]
fn unweighted_grid_10x10_4parts() {
    let (xadj, adjncy) = grid_10x10();
    compare_unweighted_kway(&xadj, &adjncy, 4, 42);
}

#[test]
fn unweighted_grid_10x10_8parts() {
    let (xadj, adjncy) = grid_10x10();
    compare_unweighted_kway(&xadj, &adjncy, 8, 42);
}

#[test]
fn unweighted_grid_10x10_16parts() {
    let (xadj, adjncy) = grid_10x10();
    compare_unweighted_kway(&xadj, &adjncy, 16, 42);
}

#[test]
fn unweighted_grid_3x5_2parts() {
    let (xadj, adjncy) = grid_3x5();
    compare_unweighted_kway(&xadj, &adjncy, 2, 42);
}

#[test]
fn unweighted_grid_3x5_4parts() {
    let (xadj, adjncy) = grid_3x5();
    compare_unweighted_kway(&xadj, &adjncy, 4, 42);
}

#[test]
fn unweighted_grid_15x15_16parts() {
    let (xadj, adjncy) = build_grid(15, 15); // 225 vertices
    compare_unweighted_kway(&xadj, &adjncy, 16, 42);
}

#[test]
fn unweighted_grid_20x20_4parts() {
    let (xadj, adjncy) = build_grid(20, 20); // 400 vertices
    compare_unweighted_kway(&xadj, &adjncy, 4, 42);
}

#[test]
fn unweighted_grid_20x20_16parts() {
    let (xadj, adjncy) = build_grid(20, 20); // 400 vertices
    compare_unweighted_kway(&xadj, &adjncy, 16, 42);
}

// ========= Diagnostic: dual graph as kway input =========

/// Build the dual graph from a mesh using Rust code, then partition it
/// via Graph::part_kway using BOTH C METIS and Rust.
/// If this passes but mesh_sphere_small_4parts fails, the dual graph
/// constructed by C METIS differs from ours.
#[test]
fn diagnostic_sphere_dual_graph_kway_4parts() {
    let (eoff, eind) = uv_sphere_mesh(8, 16);
    let ne = (eoff.len() - 1) as i32;
    let nn = eind.iter().copied().max().unwrap_or(-1) + 1;

    // Build dual graph using Rust's create_graph_dual
    let (xadj, adjacency) = graph_partitioner::create_graph_dual(ne, nn, &eoff, &eind, 1);

    for &seed in &[42, 0] {
        compare_unweighted_kway(&xadj, &adjacency, 4, seed);
    }
}

#[test]
fn diagnostic_grid10x10_dual_graph_kway_16parts() {
    let (eoff, eind) = tri_grid_mesh(10, 10);
    let ne = (eoff.len() - 1) as i32;
    let nn = eind.iter().copied().max().unwrap_or(-1) + 1;

    let (xadj, adjacency) = graph_partitioner::create_graph_dual(ne, nn, &eoff, &eind, 1);

    compare_unweighted_kway(&xadj, &adjacency, 16, 42);
}

// ========= Production use case: dual graph with ncommon=2 + part_kway + force_contiguous =========

/// Mirrors the production use case: build dual graph with ncommon=2 (shared-edge adjacency),
/// then partition with ObjType::Cut + force_contiguous(true).
/// Compares C METIS mesh_to_dual + Graph::part_kway vs Rust create_graph_dual + Graph::part_kway.
fn compare_dual_ncommon2_kway_contiguous(
    element_offsets: &[i32],
    element_indices: &[i32],
    nparts: i32,
    seed: i32,
) {
    let ne = (element_offsets.len() - 1) as i32;
    let nn = element_indices.iter().copied().max().unwrap_or(-1) + 1;

    // Build dual graph with ncommon=2 using C METIS
    let c_dual = metis::mesh_to_dual(element_offsets, element_indices, 2).unwrap();
    let c_xadj = c_dual.xadj();
    let c_adjncy = c_dual.adjncy();

    // Build dual graph with ncommon=2 using Rust
    let (r_xadj, r_adjncy) = graph_partitioner::create_graph_dual(ne, nn, element_offsets, element_indices, 2);

    // Verify dual graphs match
    assert_eq!(c_xadj, &r_xadj[..], "Dual graph xadj differs (ncommon=2, ne={}, seed={})", ne, seed);
    assert_eq!(c_adjncy, &r_adjncy[..], "Dual graph adjncy differs (ncommon=2, ne={}, seed={})", ne, seed);

    let n = c_xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    // Partition with C METIS: ObjType::Cut + Contig(true)
    let c_cut = metis::Graph::new(1, nparts, c_xadj, c_adjncy)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Contig(true))
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    // Partition with Rust: ObjType::Cut + force_contiguous(true)
    let r_cut = graph_partitioner::Graph::new(1, nparts, &r_xadj, &r_adjncy)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .force_contiguous(true)
        .seed(seed)
        .part_kway(&mut r_part)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "dual ncommon=2 kway contiguous: edge cuts differ (nparts={}, seed={}): C={}, Rust={}",
        nparts, seed, c_cut, r_cut
    );
    assert_eq!(
        c_part, r_part,
        "dual ncommon=2 kway contiguous: partitions differ (nparts={}, seed={})",
        nparts, seed
    );
}

/// Mirrors the production use case with part_recursive instead of part_kway.
fn compare_dual_ncommon2_recursive(
    element_offsets: &[i32],
    element_indices: &[i32],
    nparts: i32,
    seed: i32,
) {
    let ne = (element_offsets.len() - 1) as i32;
    let nn = element_indices.iter().copied().max().unwrap_or(-1) + 1;

    // Build dual graph with ncommon=2 using C METIS
    let c_dual = metis::mesh_to_dual(element_offsets, element_indices, 2).unwrap();
    let c_xadj = c_dual.xadj();
    let c_adjncy = c_dual.adjncy();

    // Build dual graph with ncommon=2 using Rust
    let (r_xadj, r_adjncy) = graph_partitioner::create_graph_dual(ne, nn, element_offsets, element_indices, 2);

    let n = c_xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    // Partition with C METIS part_recursive
    let c_cut = metis::Graph::new(1, nparts, c_xadj, c_adjncy)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Seed(seed))
        .part_recursive(&mut c_part)
        .unwrap();

    // Partition with Rust part_recursive
    let r_cut = graph_partitioner::Graph::new(1, nparts, &r_xadj, &r_adjncy)
        .unwrap()
        .obj_type(graph_partitioner::option::ObjType::Cut)
        .seed(seed)
        .part_recursive(&mut r_part)
        .unwrap();

    assert_eq!(
        c_cut, r_cut,
        "dual ncommon=2 recursive: edge cuts differ (nparts={}, seed={}): C={}, Rust={}",
        nparts, seed, c_cut, r_cut
    );
    assert_eq!(
        c_part, r_part,
        "dual ncommon=2 recursive: partitions differ (nparts={}, seed={})",
        nparts, seed
    );
}

#[test]
fn production_dual_ncommon2_grid5x5_4parts() {
    let (eoff, eind) = tri_grid_mesh(5, 5);
    for &seed in &[42, 0] {
        compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 4, seed);
    }
}

#[test]
fn production_dual_ncommon2_grid10x10_8parts() {
    let (eoff, eind) = tri_grid_mesh(10, 10);
    compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 8, 42);
}

#[test]
fn production_dual_ncommon2_grid10x10_16parts() {
    let (eoff, eind) = tri_grid_mesh(10, 10);
    compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 16, 42);
}

#[test]
fn production_dual_ncommon2_sphere_4parts() {
    let (eoff, eind) = uv_sphere_mesh(8, 16);
    for &seed in &[42, 0] {
        compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 4, seed);
    }
}

#[test]
fn production_dual_ncommon2_sphere_8parts() {
    let (eoff, eind) = uv_sphere_mesh(16, 32);
    compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 8, 42);
}

#[test]
fn production_dual_ncommon2_torus_4parts() {
    let (eoff, eind) = torus_mesh(16, 8);
    for &seed in &[42, 0] {
        compare_dual_ncommon2_kway_contiguous(&eoff, &eind, 4, seed);
    }
}

// ========= part_recursive on dual graphs =========

#[test]
fn production_recursive_dual_ncommon2_grid5x5_4parts() {
    let (eoff, eind) = tri_grid_mesh(5, 5);
    for &seed in &[42, 0] {
        compare_dual_ncommon2_recursive(&eoff, &eind, 4, seed);
    }
}

#[test]
fn production_recursive_dual_ncommon2_grid10x10_8parts() {
    let (eoff, eind) = tri_grid_mesh(10, 10);
    compare_dual_ncommon2_recursive(&eoff, &eind, 8, 42);
}

#[test]
fn production_recursive_dual_ncommon2_sphere_4parts() {
    let (eoff, eind) = uv_sphere_mesh(8, 16);
    for &seed in &[42, 0] {
        compare_dual_ncommon2_recursive(&eoff, &eind, 4, seed);
    }
}

// Sweep different sizes and nparts to find minimal failing case
#[test]
fn diagnostic_dual_sweep() {
    let configs = [
        (3, 3, 2), (3, 3, 4), (4, 4, 2), (4, 4, 4),
        (5, 5, 2), (5, 5, 4), (5, 5, 8),
        (6, 6, 2), (6, 6, 4), (6, 6, 8),
        (7, 7, 2), (7, 7, 4), (7, 7, 8),
        (8, 8, 2), (8, 8, 4), (8, 8, 8),
    ];
    for &(rows, cols, nparts) in &configs {
        let (eoff, eind) = tri_grid_mesh(rows, cols);
        let ne = (eoff.len() - 1) as i32;
        let nn = eind.iter().copied().max().unwrap_or(-1) + 1;
        let (xadj, adjacency) = graph_partitioner::create_graph_dual(ne, nn, &eoff, &eind, 1);

        let n = xadj.len() - 1;
        let mut c_part = vec![0i32; n];
        let mut r_part = vec![0i32; n];

        let c_cut = metis::Graph::new(1, nparts, &xadj, &adjacency)
            .unwrap()
            .set_option(metis::option::ObjType::Cut)
            .set_option(metis::option::Seed(42))
            .part_kway(&mut c_part)
            .unwrap();

        let r_cut = graph_partitioner::Graph::new(1, nparts, &xadj, &adjacency)
            .unwrap()
            .obj_type(graph_partitioner::option::ObjType::Cut)
            .seed(42)
            .part_kway(&mut r_part)
            .unwrap();

        let ne_count = ne;
        let max_deg = (0..n).map(|i| xadj[i+1] - xadj[i]).max().unwrap_or(0);
        if c_cut != r_cut || c_part != r_part {
            eprintln!("FAIL: {}x{} tri mesh, {} dual vertices, max_deg={}, nparts={}: C_cut={} Rust_cut={}",
                rows, cols, ne_count, max_deg, nparts, c_cut, r_cut);
        } else {
            eprintln!("PASS: {}x{} tri mesh, {} dual vertices, max_deg={}, nparts={}: cut={}",
                rows, cols, ne_count, max_deg, nparts, c_cut);
        }
    }
}

