# Plan: Integration Tests & Code Coverage

## Context

The existing test suite has 61 tests: 31 PQueue unit tests, 17 comparison tests, and 13 debug tests. The comparison tests cover basic graph/mesh shapes (2v–100v, 2–8 parts) but don't test the calling patterns from real production usage. There are two specific call patterns from the consuming codebase that need dedicated integration tests, plus we need code coverage tooling since the library is very branchy.

---

## Part A: Integration Tests Matching Production Call Patterns

### Call Pattern 1: `partition_mesh` (triangle mesh via `part_dual` with `ObjType::Cut`)

Production code does:
1. Takes raw triangle indices (arbitrary vertex IDs, possibly sparse)
2. Remaps vertices to compact `[0, N)` space
3. Builds uniform-stride element offsets (every 3)
4. Calls `Mesh::new(num_parts, &eptr, &remapped).set_option(ObjType::Cut).part_dual(&mut epart, &mut npart)`

This pattern is **not currently tested** — existing mesh tests don't set `ObjType::Cut` and only use tiny meshes (4 triangles, 2 quads).

### Call Pattern 2: `group_component` (weighted graph via `part_kway` with `ObjType::Cut` + edge weights)

Production code does:
1. Builds a CSR graph from an adjacency list with edge weights
2. Calls `Graph::new(1, target_groups, &xadj, &adjncy).set_option(ObjType::Cut).set_adjwgt(&adjwgt).part_kway(&mut parts)`

This pattern is **partially tested** — `compare_part_kway_weighted` exists as a helper but is `#[allow(dead_code)]` and never called.

---

### New Test File: `tests/integration.rs`

#### Fixtures to add (in `tests/fixtures.rs`)

**Grid/strip meshes:**

1. **`tri_strip_mesh(n: usize)`** — Generate a triangle strip with `n` triangles. Simple, deterministic, scalable mesh.
   ```
   0---1---2---3---4
   |\ | /|\ | /|
   | \|/ | \|/ |
   5---6---7---8---9
   ```
   Returns `(element_offsets, element_indices)` with compact vertex IDs.

2. **`tri_grid_mesh(rows: usize, cols: usize)`** — Generate a grid triangulated into `2 * rows * cols` triangles. Matches the typical game/rendering mesh topology. Returns `(element_offsets, element_indices)`.

3. **`sparse_tri_mesh(n_tris: usize)`** — Generate a triangle mesh with **non-contiguous vertex IDs** (e.g., IDs `[0, 5, 10, 15, ...]`). Tests the remapping step from the production `partition_mesh` function.

**Procedural surface meshes:**

4. **`uv_sphere_mesh(lat_segments: usize, lon_segments: usize)`** — Triangulated UV sphere. `2 * lon_segments` triangles at the poles (fans), `2 * (lat_segments - 2) * lon_segments` quad-split triangles in the middle bands. Total: `2 * lon_segments * (lat_segments - 1)` triangles. Closed manifold, no boundary, vertex valence ~6 except poles (valence = `lon_segments`). Calling with e.g. `(16, 32)` gives 960 triangles — a realistic mesh size.

5. **`torus_mesh(major_segments: usize, minor_segments: usize)`** — Triangulated torus (genus-1 surface). Each quad patch splits into 2 triangles, giving `2 * major_segments * minor_segments` triangles total. No boundary, no poles, uniform vertex valence 6. Exercises partitioning on a surface with non-trivial topology where any cut must be a closed loop. E.g. `(16, 8)` = 256 triangles.

6. **`two_disjoint_spheres_mesh(segments: usize)`** — Two separate UV spheres with non-overlapping vertex IDs (second sphere's IDs offset by first sphere's vertex count). Tests partitioning of **disconnected components** — METIS should ideally assign each component to separate partitions. Important edge case: the dual graph has two connected components.

7. **`mobius_strip_mesh(length_segments: usize, width_segments: usize)`** — Triangulated Mobius strip (non-orientable surface with a single boundary loop). Built by taking a rectangular strip of `length_segments x width_segments` quads, triangulating each, then connecting the last column to the first with a half-twist (flipping the width index). Total: `2 * length_segments * width_segments` triangles. Tests partitioning on non-orientable topology with boundary. E.g. `(16, 2)` = 64 triangles.

8. **`cylinder_mesh(height_segments: usize, radial_segments: usize)`** — Open cylinder (tube). Like a torus but without closing the second axis — has two circular boundary loops (top and bottom edges). `2 * height_segments * radial_segments` triangles. Provides a contrast to the torus: same local connectivity but with boundary edges that create natural cut locations.

**Weighted graph fixture:**

9. **`weighted_grid_graph(rows: usize, cols: usize)`** — Returns `(xadj, adjncy, adjwgt)` where edge weights vary (e.g., horizontal edges weight 1, vertical edges weight 2). Matches the weighted adjacency pattern from `group_component`.

#### Test helpers

Every test calls **both** C METIS and metis_clone with identical inputs and asserts **bit-exact** output — same edge cut value, same partition arrays element-for-element. No validity-only tests.

```rust
/// Mirrors the production partition_mesh function exactly.
/// Remaps sparse vertex indices, builds eptr, calls part_dual with ObjType::Cut
/// on BOTH C metis and metis_clone, asserts bit-identical results.
fn partition_mesh_compare(indices: &[i32], num_parts: i32, seed: i32) {
    let num_tris = indices.len() / 3;
    let mut eptr: Vec<i32> = (0..=num_tris).map(|i| (i * 3) as i32).collect();
    let nn = indices.iter().copied().max().unwrap_or(-1) + 1;

    let mut c_epart = vec![0i32; num_tris];
    let mut c_npart = vec![0i32; nn as usize];
    let mut r_epart = vec![0i32; num_tris];
    let mut r_npart = vec![0i32; nn as usize];

    let c_cut = metis::Mesh::new(num_parts, &eptr, indices)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_option(metis::option::Seed(seed))
        .part_dual(&mut c_epart, &mut c_npart)
        .unwrap();

    let r_cut = metis_clone::Mesh::new(num_parts, &eptr, indices)
        .unwrap()
        .set_option(metis_clone::option::ObjType::Cut)
        .set_option(metis_clone::option::Seed(seed))
        .part_dual(&mut r_epart, &mut r_npart)
        .unwrap();

    assert_eq!(c_cut, r_cut, "edge cuts differ");
    assert_eq!(c_epart, r_epart, "element partitions differ");
    assert_eq!(c_npart, r_npart, "node partitions differ");
}

/// Mirrors the production group_component function exactly.
/// Weighted CSR graph, part_kway with ObjType::Cut + edge weights on BOTH
/// C metis and metis_clone, asserts bit-identical results.
fn group_component_compare(
    xadj: &[i32], adjncy: &[i32], adjwgt: &[i32],
    target_groups: i32, seed: i32,
) {
    let n = xadj.len() - 1;
    let mut c_part = vec![0i32; n];
    let mut r_part = vec![0i32; n];

    let c_cut = metis::Graph::new(1, target_groups, xadj, adjncy)
        .unwrap()
        .set_option(metis::option::ObjType::Cut)
        .set_adjwgt(adjwgt)
        .set_option(metis::option::Seed(seed))
        .part_kway(&mut c_part)
        .unwrap();

    let r_cut = metis_clone::Graph::new(1, target_groups, xadj, adjncy)
        .unwrap()
        .set_option(metis_clone::option::ObjType::Cut)
        .set_edge_weights(adjwgt)
        .set_option(metis_clone::option::Seed(seed))
        .part_kway(&mut r_part)
        .unwrap();

    assert_eq!(c_cut, r_cut, "edge cuts differ");
    assert_eq!(c_part, r_part, "partitions differ");
}
```

#### Test cases

**Mesh partitioning (partition_mesh pattern):**

Grid/strip meshes:

| Test | Mesh | num_parts | Seeds | What it exercises |
|------|------|-----------|-------|-------------------|
| `mesh_strip_20_2parts` | tri_strip(20) | 2 | 42, 0 | Basic mesh partition with ObjType::Cut |
| `mesh_strip_20_4parts` | tri_strip(20) | 4 | 42, 0 | Multi-way mesh partition |
| `mesh_grid_5x5_4parts` | tri_grid(5,5) = 50 tris | 4 | 42, 0 | Larger mesh, 4-way |
| `mesh_grid_5x5_8parts` | tri_grid(5,5) | 8 | 42 | Many parts relative to elements |
| `mesh_grid_10x10_16parts` | tri_grid(10,10) = 200 tris | 16 | 42 | Large mesh, high part count (exercises deeper coarsening) |
| `mesh_sparse_ids` | sparse_tri_mesh(30) | 4 | 42, 0 | Non-contiguous vertex IDs (tests the remapping path) |
| `mesh_single_tri` | 1 triangle | 2 | 42 | Edge case: fewer elements than parts |
| `mesh_two_tris` | 2 triangles | 2 | 42 | Minimal non-trivial case |

Procedural surface meshes:

| Test | Mesh | num_parts | Seeds | What it exercises |
|------|------|-----------|-------|-------------------|
| `mesh_sphere_small_4parts` | uv_sphere(8, 16) = 240 tris | 4 | 42, 0 | Closed manifold, pole singularities (high-valence vertices) |
| `mesh_sphere_large_8parts` | uv_sphere(16, 32) = 960 tris | 8 | 42 | Realistic mesh size, deep coarsening hierarchy |
| `mesh_sphere_large_16parts` | uv_sphere(16, 32) | 16 | 42 | High part count on closed manifold |
| `mesh_torus_4parts` | torus(16, 8) = 256 tris | 4 | 42, 0 | Genus-1 surface — any partition cut must be a loop, no natural "pinch points" |
| `mesh_torus_8parts` | torus(24, 12) = 576 tris | 8 | 42 | Larger torus, more parts |
| `mesh_disjoint_spheres_2parts` | two_disjoint_spheres(8) = 480 tris | 2 | 42, 0 | Disconnected components — ideal result: one sphere per part |
| `mesh_disjoint_spheres_4parts` | two_disjoint_spheres(8) | 4 | 42 | Disconnected components with more parts than components |
| `mesh_mobius_4parts` | mobius_strip(16, 2) = 64 tris | 4 | 42, 0 | Non-orientable surface with boundary |
| `mesh_mobius_8parts` | mobius_strip(32, 3) = 192 tris | 8 | 42 | Larger non-orientable mesh |
| `mesh_cylinder_4parts` | cylinder(8, 16) = 256 tris | 4 | 42, 0 | Surface with two boundary loops (natural cut locations along height) |
| `mesh_cylinder_8parts` | cylinder(16, 24) = 768 tris | 8 | 42 | Larger cylinder |

**Weighted graph partitioning (group_component pattern):**

| Test | Graph | num_parts | Seeds | What it exercises |
|------|-------|-----------|-------|-------------------|
| `weighted_grid_3x5_2parts` | weighted_grid(3,5) | 2 | 42, 0 | Edge weights + ObjType::Cut |
| `weighted_grid_3x5_3parts` | weighted_grid(3,5) | 3 | 42, 0 | Multi-way weighted |
| `weighted_grid_5x5_4parts` | weighted_grid(5,5) | 4 | 42 | Larger weighted graph |
| `weighted_grid_10x10_8parts` | weighted_grid(10,10) | 8 | 42 | Large weighted, many parts |
| `weighted_irregular_adjwgt` | irregular_6v + custom weights | 2, 3 | 42, 0 | Non-uniform edge weights |
| `weighted_single_edge` | 2 vertices, 1 weighted edge | 2 | 42 | Minimal weighted case |

**ObjType::Cut specific tests** (verify the option actually changes behavior):

| Test | What it checks |
|------|----------------|
| `objtype_cut_vs_default_graph` | Compare partition with and without ObjType::Cut on a graph that produces different results |
| `objtype_cut_vs_default_mesh` | Same for mesh |

---

### Files to modify

| File | Change |
|------|--------|
| `tests/fixtures.rs` | Add `tri_strip_mesh`, `tri_grid_mesh`, `sparse_tri_mesh`, `uv_sphere_mesh`, `torus_mesh`, `two_disjoint_spheres_mesh`, `mobius_strip_mesh`, `cylinder_mesh`, `weighted_grid_graph` generators |
| `tests/integration.rs` | New file: all integration test cases with `partition_mesh_compare` and `group_component_compare` helpers |

---

## Part B: Code Coverage Setup

### Tool: `cargo-llvm-cov`

Use `cargo-llvm-cov` over `cargo-tarpaulin` because:
- It uses LLVM's native instrumentation (same as `rustc`), so it handles branches, generics, and inlining correctly
- Better branch coverage support (important for this "very branchy" codebase)
- Works on Windows (tarpaulin is Linux-only)
- Can output LCOV format for IDE integration and HTML reports

### Setup

**Install:**
```
cargo install cargo-llvm-cov
```

Or with rustup:
```
rustup component add llvm-tools-preview
cargo install cargo-llvm-cov
```

**Add to `Cargo.toml`:**
No changes needed — `cargo-llvm-cov` works as an external cargo subcommand.

**Add `.cargo/llvm-cov-config` or use command-line flags.**

### Usage commands

```bash
# Run tests with coverage, print summary
cargo llvm-cov

# Generate HTML report (opens in browser)
cargo llvm-cov --html
# Output goes to target/llvm-cov/html/index.html

# Generate LCOV format (for CI/IDE integration)
cargo llvm-cov --lcov --output-path lcov.info

# Branch coverage (the key feature for a branchy codebase)
cargo llvm-cov --branch

# Branch coverage with HTML report
cargo llvm-cov --branch --html

# Show only uncovered lines
cargo llvm-cov --html --show-missing-lines
```

### Configuration file: `codecov.toml`

```toml
[report]
# Files to exclude from coverage
exclude = [
    "tests/*",
    "src/rng.rs",     # Simple LCG, already verified via comparison tests
]
```

### .gitignore additions

```
lcov.info
target/llvm-cov/
```

### Coverage targets

Given the branchy nature of the library, focus coverage analysis on:

1. **`src/partition/kwayfm.rs`** — Greedy k-way FM refinement (most branches: balance checks, gain computation, vertex selection)
2. **`src/partition/fm.rs`** — 2-way FM refinement (moved/unmoved vertex logic, gain tracking)
3. **`src/partition/pmetis.rs`** — Recursive bisection (different split strategies based on graph size)
4. **`src/graph/coarsen.rs`** — Coarsening strategies (2-hop matching, different matching policies)
5. **`src/contig.rs`** — Contiguity enforcement (only triggered when `force_contiguous` option is set)
6. **`src/minconn.rs`** — Connectivity minimization (only triggered when `minimize_connectivity` option is set)
7. **`src/balance_kway.rs`** — K-way balancing (triggered when partition is imbalanced)

### Coverage-driven test additions

After getting the initial coverage report, add targeted tests to hit uncovered branches. Likely gaps:

| Area | Why it might be uncovered | Test to add |
|------|--------------------------|-------------|
| `contig.rs` | No tests set `Contig(true)` | Test with `.set_option(Contig(true))` on a graph where partitions aren't naturally contiguous |
| `minconn.rs` | No tests set `MinConn(true)` | Test with `.set_option(MinConn(true))` |
| Volume objective | No tests use `ObjType::Vol` | Test with `ObjType::Vol` (different code paths in kwayfm) |
| `IpType` variants | Only default `IpType` tested | Test `IpType::Random`, `IpType::Edge`, `IpType::Grow` |
| `CType` variants | Only default tested | Test `CType::Rm` vs `CType::Shem` |
| `compress` path | Never tested | Test `Compress(true)` on graphs with many degree-1 vertices |
| `No2Hop` path | Never tested | Test `No2Hop(true)` to disable 2-hop matching |
| Multi-constraint | All tests use `ncon=1` | Test with `ncon=2` and appropriate vertex weight arrays |
| `nparts > vertices` | Edge case | Test where nparts exceeds the number of vertices |
| `UFactor` | Never set | Test with different `UFactor` values to exercise balance paths |

These tests go in a new file `tests/coverage_options.rs` — each test calls both C metis and metis_clone to verify bit-identical behavior across all option combinations.

---

## Execution Order

1. Add fixture generators to `tests/fixtures.rs`
2. Create `tests/integration.rs` with production-pattern tests
3. Install `cargo-llvm-cov` and run initial coverage report
4. Identify uncovered branches from the HTML report
5. Create `tests/coverage_options.rs` with targeted tests for each uncovered option/path
6. Re-run coverage to verify improvement
7. Add `.gitignore` entries for coverage artifacts

## Verification

```bash
# All tests pass
cargo test

# Coverage report
cargo llvm-cov --branch --html
# Review target/llvm-cov/html/index.html
```
