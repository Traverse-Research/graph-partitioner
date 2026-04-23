# Plan: Idiomatic Rust Cleanup

## Context

The codebase is a working pure-Rust METIS clone (61 tests passing, bit-identical output to C METIS). It was ported from C and still has several C-isms: a fixed-size `[i32; 40]` options array with `-1` sentinels, raw pointer chains for the coarsening hierarchy, and `unsafe` blocks for pointer traversal. This plan replaces those patterns with idiomatic Rust while preserving identical test output.

---

## Phase 1: Options System Rewrite (Builder Pattern)

**Risk: Medium. Touches public API, internal construction, and one internal call site.**

Replace the C-style `[Idx; 40]` options array with a proper `Options` struct and traditional builder methods. Remove the `Opt` trait, `SetOption` trait, and all the sealed-trait machinery entirely.

### Current flow

1. Public API (`Graph`/`Mesh` in `src/lib.rs`): stores `options: [Idx; NOPTIONS]`, filled via `set_option<O: Opt>` which writes `options[O::INDEX] = option.value()`
2. Entry points (`src/partition/kway.rs`, `src/mesh/meshpart.rs`): pass `&graph.options` to `Control::new()`
3. `Control::new()` in `src/ctrl.rs`: reads each option by index with `-1` = "use default" sentinel via `get(idx, default)` closure
4. Internal recursive bisection (`src/partition/pmetis.rs:35-54`): manually constructs `[-1; NOPTIONS]` and sets `options[8]` and `options[7]` by magic index

### New design

**`Options` struct** (replaces `[Idx; NOPTIONS]`):

```rust
// src/option.rs
#[derive(Debug, Clone, Default)]
pub struct Options {
    pub ptype: Option<PType>,
    pub obj_type: Option<ObjType>,
    pub coarsen_type: Option<CType>,
    pub init_part_type: Option<IpType>,
    pub refine_type: Option<RType>,
    pub debug_level: Option<DbgLvl>,
    pub num_iter: Option<i32>,
    pub num_cuts: Option<i32>,
    pub seed: Option<i32>,
    pub minimize_connectivity: Option<bool>,
    pub force_contiguous: Option<bool>,
    pub compress: Option<bool>,
    pub cc_order: Option<bool>,
    pub prune_factor: Option<i32>,
    pub num_separators: Option<i32>,
    pub imbalance_factor: Option<i32>,
    pub(crate) numbering: Option<Numbering>,
    pub disable_2hop: Option<bool>,
}
```

Each field uses `Option<T>` — `None` means "use default" (replaces the `-1` sentinel). `#[derive(Default)]` gives us all-`None` for free.

**No trait at all.** The enum types (`PType`, `ObjType`, `CType`, etc.) and newtype structs (`Seed`, `NCuts`, etc.) are kept as value types but no longer implement any trait. They're just used as field values on `Options`.

**Builder methods on `Graph` and `Mesh`** — one method per option, each taking the appropriate type and returning `Self`:

```rust
// src/lib.rs
impl<'a> Graph<'a> {
    pub fn seed(mut self, seed: i32) -> Self {
        self.options.seed = Some(seed);
        self
    }
    pub fn obj_type(mut self, obj_type: ObjType) -> Self {
        self.options.obj_type = Some(obj_type);
        self
    }
    pub fn num_cuts(mut self, n: i32) -> Self {
        self.options.num_cuts = Some(n);
        self
    }
    // ... etc for all options
}
```

**Usage becomes:**

```rust
// Before:
Graph::new(1, 2, &xadj, &adj).unwrap()
    .set_option(metis_clone::option::Seed(42))
    .part_kway(&mut part)

// After:
Graph::new(1, 2, &xadj, &adj).unwrap()
    .seed(42)
    .part_kway(&mut part)
```

**`Control::new`** takes `&Options` instead of `&[Idx; NOPTIONS]`:

```rust
let seed = options.seed.unwrap_or(-1);
let ncuts = options.num_cuts.unwrap_or(1);
```

**pmetis.rs internal construction** becomes:

```rust
let mut options = Options::default();
options.num_cuts = Some(ncuts);
options.num_iter = Some(ctrl.num_iter);
let mut pctrl = Control::new(&options, graph.num_constraints, nparts, false);
```

### What gets removed

- `Opt` trait, `private::Sealed` trait, and all `impl Opt for X` blocks in `src/option.rs`
- `NOPTIONS` constant from `src/types.rs`
- `set_option<O: Opt>` and `set_options(&[Idx; NOPTIONS])` from `Graph`/`Mesh`
- `NOPTIONS` from `pub use` in `src/lib.rs`
- The `Seed(i32)`, `NCuts(i32)`, etc. newtype structs can be removed since the builder methods take primitive types directly. The enum types (`PType`, `ObjType`, etc.) stay since they have real variants.

### Files to modify

| File | Change |
|------|--------|
| `src/option.rs` | Remove `Opt` trait, sealed trait, all `impl Opt` blocks. Remove newtype structs (`Seed`, `NCuts`, `NIter`, etc.). Add `Options` struct with `#[derive(Default)]`. Keep enum types (`PType`, `ObjType`, `CType`, `IpType`, `RType`, `Numbering`) and `DbgLvl`. |
| `src/types.rs` | Remove `NOPTIONS` constant |
| `src/lib.rs` | Change `Graph.options` and `Mesh.options` from `[Idx; NOPTIONS]` to `Options`. Replace `set_option`/`set_options` with per-option builder methods. Remove `NOPTIONS` from `pub use`. Remove `use option::Opt`. |
| `src/ctrl.rs` | Change `Control::new` signature from `&[Idx; NOPTIONS]` to `&Options`. Replace `get(idx, default)` closure with `options.field.unwrap_or(default)`. Remove `use crate::option::{self, Opt}`. |
| `src/partition/kway.rs` | Pass `&graph.options` (now `Options` type) to `Control::new` |
| `src/partition/pmetis.rs` | Replace `[-1; NOPTIONS]` + magic index writes with `Options::default()` + field assignments. Remove `use crate::types::NOPTIONS`. |
| `src/mesh/meshpart.rs` | Pass `&mesh.options` (now `Options` type) to `Control::new` |
| `tests/comparison.rs` | Update test calls from `.set_option(metis_clone::option::Seed(seed))` to `.seed(seed)` |

---

## Phase 2: Eliminate Unsafe Code (Arena-Based Coarsening)

**Risk: High. Changes the core data structure threading through coarsen/refine. Must be carefully tested.**

Replace the linked-list coarsening chain (`coarser: Option<Box<GraphData>>` + `finer: *mut GraphData`) with a `Vec<GraphData>` arena. This eliminates all `unsafe` blocks except the two `new_unchecked` methods.

### Current pattern (6 unsafe blocks)

The coarsening chain is a linked list where each graph owns its coarser neighbor (`Box`) and has a raw back-pointer to its finer neighbor. Walking the chain requires `unsafe` pointer dereferences:

- `src/graph/coarsen.rs:31` — building the chain
- `src/partition/mod.rs:42,112-117` — finding coarsest, walking chain
- `src/partition/pmetis.rs:241,334-340` — same pattern for recursive bisection
- `src/partition/kwayrefine.rs:24,38,57,64` — walking coarser-to-finer
- `src/partition/refine2way.rs:25,41` — walking coarser-to-finer

### New design

**`Vec<GraphData>` arena**: `coarsen_graph` returns a `Vec<GraphData>` where index 0 is one level coarser than the original and the last index is the coarsest. The coarsening chain is just indices into this vec.

**Key changes:**

1. **`coarsen_graph`** returns `Vec<GraphData>` (the coarsened levels, excluding the original). The original graph is mutated in place (it gets `coarse_map` and `matching` set on it). The returned vec is ordered finest-first (index 0 = one level coarser than original, last = coarsest).

2. **Remove from `GraphData`**: `coarser`, `finer` fields.

3. **`refine_kway` / `refine_2way`**: Take `&mut [GraphData]` (the arena slice) plus `&mut GraphData` (the original graph). Walk from coarsest to finest using index iteration. Use `split_at_mut` when simultaneous access to coarser and finer levels is needed during projection.

4. **`get_coarsest_ptr`**: Replaced by `arena.last()`.

5. **`mlevel_kway_partitioning`**: The coarsen->init->refine flow becomes:
   ```rust
   let arena = coarsen_graph(ctrl, graph);
   let coarsest = arena.last_mut().unwrap();
   init_kway_partition(ctrl, coarsest, nparts);
   refine_kway(ctrl, graph, &mut arena);
   // graph now has the refined partition
   ```

6. **`multilevel_bisect`** (in pmetis.rs): Same pattern — returns arena, refines through it.

### Files to modify

| File | Change |
|------|--------|
| `src/graph/mod.rs` | Remove `coarser` and `finer` fields from `GraphData` |
| `src/graph/coarsen.rs` | Return `Vec<GraphData>` instead of building linked list |
| `src/graph/contract.rs` | Return `GraphData` directly instead of `Box<GraphData>` |
| `src/partition/mod.rs` | Rewrite `mlevel_kway_partitioning` to use arena. Remove `get_coarsest_ptr`. |
| `src/partition/kwayrefine.rs` | Rewrite `refine_kway` to take arena slice. Remove all `unsafe`. |
| `src/partition/refine2way.rs` | Rewrite `refine_2way` to take arena slice. Remove all `unsafe`. |
| `src/partition/pmetis.rs` | Rewrite `multilevel_bisect` to use arena. Remove `get_coarsest_ptr`. |

### Critical constraint

The projection step needs simultaneous mutable access to two adjacent levels (coarser reads `coarse_map`, finer writes `partition`). This is handled by `split_at_mut`:

```rust
// Walk from coarsest (last) to finest (first), then to original graph
for i in (1..arena.len()).rev() {
    let (left, right) = arena.split_at_mut(i);
    let finer = &mut left[i - 1];
    let coarser = &mut right[0];
    project_kway_partition(ctrl, finer, coarser);
    // ... balance and refine finer ...
}
// Final projection to original graph
let coarser = &arena[0];
project_kway_partition_to_original(ctrl, graph, coarser);
```

---

## Phase 3: Clean Up Sentinel Values

**Risk: Low-Medium. Mechanical replacements within internal code.**

Replace the most prominent `-1` sentinel patterns with `Option` types or dedicated constants.

### Targets (prioritized)

1. **`boundary_map: Vec<Idx>`** — Uses `-1` to mean "not on boundary". Replace with named constant `const NOT_ON_BOUNDARY: Idx = -1;` in `src/graph/mod.rs`. Replace bare `-1` comparisons with `NOT_ON_BOUNDARY`.

2. **`neighbor_offset: i32`** in `KwayCutInfo` — Uses `-1` to mean "no neighbors allocated". Replace with named constant `const NO_NEIGHBOR_OFFSET: i32 = -1;`.

3. **`moved: Vec<Idx>`** in FM refinement — Uses `-1` as "not moved". Named constant `const NOT_MOVED: Idx = -1;`.

4. **`update_index: Vec<Idx>`** in kway FM — Uses `-1` as "not in update list". Named constant `const NOT_IN_UPDATE_LIST: Idx = -1;`.

5. **`num_init_parts: Idx = -1`** in Control — Uses `-1` as "not yet computed". Change to `Option<Idx>`.

### Files to modify

- `src/graph/mod.rs` — Constants, use in `add_to_boundary`/`remove_from_boundary`/`alloc_2way`
- `src/partition/fm.rs` — `NOT_MOVED` constant
- `src/partition/kwayfm.rs` — `NOT_IN_UPDATE_LIST` constant
- `src/ctrl.rs` — `num_init_parts: Option<Idx>`

---

## Phase 4: Remove `unsafe` from Public API

**Risk: Minimal.**

The two `new_unchecked` methods on `Graph` and `Mesh` are marked `unsafe` but don't actually perform any unsafe operations — they just skip validation. Remove the `unsafe` keyword and document them as "unchecked" without requiring an unsafe block.

### Files to modify

- `src/lib.rs` — Remove `unsafe` from `Graph::new_unchecked` and `Mesh::new_unchecked`

---

## Execution Order

1. **Phase 1** (Options Rewrite) — independent, medium risk
2. **Phase 3** (Sentinel Cleanup) — independent, low risk
3. **Phase 4** (Remove unsafe from public API) — trivial, do anytime
4. **Phase 2** (Arena Coarsening) — highest risk, do last when everything else is stable

---

## Verification

After each phase, run:
```
cargo test
```

All 61 tests must pass (31 PQueue + 17 comparison + 13 debug). The comparison tests validate bit-identical output against C METIS, so any behavioral regression will be caught immediately.
