# Variable Rename Plan

Mark items with `[x]` to approve, `[ ]` to skip, or edit the "New" column.

---

## Struct & Type Renames

| Current | New | Notes |
|---------|-----|-------|
| `Ctrl` | `Control` | Used everywhere as the central config/state struct |
| `CKRInfo` | `KwayCutInfo` | Per-vertex cut refinement info |
| `CnbrInfo` | `NeighborPartInfo` | Neighbor-partition edge weight |
| `GraphData` | `GraphData` | Keep as-is (already descriptive enough) |

---

## GraphData Fields

| Current | New | Meaning |
|---------|-----|---------|
| `nvtxs` | `num_vertices` | Number of vertices |
| `nedges` | `num_edges` | Number of edges (total adjacency entries) |
| `ncon` | `num_constraints` | Number of balance constraints |
| `xadj` | `xadj` | Keep (standard CSR terminology, matches public API) |
| `adjncy` | `adjacency` | Adjacency list (CSR column indices) |
| `vwgt` | `vertex_weights` | Vertex weights |
| `vsize` | `vertex_sizes` | Vertex communication sizes |
| `adjwgt` | `edge_weights` | Edge weights |
| `where_` | `partition` | Partition assignment per vertex |
| `pwgts` | `part_weights` | Partition weight sums |
| `bndptr` | `boundary_map` | Reverse lookup: vertex -> position in boundary_list (-1 if not) |
| `bndind` | `boundary_list` | Dense list of boundary vertex IDs |
| `nbnd` | `num_boundary` | Number of boundary vertices |
| `id` | `internal_degree` | Internal degree (sum of edge weights to same partition) |
| `ed` | `external_degree` | External degree (sum of edge weights to other partitions) |
| `mincut` | `edge_cut` | Total edge cut of current partition |
| `ckrinfo` | `kway_refinement_info` | Per-vertex k-way refinement info |
| `cmap` | `coarse_map` | Vertex-to-coarse-vertex mapping |
| `match_` | `matching` | Matching partner for coarsening |
| `tvwgt` | `total_vertex_weight` | Total vertex weight per constraint |
| `invtvwgt` | `inv_total_vertex_weight` | Reciprocal of total vertex weight |
| `label` | `label` | Keep (already clear) |
| `coarser` | `coarser` | Keep (already clear) |
| `finer` | `finer` | Keep (already clear) |

---

## CKRInfo Fields

| Current | New | Meaning |
|---------|-----|---------|
| `id` | `internal_degree` | Internal degree |
| `ed` | `external_degree` | External degree |
| `nnbrs` | `num_neighbors` | Number of neighbor partitions |
| `inbr` | `neighbor_offset` | Offset into neighbor pool |

---

## CnbrInfo Fields

| Current | New | Meaning |
|---------|-----|---------|
| `pid` | `part_id` | Partition ID |
| `ed` | `external_degree` | External degree to this partition |

---

## Ctrl Fields

| Current | New | Meaning |
|---------|-----|---------|
| `optype` | `op_type` | Operation type (kway vs recursive) |
| `objtype` | `obj_type` | Objective type (cut vs volume) |
| `ctype` | `coarsen_type` | Coarsening algorithm (RM vs SHEM) |
| `iptype` | `init_part_type` | Initial partition method (grow vs random) |
| `rtype` | `refine_type` | Refinement algorithm (FM vs greedy) |
| `ncon` | `num_constraints` | Number of balance constraints |
| `nparts` | `num_parts` | Number of target partitions |
| `ncuts` | `num_cuts` | Number of bisection attempts |
| `nseps` | `num_separators` | Number of separators (not used in our flow) |
| `niter` | `num_iter` | Number of refinement iterations |
| `niparts` | `num_init_parts` | Number of initial partition attempts |
| `seed` | `seed` | Keep |
| `minconn` | `minimize_connectivity` | Minimize subdomain connectivity |
| `contig` | `force_contiguous` | Force contiguous partitions |
| `compress` | `compress` | Keep |
| `ccorder` | `cc_order` | Connected-component ordering |
| `pfactor` | `prune_factor` | Pruning factor |
| `ufactor` | `imbalance_factor` | Imbalance factor (integer, used to compute tolerance) |
| `no2hop` | `disable_2hop` | Disable 2-hop matching |
| `dbglvl` | `debug_level` | Debug verbosity level |
| `numflag` | `base_numbering` | 0-based vs 1-based numbering |
| `coarsen_to` | `coarsen_to` | Keep (already clear) |
| `ubfactors` | `imbalance_tols` | Imbalance tolerance per constraint |
| `tpwgts` | `target_part_weights` | Target partition weight fractions |
| `maxvwgt` | `max_vertex_weight` | Max vertex weight during coarsening |
| `pijbm` | `partition_ij_balance_multipliers` | Partition-constraint balance multipliers |
| `rng` | `rng` | Keep |
| `cnbrpool` | `neighbor_pool` | K-way neighbor info pool |
| `cnbrpool_pos` | `neighbor_pool_pos` | Next free position in neighbor pool |

---

## Function Parameters (heavily-used)

| Current | New | Where |
|---------|-----|-------|
| `ne` | `num_elements` | mesh functions |
| `nn` | `num_nodes` | mesh functions |
| `eptr` | `element_offsets` | Element CSR row offsets (matches public API param name) |
| `eind` | `element_indices` | Element node indices (matches public API param name) |
| `ncommon` | `min_common_nodes` | Minimum shared nodes for dual graph adjacency |
| `tpwgts` | `target_weights` | partition functions |
| `nparts` | `num_parts` | partition functions |
| `niparts` | `num_init_parts` | init partition |
| `niter` | `num_iter` | refinement functions |
| `ncon` | `num_constraints` | setup functions |
| `cnvtxs` | `coarse_num_vertices` | coarsening functions |
| `pijbm` | `partition_ij_balance_multipliers` | 2-way balance functions |

---

## Key Local Variables (across many files)

| Current | New | Context |
|---------|-----|---------|
| `higain` | `best_vertex` | FM refinement — vertex being moved |
| `tid` / `ted` | `internal_deg` / `external_deg` | Temporary degree accumulators |
| `rgain` | `scaled_gain` | Priority = ed/sqrt(nnbrs) - id |
| `itpwgts` | `int_target_weights` | Integer-floored target weights |
| `cg` | `coarse_graph` | Contracted graph being built |
| `hk` | `hash_slot` | Hash table probe position |
| `nbrmrk` | `neighbor_mark` | Neighbor-partition marker array |
| `nbrdom` | `neighbor_partitions` | Neighbor partition IDs |
| `nbrwgt` | `neighbor_weights` | Neighbor partition edge counts |
| `onemaxpwgt` | `max_target_weight` | Max weight for partition 1 in bisection |
| `oneminpwgt` | `min_target_weight` | Min weight for partition 1 in bisection |
| `zeromaxpwgt` | `max_zero_weight` | Max weight for partition 0 in random bisection |
| `kwgt` | `weight_delta` | Edge weight sign-flipped by partition side |
| `vstatus` | `vertex_status` | PQ status in greedy kway FM |
| `updptr` | `update_index` | Reverse-index: vertex -> position in update_list (-1 if not) |
| `updind` | `update_list` | List of vertices touched during pass |
| `nupd` | `num_updates` | Count of touched vertices |
| `ewgt` | `edge_weight` | Single edge weight in adjacency loop |

---

## Ctrl Methods

| Current | New |
|---------|-----|
| `cnbrpool_init` | `init_neighbor_pool` |
| `cnbrpool_reset` | `reset_neighbor_pool` |
| `cnbrpool_get_next` | `alloc_neighbor_info` |
| `setup_2way_bal_multipliers` | `setup_2way_balance_multipliers` |
| `setup_kway_bal_multipliers` | `setup_kway_balance_multipliers` |

---

## GraphData Methods

| Current | New |
|---------|-----|
| `bnd_insert` | `add_to_boundary` |
| `bnd_delete` | `remove_from_boundary` |
| `alloc_2way` | `alloc_2way` | Keep (already clear enough) |

---

## Constants

| Current | New |
|---------|-----|
| `SMALLNIPARTS` | `SMALL_NUM_INIT_PARTS` |
| `LARGENIPARTS` | `LARGE_NUM_INIT_PARTS` |
| `UNMATCHEDFOR2HOP` | `UNMATCHED_THRESHOLD_2HOP` |
| `BNDTYPE_REFINE` | `BOUNDARY_REFINE` |
| `BNDTYPE_BALANCE` | `BOUNDARY_BALANCE` |
| `VPQSTATUS_PRESENT` | `PQ_PRESENT` |
| `VPQSTATUS_EXTRACTED` | `PQ_EXTRACTED` |
| `VPQSTATUS_NOTPRESENT` | `PQ_NOT_PRESENT` |
| `OMODE_REFINE` | `MODE_REFINE` |
| `OMODE_BALANCE` | `MODE_BALANCE` |
