`graph-partitioner`
========

[![Actions Status](https://github.com/Traverse-Research/graph-partitioner/workflows/Continuous%20integration/badge.svg)](https://github.com/Traverse-Research/graph-partitioner/actions)
[![Latest version](https://img.shields.io/crates/v/graph-partitioner.svg)](https://crates.io/crates/graph-partitioner)
[![Documentation](https://docs.rs/graph-partitioner/badge.svg)](https://docs.rs/graph-partitioner)
![MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![Apache](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
[![Contributor Covenant](https://img.shields.io/badge/contributor%20covenant-v1.4%20adopted-ff69b4.svg)](../master/CODE_OF_CONDUCT.md)

[![Banner](banner.png)](https://traverseresearch.nl)

A pure Rust reimplementation of the [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) graph partitioning algorithms, generated with the help of an LLM from the original [C source code](https://github.com/KarypisLab/METIS). Produces bit-identical output to C METIS 5.1.0.

Compared to the C original, this crate:

- Is written in safe, idiomatic Rust with no C dependencies
- Fixes a thread-safety issue in the original METIS RNG, which uses global mutable state
- Provides a builder API for ergonomic graph and mesh partitioning

- [Documentation](https://docs.rs/graph-partitioner)

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
graph-partitioner = "0.1.0"
```

### Graph partitioning

```rust
let mut partition = vec![0i32; num_vertices];
graph_partitioner::Graph::new(1, num_parts, &xadj, &adjacency)
    .unwrap()
    .seed(42)
    .part_kway(&mut partition)
    .unwrap();
```

### Mesh partitioning

```rust
let mut epart = vec![0i32; num_elements];
let mut npart = vec![0i32; num_nodes];
graph_partitioner::Mesh::new(num_parts, &element_offsets, &element_indices)
    .unwrap()
    .seed(42)
    .part_dual(&mut epart, &mut npart)
    .unwrap();
```

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
