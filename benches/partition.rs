#[path = "../tests/fixtures.rs"]
mod fixtures;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput};
use fixtures::*;

// ========= Graph k-way benchmarks =========

fn bench_grid_100x100(c: &mut Criterion) {
    let (xadj, adjacency) = build_grid(100, 100);
    let nvertices = (xadj.len() - 1) as u64;

    let mut group = c.benchmark_group("grid_100x100");
    group.throughput(Throughput::Elements(nvertices));

    for &nparts in &[4, 16] {
        group.bench_with_input(BenchmarkId::new("rust", nparts), &nparts, |b, &nparts| {
            b.iter_batched(
                || vec![0i32; nvertices as usize],
                |mut part| {
                    graph_partitioner::Graph::new(1, nparts, &xadj, &adjacency)
                        .unwrap()
                        .seed(42)
                        .part_kway(&mut part)
                        .unwrap();
                },
                BatchSize::LargeInput,
            );
        });

        group.bench_with_input(
            BenchmarkId::new("c_metis", nparts),
            &nparts,
            |b, &nparts| {
                b.iter_batched(
                    || vec![0i32; nvertices as usize],
                    |mut part| {
                        metis::Graph::new(1, nparts, &xadj, &adjacency)
                            .unwrap()
                            .set_option(metis::option::Seed(42))
                            .part_kway(&mut part)
                            .unwrap();
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_grid_300x300(c: &mut Criterion) {
    let (xadj, adjacency) = build_grid(300, 300);
    let nvertices = (xadj.len() - 1) as u64;

    let mut group = c.benchmark_group("grid_300x300");
    group.throughput(Throughput::Elements(nvertices));
    group.sample_size(10);

    for &nparts in &[4, 16, 64] {
        group.bench_with_input(BenchmarkId::new("rust", nparts), &nparts, |b, &nparts| {
            b.iter_batched(
                || vec![0i32; nvertices as usize],
                |mut part| {
                    graph_partitioner::Graph::new(1, nparts, &xadj, &adjacency)
                        .unwrap()
                        .seed(42)
                        .part_kway(&mut part)
                        .unwrap();
                },
                BatchSize::LargeInput,
            );
        });

        group.bench_with_input(
            BenchmarkId::new("c_metis", nparts),
            &nparts,
            |b, &nparts| {
                b.iter_batched(
                    || vec![0i32; nvertices as usize],
                    |mut part| {
                        metis::Graph::new(1, nparts, &xadj, &adjacency)
                            .unwrap()
                            .set_option(metis::option::Seed(42))
                            .part_kway(&mut part)
                            .unwrap();
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_weighted_grid_200x200(c: &mut Criterion) {
    let (xadj, adjncy, adjwgt) = weighted_grid_graph(200, 200);
    let nvertices = (xadj.len() - 1) as u64;

    let mut group = c.benchmark_group("weighted_grid_200x200");
    group.throughput(Throughput::Elements(nvertices));

    for &nparts in &[4, 16] {
        group.bench_with_input(BenchmarkId::new("rust", nparts), &nparts, |b, &nparts| {
            b.iter_batched(
                || vec![0i32; nvertices as usize],
                |mut part| {
                    graph_partitioner::Graph::new(1, nparts, &xadj, &adjncy)
                        .unwrap()
                        .set_edge_weights(&adjwgt)
                        .seed(42)
                        .part_kway(&mut part)
                        .unwrap();
                },
                BatchSize::LargeInput,
            );
        });

        group.bench_with_input(
            BenchmarkId::new("c_metis", nparts),
            &nparts,
            |b, &nparts| {
                b.iter_batched(
                    || vec![0i32; nvertices as usize],
                    |mut part| {
                        metis::Graph::new(1, nparts, &xadj, &adjncy)
                            .unwrap()
                            .set_adjwgt(&adjwgt)
                            .set_option(metis::option::Seed(42))
                            .part_kway(&mut part)
                            .unwrap();
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    group.finish();
}

// ========= Mesh dual benchmarks =========

fn bench_mesh_dual(
    c: &mut Criterion,
    group_name: &str,
    eptr: &[i32],
    eind: &[i32],
    nparts_list: &[i32],
    sample_size: Option<usize>,
) {
    let ne = (eptr.len() - 1) as u64;
    let nn = (*eind.iter().max().unwrap() + 1) as usize;

    let mut group = c.benchmark_group(group_name);
    group.throughput(Throughput::Elements(ne));
    if let Some(s) = sample_size {
        group.sample_size(s);
    }

    for &nparts in nparts_list {
        group.bench_with_input(BenchmarkId::new("rust", nparts), &nparts, |b, &nparts| {
            b.iter_batched(
                || (vec![0i32; ne as usize], vec![0i32; nn]),
                |(mut epart, mut npart)| {
                    graph_partitioner::Mesh::new(nparts, eptr, eind)
                        .unwrap()
                        .seed(42)
                        .part_dual(&mut epart, &mut npart)
                        .unwrap();
                },
                BatchSize::LargeInput,
            );
        });

        group.bench_with_input(
            BenchmarkId::new("c_metis", nparts),
            &nparts,
            |b, &nparts| {
                b.iter_batched(
                    || (vec![0i32; ne as usize], vec![0i32; nn]),
                    |(mut epart, mut npart)| {
                        metis::Mesh::new(nparts, eptr, eind)
                            .unwrap()
                            .set_option(metis::option::Seed(42))
                            .part_dual(&mut epart, &mut npart)
                            .unwrap();
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_tri_grid_100x100(c: &mut Criterion) {
    let (eptr, eind) = tri_grid_mesh(100, 100);
    bench_mesh_dual(c, "tri_grid_100x100", &eptr, &eind, &[4, 16], None);
}

fn bench_tri_grid_300x300(c: &mut Criterion) {
    let (eptr, eind) = tri_grid_mesh(300, 300);
    bench_mesh_dual(c, "tri_grid_300x300", &eptr, &eind, &[4, 16, 64], Some(10));
}

fn bench_sphere_64x128(c: &mut Criterion) {
    let (eptr, eind) = uv_sphere_mesh(64, 128);
    bench_mesh_dual(c, "sphere_64x128", &eptr, &eind, &[4, 16], None);
}

fn bench_sphere_128x256(c: &mut Criterion) {
    let (eptr, eind) = uv_sphere_mesh(128, 256);
    bench_mesh_dual(c, "sphere_128x256", &eptr, &eind, &[4, 16], None);
}

fn bench_torus_128x128(c: &mut Criterion) {
    let (eptr, eind) = torus_mesh(128, 128);
    bench_mesh_dual(c, "torus_128x128", &eptr, &eind, &[4, 16, 64], None);
}

criterion_group!(
    graph_benches,
    bench_grid_100x100,
    bench_grid_300x300,
    bench_weighted_grid_200x200,
);

criterion_group!(
    mesh_benches,
    bench_tri_grid_100x100,
    bench_tri_grid_300x300,
    bench_sphere_64x128,
    bench_sphere_128x256,
    bench_torus_128x128,
);

criterion_main!(graph_benches, mesh_benches);
