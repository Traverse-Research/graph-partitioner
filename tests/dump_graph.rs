mod fixtures;
use fixtures::*;

#[test]
fn dump_sphere_dual_graph() {
    let (eoff, eind) = uv_sphere_mesh(8, 16);
    let ne = (eoff.len() - 1) as i32;
    let nn = eind.iter().copied().max().unwrap_or(-1) + 1;

    let (xadj, adjacency) = graph_partitioner::create_graph_dual(ne, nn, &eoff, &eind, 1);

    let n = xadj.len() - 1;
    eprintln!("// nvtxs={}, nedges={}", n, adjacency.len());
    eprintln!("static idx_t nvtxs = {};", n);
    eprintln!("static idx_t xadj[] = {{{}}};",
        xadj.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
    eprintln!("static idx_t adjncy[] = {{{}}};",
        adjacency.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
    eprintln!("// xadj len={}, adjncy len={}", xadj.len(), adjacency.len());
}
