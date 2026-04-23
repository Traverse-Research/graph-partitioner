/// 2 vertices, 1 edge: 0--1
pub fn trivial_2v() -> (Vec<i32>, Vec<i32>) {
    let xadj = vec![0, 1, 2];
    let adjncy = vec![1, 0];
    (xadj, adjncy)
}

/// Path graph: 0--1--2--3--4
pub fn path_5v() -> (Vec<i32>, Vec<i32>) {
    let xadj = vec![0, 1, 3, 5, 7, 8];
    let adjncy = vec![1, 0, 2, 1, 3, 2, 4, 3];
    (xadj, adjncy)
}

/// 3x5 grid graph (15 vertices)
///  0- 1- 2- 3- 4
///  |  |  |  |  |
///  5- 6- 7- 8- 9
///  |  |  |  |  |
/// 10-11-12-13-14
pub fn grid_3x5() -> (Vec<i32>, Vec<i32>) {
    let rows = 3;
    let cols = 5;
    build_grid(rows, cols)
}

/// 10x10 grid graph (100 vertices)
pub fn grid_10x10() -> (Vec<i32>, Vec<i32>) {
    build_grid(10, 10)
}

fn build_grid(rows: usize, cols: usize) -> (Vec<i32>, Vec<i32>) {
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjncy = Vec::new();

    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if r > 0 {
                adjncy.push((v - cols) as i32);
            }
            if c > 0 {
                adjncy.push((v - 1) as i32);
            }
            if c + 1 < cols {
                adjncy.push((v + 1) as i32);
            }
            if r + 1 < rows {
                adjncy.push((v + cols) as i32);
            }
            xadj[v + 1] = adjncy.len() as i32;
        }
    }

    (xadj, adjncy)
}

/// The metis-rs example graph (6 vertices, irregular)
/// Adjacency:
///   0: [1, 3]
///   1: [0, 2, 3]
///   2: [1, 4]
///   3: [0, 1, 4, 5]
///   4: [2, 3, 5]
///   5: [3, 4]
pub fn irregular_6v() -> (Vec<i32>, Vec<i32>) {
    let xadj = vec![0, 2, 5, 7, 11, 14, 16];
    let adjncy = vec![1, 3, 0, 2, 3, 1, 4, 0, 1, 4, 5, 2, 3, 5, 3, 4];
    (xadj, adjncy)
}

/// Simple 2D triangle mesh (4 triangles forming a square)
///  0---1---2
///  |\ | /|
///  | \|/ |
///  3---4---5
///
/// Elements: [0,1,4], [0,4,3], [1,2,4], [2,5,4]
pub fn tri_mesh() -> (Vec<i32>, Vec<i32>, i32) {
    let eptr = vec![0, 3, 6, 9, 12];
    let eind = vec![0, 1, 4, 0, 4, 3, 1, 2, 4, 2, 5, 4];
    let nn = 6;
    (eptr, eind, nn)
}

/// Simple 2D quad mesh (2 quads)
///  0---1---2
///  |   |   |
///  3---4---5
///
/// Elements: [0,1,4,3], [1,2,5,4]
pub fn quad_mesh() -> (Vec<i32>, Vec<i32>, i32) {
    let eptr = vec![0, 4, 8];
    let eind = vec![0, 1, 4, 3, 1, 2, 5, 4];
    let nn = 6;
    (eptr, eind, nn)
}
