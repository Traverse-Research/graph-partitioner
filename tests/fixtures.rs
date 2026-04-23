/// 2 vertices, 1 edge: 0--1
pub fn trivial_2v() -> (Vec<i32>, Vec<i32>) {
    let xadj = vec![0, 1, 2];
    let adjacency = vec![1, 0];
    (xadj, adjacency)
}

/// Path graph: 0--1--2--3--4
pub fn path_5v() -> (Vec<i32>, Vec<i32>) {
    let xadj = vec![0, 1, 3, 5, 7, 8];
    let adjacency = vec![1, 0, 2, 1, 3, 2, 4, 3];
    (xadj, adjacency)
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

pub fn build_grid(rows: usize, cols: usize) -> (Vec<i32>, Vec<i32>) {
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjacency = Vec::new();

    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if r > 0 {
                adjacency.push((v - cols) as i32);
            }
            if c > 0 {
                adjacency.push((v - 1) as i32);
            }
            if c + 1 < cols {
                adjacency.push((v + 1) as i32);
            }
            if r + 1 < rows {
                adjacency.push((v + cols) as i32);
            }
            xadj[v + 1] = adjacency.len() as i32;
        }
    }

    (xadj, adjacency)
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
    let adjacency = vec![1, 3, 0, 2, 3, 1, 4, 0, 1, 4, 5, 2, 3, 5, 3, 4];
    (xadj, adjacency)
}

/// Simple 2D triangle mesh (4 triangles forming a square)
///  0---1---2
///  |\ | /|
///  | \|/ |
///  3---4---5
///
/// Elements: [0,1,4], [0,4,3], [1,2,4], [2,5,4]
pub fn tri_mesh() -> (Vec<i32>, Vec<i32>, i32) {
    let element_offsets = vec![0, 3, 6, 9, 12];
    let element_indices = vec![0, 1, 4, 0, 4, 3, 1, 2, 4, 2, 5, 4];
    let nn = 6;
    (element_offsets, element_indices, nn)
}

/// Simple 2D quad mesh (2 quads)
///  0---1---2
///  |   |   |
///  3---4---5
///
/// Elements: [0,1,4,3], [1,2,5,4]
pub fn quad_mesh() -> (Vec<i32>, Vec<i32>, i32) {
    let element_offsets = vec![0, 4, 8];
    let element_indices = vec![0, 1, 4, 3, 1, 2, 5, 4];
    let nn = 6;
    (element_offsets, element_indices, nn)
}

// ========= Procedural mesh fixtures =========

/// Triangulated grid mesh: rows x cols quads, each split into 2 triangles.
/// Vertices are laid out row-major: (rows+1) * (cols+1) vertices.
/// Returns (element_offsets, element_indices).
pub fn tri_grid_mesh(rows: usize, cols: usize) -> (Vec<i32>, Vec<i32>) {
    let num_tris = 2 * rows * cols;
    let mut element_offsets = Vec::with_capacity(num_tris + 1);
    let mut element_indices = Vec::with_capacity(num_tris * 3);
    element_offsets.push(0);

    let vcols = cols + 1;
    for r in 0..rows {
        for c in 0..cols {
            let v00 = (r * vcols + c) as i32;
            let v10 = v00 + 1;
            let v01 = ((r + 1) * vcols + c) as i32;
            let v11 = v01 + 1;
            // Lower-left triangle
            element_indices.extend_from_slice(&[v00, v10, v01]);
            element_offsets.push(element_indices.len() as i32);
            // Upper-right triangle
            element_indices.extend_from_slice(&[v10, v11, v01]);
            element_offsets.push(element_indices.len() as i32);
        }
    }
    (element_offsets, element_indices)
}

/// Triangle strip mesh: n triangles in a zigzag strip.
/// Vertices: 0..=(n+1) along top, (n+2)..=(2n+2) along bottom (for even-width strip).
/// Returns (element_offsets, element_indices).
pub fn tri_strip_mesh(n: usize) -> (Vec<i32>, Vec<i32>) {
    // Use a 1-row x (n/2) grid for even n, or build explicitly for any n.
    // Simpler: build a strip of n triangles sharing edges.
    // Vertices: row 0 has ceil((n+2)/2) verts, row 1 has floor((n+2)/2) verts.
    // Actually, simplest: use a 2-row strip = tri_grid_mesh(1, cols) where cols = n/2
    // But n might be odd. Let's just build it directly.
    let mut element_offsets = Vec::with_capacity(n + 1);
    let mut element_indices = Vec::with_capacity(n * 3);
    element_offsets.push(0);

    // Vertices numbered 0, 1, 2, ... along a zigzag.
    // Triangle i uses vertices i, i+1, i+2 but alternating orientation.
    // This is a classic triangle strip topology.
    for i in 0..n {
        if i % 2 == 0 {
            element_indices.extend_from_slice(&[i as i32, (i + 1) as i32, (i + 2) as i32]);
        } else {
            element_indices.extend_from_slice(&[(i + 1) as i32, i as i32, (i + 2) as i32]);
        }
        element_offsets.push(element_indices.len() as i32);
    }
    (element_offsets, element_indices)
}

/// Triangle mesh with sparse (non-contiguous) vertex IDs.
/// Generates n_tris triangles but vertex IDs are multiplied by `stride` (default 7),
/// so there are gaps in the ID space. Tests the remapping path.
/// Returns raw indices as u32 (matching the production partition_mesh signature).
pub fn sparse_tri_mesh(n_tris: usize) -> Vec<u32> {
    let stride: u32 = 7;
    let (_eptr, eindices) = tri_strip_mesh(n_tris);
    eindices.iter().map(|&v| (v as u32) * stride).collect()
}

/// UV sphere mesh. lat_segments rings from pole to pole, lon_segments divisions around.
/// Top pole = vertex 0, bottom pole = last vertex.
/// Returns (element_offsets, element_indices).
pub fn uv_sphere_mesh(lat_segments: usize, lon_segments: usize) -> (Vec<i32>, Vec<i32>) {
    assert!(lat_segments >= 2 && lon_segments >= 3);

    // Vertices:
    //   0 = north pole
    //   1..=(lon * (lat-1)) = ring vertices, row-major
    //   last = south pole
    let num_ring_verts = lon_segments * (lat_segments - 1);
    let south_pole = (1 + num_ring_verts) as i32;

    let num_tris = 2 * lon_segments * (lat_segments - 1);
    let mut eoff = Vec::with_capacity(num_tris + 1);
    let mut eind = Vec::with_capacity(num_tris * 3);
    eoff.push(0);

    // Helper: ring vertex at (ring r [0-based from top], segment s)
    let rv = |r: usize, s: usize| -> i32 { (1 + r * lon_segments + (s % lon_segments)) as i32 };

    // Top cap: triangles connecting north pole to first ring
    for s in 0..lon_segments {
        eind.extend_from_slice(&[0, rv(0, s), rv(0, s + 1)]);
        eoff.push(eind.len() as i32);
    }

    // Middle bands: quads split into 2 triangles each
    for r in 0..(lat_segments - 2) {
        for s in 0..lon_segments {
            // Lower-left triangle
            eind.extend_from_slice(&[rv(r, s), rv(r + 1, s), rv(r, s + 1)]);
            eoff.push(eind.len() as i32);
            // Upper-right triangle
            eind.extend_from_slice(&[rv(r, s + 1), rv(r + 1, s), rv(r + 1, s + 1)]);
            eoff.push(eind.len() as i32);
        }
    }

    // Bottom cap: triangles connecting last ring to south pole
    let last_ring = lat_segments - 2;
    for s in 0..lon_segments {
        eind.extend_from_slice(&[rv(last_ring, s), south_pole, rv(last_ring, s + 1)]);
        eoff.push(eind.len() as i32);
    }

    (eoff, eind)
}

/// Torus mesh. major_segments around the ring, minor_segments around the tube cross-section.
/// Vertex (i, j) = i * minor_segments + j, both indices wrap around.
/// Returns (element_offsets, element_indices).
pub fn torus_mesh(major_segments: usize, minor_segments: usize) -> (Vec<i32>, Vec<i32>) {
    assert!(major_segments >= 3 && minor_segments >= 3);

    let num_tris = 2 * major_segments * minor_segments;
    let mut eoff = Vec::with_capacity(num_tris + 1);
    let mut eind = Vec::with_capacity(num_tris * 3);
    eoff.push(0);

    let v = |i: usize, j: usize| -> i32 {
        ((i % major_segments) * minor_segments + (j % minor_segments)) as i32
    };

    for i in 0..major_segments {
        for j in 0..minor_segments {
            eind.extend_from_slice(&[v(i, j), v(i + 1, j), v(i, j + 1)]);
            eoff.push(eind.len() as i32);
            eind.extend_from_slice(&[v(i, j + 1), v(i + 1, j), v(i + 1, j + 1)]);
            eoff.push(eind.len() as i32);
        }
    }

    (eoff, eind)
}

/// Two disjoint UV spheres. Second sphere's vertex IDs are offset so they don't overlap.
/// Returns (element_offsets, element_indices).
pub fn two_disjoint_spheres_mesh(segments: usize) -> (Vec<i32>, Vec<i32>) {
    let (eoff1, eind1) = uv_sphere_mesh(segments, segments);
    let (eoff2, eind2) = uv_sphere_mesh(segments, segments);

    let vertex_offset = *eind1.iter().max().unwrap() + 1;

    let mut eoff = eoff1;
    let mut eind = eind1;

    let base = *eoff.last().unwrap();
    for &o in &eoff2[1..] {
        eoff.push(base + o);
    }
    for &v in &eind2 {
        eind.push(v + vertex_offset);
    }

    (eoff, eind)
}

/// Mobius strip mesh. A rectangular strip of length_segments x width_segments quads,
/// with the last column connected to the first with a half-twist.
/// Returns (element_offsets, element_indices).
pub fn mobius_strip_mesh(length_segments: usize, width_segments: usize) -> (Vec<i32>, Vec<i32>) {
    assert!(length_segments >= 3 && width_segments >= 1);

    // Vertices: length_segments * (width_segments + 1)
    // Row i, column j => vertex i * (width_segments + 1) + j
    // The twist: when wrapping from column (length_segments-1) to column 0,
    // we flip the width index: j -> width_segments - j
    let w = width_segments + 1;
    let num_tris = 2 * length_segments * width_segments;
    let mut eoff = Vec::with_capacity(num_tris + 1);
    let mut eind = Vec::with_capacity(num_tris * 3);
    eoff.push(0);

    let v = |i: usize, j: usize| -> i32 {
        let col = i % length_segments;
        // For the wrap-around edge (i == length_segments), we use column 0 with flipped j
        if i >= length_segments {
            let flipped_j = width_segments - j;
            (0 * w + flipped_j) as i32
        } else {
            (col * w + j) as i32
        }
    };

    for i in 0..length_segments {
        for j in 0..width_segments {
            eind.extend_from_slice(&[v(i, j), v(i + 1, j), v(i, j + 1)]);
            eoff.push(eind.len() as i32);
            eind.extend_from_slice(&[v(i, j + 1), v(i + 1, j), v(i + 1, j + 1)]);
            eoff.push(eind.len() as i32);
        }
    }

    (eoff, eind)
}

/// Open cylinder mesh. height_segments rows, radial_segments columns, wrapping around the axis.
/// Has two boundary loops (top and bottom edges).
/// Returns (element_offsets, element_indices).
pub fn cylinder_mesh(height_segments: usize, radial_segments: usize) -> (Vec<i32>, Vec<i32>) {
    assert!(height_segments >= 1 && radial_segments >= 3);

    // Vertices: (height_segments + 1) * radial_segments
    // Row r, segment s => vertex r * radial_segments + (s % radial_segments)
    let num_tris = 2 * height_segments * radial_segments;
    let mut eoff = Vec::with_capacity(num_tris + 1);
    let mut eind = Vec::with_capacity(num_tris * 3);
    eoff.push(0);

    let v = |r: usize, s: usize| -> i32 {
        (r * radial_segments + (s % radial_segments)) as i32
    };

    for r in 0..height_segments {
        for s in 0..radial_segments {
            eind.extend_from_slice(&[v(r, s), v(r + 1, s), v(r, s + 1)]);
            eoff.push(eind.len() as i32);
            eind.extend_from_slice(&[v(r, s + 1), v(r + 1, s), v(r + 1, s + 1)]);
            eoff.push(eind.len() as i32);
        }
    }

    (eoff, eind)
}

/// Weighted grid graph: rows x cols vertices with non-uniform edge weights.
/// Horizontal edges have weight 1, vertical edges have weight 3.
/// Returns (xadj, adjncy, adjwgt).
pub fn weighted_grid_graph(rows: usize, cols: usize) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjncy = Vec::new();
    let mut adjwgt = Vec::new();

    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if r > 0 {
                adjncy.push((v - cols) as i32);
                adjwgt.push(3); // vertical = heavy
            }
            if c > 0 {
                adjncy.push((v - 1) as i32);
                adjwgt.push(1); // horizontal = light
            }
            if c + 1 < cols {
                adjncy.push((v + 1) as i32);
                adjwgt.push(1);
            }
            if r + 1 < rows {
                adjncy.push((v + cols) as i32);
                adjwgt.push(3);
            }
            xadj[v + 1] = adjncy.len() as i32;
        }
    }

    (xadj, adjncy, adjwgt)
}
