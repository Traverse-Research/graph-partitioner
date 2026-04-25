#![allow(dead_code)]

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

/// 16×16 quad-grid mesh dual graph (`ncommon = 2` shared-edge dual), captured
/// from the `breda-mesh` meshlet builder's `partition_mesh` call with the test
/// fixture `grid_mesh(16, 16)` after `spatial_sort_triangles` +
/// `meshopt::optimize_vertex_cache_in_place`. The first connected component
/// (450 triangles) is partitioned into 4 parts via `Graph::part_recursive`.
///
/// On this input, `graph_partitioner::Graph::part_recursive` and
/// `metis::Graph::part_recursive` produce different partition assignments,
/// which propagates downstream as a one-triangle drift in the Nanite golden
/// snapshot. The port should be bit-identical with C metis.
#[allow(dead_code)]
pub fn meshlet_grid_dual_450t_4p() -> (Vec<i32>, Vec<i32>) {
    let xadj: Vec<i32> = vec![
        0, 1, 4, 6, 8, 11, 13, 16, 19, 21, 24, 26, 29, 32, 35, 38, 40, 43, 45, 48, 51, 54, 57, 60,
        63, 65, 68, 70, 73, 76, 79, 81, 84, 87, 90, 93, 96, 99, 101, 104, 107, 110, 113, 116, 119,
        122, 124, 127, 130, 133, 136, 139, 142, 145, 147, 150, 153, 156, 159, 162, 165, 168, 170,
        173, 176, 179, 182, 185, 188, 191, 193, 196, 199, 202, 205, 208, 211, 214, 216, 219, 222,
        225, 228, 231, 234, 237, 239, 242, 245, 248, 251, 254, 257, 260, 262, 265, 268, 271, 274,
        277, 280, 283, 285, 287, 290, 293, 296, 299, 301, 304, 307, 310, 313, 316, 319, 321, 324,
        327, 330, 333, 336, 339, 342, 345, 347, 350, 352, 355, 358, 361, 364, 367, 370, 373, 375,
        378, 381, 384, 387, 390, 393, 396, 398, 401, 404, 407, 410, 413, 416, 419, 421, 424, 427,
        430, 433, 436, 439, 442, 444, 447, 450, 453, 456, 459, 462, 465, 467, 470, 473, 476, 479,
        482, 485, 488, 490, 493, 496, 499, 502, 505, 508, 511, 513, 516, 519, 522, 525, 528, 531,
        534, 536, 539, 542, 545, 548, 551, 554, 557, 559, 562, 563, 565, 568, 571, 574, 576, 579,
        582, 585, 587, 590, 593, 596, 598, 601, 604, 606, 609, 612, 615, 617, 620, 623, 625, 628,
        631, 634, 637, 640, 643, 645, 648, 651, 654, 657, 660, 663, 666, 669, 671, 674, 677, 680,
        683, 686, 689, 692, 695, 698, 701, 704, 707, 710, 713, 716, 719, 722, 725, 728, 731, 734,
        737, 740, 743, 746, 749, 752, 755, 758, 761, 764, 767, 770, 773, 776, 779, 782, 785, 788,
        791, 794, 797, 800, 803, 806, 809, 812, 815, 818, 821, 824, 827, 830, 833, 836, 839, 842,
        845, 848, 851, 854, 857, 860, 863, 866, 869, 872, 875, 878, 881, 884, 887, 890, 893, 896,
        899, 902, 905, 908, 911, 914, 917, 920, 923, 926, 929, 932, 935, 938, 941, 944, 947, 950,
        953, 956, 959, 962, 965, 968, 971, 974, 977, 980, 983, 986, 988, 991, 994, 996, 999, 1002,
        1005, 1008, 1011, 1013, 1016, 1019, 1022, 1025, 1028, 1030, 1033, 1036, 1039, 1042, 1044,
        1047, 1050, 1053, 1056, 1059, 1062, 1064, 1067, 1070, 1073, 1076, 1079, 1082, 1085, 1088,
        1090, 1093, 1096, 1099, 1102, 1105, 1107, 1110, 1113, 1116, 1119, 1122, 1125, 1127, 1129,
        1132, 1135, 1138, 1141, 1143, 1146, 1149, 1152, 1155, 1158, 1161, 1163, 1166, 1169, 1172,
        1175, 1178, 1181, 1184, 1187, 1189, 1192, 1195, 1198, 1201, 1204, 1207, 1210, 1213, 1216,
        1219, 1221, 1224, 1227, 1230, 1233, 1236, 1239, 1242, 1245, 1248, 1251, 1254, 1257, 1260,
        1263, 1266, 1269, 1272, 1275, 1278, 1281, 1284, 1287, 1290,
    ];
    let adjacency: Vec<i32> = vec![
        1, 0, 3, 2, 1, 4, 1, 6, 2, 7, 5, 4, 9, 3, 8, 7, 4, 6, 11, 6, 13, 5, 12, 10, 9, 16, 7, 14,
        12, 9, 11, 18, 8, 15, 14, 11, 13, 20, 13, 22, 10, 19, 17, 16, 25, 12, 21, 19, 16, 18, 27,
        14, 23, 21, 18, 20, 33, 15, 24, 23, 20, 22, 32, 22, 29, 17, 28, 26, 25, 343, 19, 34, 28,
        25, 27, 338, 24, 30, 31, 29, 36, 29, 32, 38, 23, 31, 35, 21, 35, 34, 27, 33, 43, 32, 33,
        40, 30, 37, 39, 36, 44, 31, 39, 41, 36, 38, 46, 35, 41, 42, 38, 40, 48, 40, 43, 51, 34,
        42, 337, 37, 45, 47, 44, 52, 39, 47, 49, 44, 46, 54, 41, 49, 50, 46, 48, 56, 48, 51, 59,
        42, 50, 329, 45, 53, 55, 52, 60, 47, 55, 57, 52, 54, 62, 49, 57, 58, 54, 56, 64, 56, 59,
        67, 50, 58, 321, 53, 61, 63, 60, 68, 55, 63, 65, 60, 62, 70, 57, 65, 66, 62, 64, 72, 64,
        67, 75, 58, 66, 313, 61, 69, 71, 68, 76, 63, 71, 73, 68, 70, 78, 65, 73, 74, 70, 72, 80,
        72, 75, 83, 66, 74, 305, 69, 77, 79, 76, 84, 71, 79, 81, 76, 78, 86, 73, 81, 82, 78, 80,
        88, 80, 83, 91, 74, 82, 303, 77, 85, 87, 84, 92, 79, 87, 89, 84, 86, 94, 81, 89, 90, 86,
        88, 96, 88, 91, 99, 82, 90, 294, 85, 93, 95, 92, 100, 87, 95, 97, 92, 94, 104, 89, 97, 98,
        94, 96, 109, 96, 99, 116, 90, 98, 287, 93, 101, 103, 100, 102, 101, 105, 100, 104, 106,
        95, 103, 108, 102, 106, 107, 103, 105, 110, 105, 112, 104, 109, 111, 97, 108, 115, 106,
        111, 113, 108, 110, 117, 107, 113, 114, 110, 112, 119, 112, 121, 109, 116, 118, 98, 115,
        282, 111, 118, 120, 115, 117, 131, 113, 120, 122, 117, 119, 128, 114, 122, 123, 119, 121,
        126, 121, 124, 123, 127, 125, 124, 132, 122, 129, 127, 124, 126, 134, 120, 130, 129, 126,
        128, 136, 128, 131, 138, 118, 281, 130, 125, 135, 133, 132, 140, 127, 137, 135, 132, 134,
        142, 129, 139, 137, 134, 136, 144, 130, 273, 139, 136, 138, 146, 133, 143, 141, 140, 148,
        135, 145, 143, 140, 142, 150, 137, 147, 145, 142, 144, 152, 139, 265, 147, 144, 146, 154,
        141, 151, 149, 148, 156, 143, 153, 151, 148, 150, 158, 145, 155, 153, 150, 152, 160, 147,
        257, 155, 152, 154, 162, 149, 159, 157, 156, 164, 151, 161, 159, 156, 158, 166, 153, 163,
        161, 158, 160, 168, 155, 249, 163, 160, 162, 170, 157, 167, 165, 164, 172, 159, 169, 167,
        164, 166, 174, 161, 171, 169, 166, 168, 176, 163, 247, 171, 168, 170, 178, 165, 175, 173,
        172, 180, 167, 177, 175, 172, 174, 182, 169, 179, 177, 174, 176, 184, 171, 236, 179, 176,
        178, 186, 173, 183, 181, 180, 188, 175, 185, 183, 180, 182, 190, 177, 187, 185, 182, 184,
        192, 179, 227, 187, 184, 186, 194, 181, 191, 189, 188, 196, 183, 193, 191, 188, 190, 203,
        185, 195, 193, 190, 192, 207, 187, 220, 195, 192, 194, 211, 189, 201, 197, 196, 198, 197,
        200, 199, 198, 198, 202, 196, 203, 202, 201, 204, 200, 191, 205, 201, 202, 206, 203, 207,
        206, 205, 208, 204, 193, 209, 205, 206, 210, 207, 211, 210, 209, 212, 208, 195, 214, 209,
        210, 213, 214, 215, 212, 211, 216, 213, 213, 217, 218, 220, 214, 218, 219, 215, 216, 223,
        217, 217, 221, 194, 225, 216, 222, 224, 219, 221, 228, 224, 226, 218, 223, 230, 221, 226,
        227, 220, 225, 232, 223, 186, 234, 225, 229, 231, 222, 228, 237, 231, 233, 224, 230, 239,
        228, 233, 235, 226, 232, 241, 230, 235, 236, 227, 234, 243, 232, 178, 245, 234, 238, 240,
        229, 237, 427, 240, 242, 231, 239, 430, 237, 242, 244, 233, 241, 255, 239, 244, 246, 235,
        243, 252, 241, 246, 247, 236, 245, 250, 243, 170, 248, 245, 249, 251, 247, 162, 256, 248,
        251, 253, 246, 248, 258, 250, 253, 254, 244, 250, 260, 252, 252, 262, 255, 254, 431, 242,
        257, 259, 249, 154, 264, 256, 259, 261, 251, 256, 266, 258, 261, 263, 253, 258, 268, 260,
        263, 435, 254, 260, 270, 262, 265, 267, 257, 146, 272, 264, 267, 269, 259, 264, 274, 266,
        269, 271, 261, 266, 276, 268, 271, 439, 263, 268, 278, 270, 273, 275, 265, 138, 280, 272,
        275, 277, 267, 272, 284, 274, 277, 279, 269, 274, 289, 276, 279, 443, 271, 276, 296, 278,
        281, 283, 273, 131, 282, 280, 116, 285, 281, 280, 286, 284, 283, 288, 275, 282, 287, 286,
        285, 290, 283, 99, 292, 285, 284, 291, 289, 288, 295, 277, 286, 293, 291, 290, 297, 288,
        287, 294, 293, 292, 299, 290, 91, 301, 292, 289, 298, 296, 295, 447, 279, 291, 300, 298,
        297, 311, 295, 293, 302, 300, 299, 308, 297, 294, 303, 302, 301, 306, 299, 83, 304, 301,
        303, 305, 307, 75, 312, 304, 302, 307, 309, 304, 314, 306, 300, 309, 310, 306, 316, 308,
        308, 318, 311, 298, 310, 449, 305, 313, 315, 67, 320, 312, 307, 315, 317, 312, 322, 314,
        309, 317, 319, 314, 324, 316, 310, 319, 382, 316, 326, 318, 313, 321, 323, 59, 328, 320,
        315, 323, 325, 320, 330, 322, 317, 325, 327, 322, 332, 324, 319, 327, 373, 324, 334, 326,
        321, 329, 331, 51, 336, 328, 323, 331, 333, 328, 340, 330, 325, 333, 335, 330, 350, 332,
        327, 335, 366, 332, 356, 334, 329, 337, 339, 43, 338, 336, 28, 337, 341, 336, 342, 340,
        331, 339, 348, 338, 343, 342, 339, 341, 346, 26, 341, 344, 343, 345, 344, 346, 347, 342,
        345, 349, 345, 351, 340, 349, 350, 346, 348, 352, 333, 348, 354, 347, 352, 353, 349, 351,
        355, 351, 357, 350, 355, 356, 352, 354, 358, 335, 354, 360, 353, 358, 359, 355, 357, 361,
        357, 362, 356, 361, 365, 358, 360, 363, 359, 363, 364, 361, 362, 367, 362, 369, 360, 366,
        368, 334, 365, 372, 363, 368, 370, 365, 367, 374, 364, 370, 371, 367, 369, 376, 369, 378,
        366, 373, 375, 326, 372, 381, 368, 375, 377, 372, 374, 383, 370, 377, 379, 374, 376, 389,
        371, 379, 380, 376, 378, 388, 378, 385, 373, 382, 384, 318, 381, 448, 375, 384, 390, 381,
        383, 417, 380, 387, 386, 385, 392, 385, 388, 396, 379, 391, 387, 377, 390, 391, 383, 389,
        408, 388, 389, 401, 386, 395, 393, 392, 394, 393, 397, 392, 396, 398, 387, 400, 395, 394,
        398, 399, 395, 402, 397, 397, 404, 396, 401, 403, 391, 407, 400, 398, 403, 405, 400, 409,
        402, 399, 405, 406, 402, 411, 404, 404, 413, 401, 408, 410, 390, 416, 407, 403, 410, 412,
        407, 418, 409, 405, 412, 414, 409, 420, 411, 406, 414, 415, 411, 422, 413, 413, 424, 408,
        417, 419, 384, 446, 416, 410, 419, 421, 416, 442, 418, 412, 421, 423, 418, 438, 420, 414,
        423, 425, 420, 434, 422, 415, 425, 426, 422, 429, 424, 424, 427, 426, 428, 238, 429, 430,
        427, 425, 432, 428, 428, 431, 240, 255, 433, 430, 433, 434, 429, 432, 435, 431, 423, 436,
        432, 262, 437, 433, 437, 438, 434, 436, 439, 435, 421, 440, 436, 270, 441, 437, 441, 442,
        438, 440, 443, 439, 419, 444, 440, 278, 445, 441, 445, 446, 442, 444, 447, 443, 417, 448,
        444, 296, 449, 445, 382, 449, 446, 311, 448, 447,
    ];
    (xadj, adjacency)
}

/// Second connected component of the same 16×16 quad-grid mesh as
/// `meshlet_grid_dual_450t_4p` — the welder splits the grid into two
/// components, and `partition_mesh` is invoked once per component. This
/// component is 247 elements partitioned into 2 parts via
/// `Graph::part_recursive`. Captured at the same call site.
#[allow(dead_code)]
pub fn meshlet_grid_dual_247t_2p() -> (Vec<i32>, Vec<i32>) {
    let xadj: Vec<i32> = vec![
        0, 3, 5, 8, 10, 13, 16, 18, 21, 24, 27, 29, 32, 35, 38, 40, 43, 46, 49, 51, 54, 57, 58, 60,
        63, 65, 68, 71, 74, 76, 79, 81, 84, 87, 89, 92, 94, 97, 100, 103, 105, 108, 111, 114, 117,
        120, 123, 126, 129, 132, 135, 138, 140, 142, 145, 147, 150, 152, 155, 158, 161, 164, 167,
        170, 172, 175, 177, 179, 181, 184, 187, 190, 193, 196, 199, 202, 205, 208, 211, 214, 217,
        220, 223, 226, 229, 232, 235, 238, 241, 244, 247, 250, 253, 256, 259, 262, 265, 268, 271,
        274, 277, 280, 283, 286, 289, 290, 293, 295, 297, 300, 302, 305, 308, 311, 313, 315, 318,
        320, 323, 326, 329, 331, 334, 337, 340, 342, 345, 348, 351, 353, 356, 359, 362, 364, 367,
        370, 372, 375, 378, 380, 382, 385, 388, 391, 394, 397, 400, 403, 406, 409, 412, 415, 418,
        421, 423, 426, 429, 432, 435, 438, 441, 444, 447, 450, 453, 456, 459, 462, 465, 468, 471,
        474, 477, 480, 483, 485, 488, 491, 493, 496, 499, 502, 505, 507, 510, 513, 516, 519, 522,
        525, 527, 530, 533, 536, 539, 542, 545, 548, 551, 553, 556, 559, 562, 565, 568, 570, 573,
        576, 579, 582, 585, 588, 590, 592, 595, 598, 601, 604, 606, 609, 612, 615, 618, 621, 624,
        626, 629, 632, 635, 638, 641, 644, 647, 650, 652, 655, 658, 661, 664, 667, 670, 673, 676,
        679, 682, 684, 687, 690,
    ];
    let adjacency: Vec<i32> = vec![
        67, 72, 1, 0, 2, 1, 4, 3, 2, 5, 2, 75, 8, 3, 12, 6, 5, 9, 76, 77, 8, 4, 7, 11, 6, 16, 10,
        9, 13, 8, 41, 12, 5, 11, 15, 10, 20, 14, 13, 17, 12, 38, 16, 9, 15, 19, 14, 26, 18, 17, 22,
        16, 32, 20, 13, 19, 27, 23, 18, 23, 22, 24, 21, 23, 25, 26, 28, 24, 17, 27, 25, 20, 29, 26,
        25, 29, 27, 30, 28, 31, 29, 32, 33, 30, 19, 36, 31, 31, 34, 35, 37, 33, 34, 39, 37, 38, 32,
        36, 44, 34, 15, 40, 36, 35, 46, 41, 43, 38, 11, 42, 40, 77, 81, 41, 40, 45, 44, 43, 47, 37,
        81, 98, 43, 47, 244, 39, 44, 48, 46, 47, 97, 243, 139, 143, 50, 49, 51, 146, 50, 53, 53,
        55, 51, 52, 57, 55, 56, 52, 54, 59, 54, 61, 53, 58, 60, 57, 146, 92, 55, 60, 62, 57, 59,
        90, 56, 62, 63, 59, 61, 69, 61, 64, 63, 70, 65, 64, 68, 68, 67, 66, 0, 65, 71, 66, 62, 86,
        70, 64, 69, 73, 68, 74, 72, 0, 71, 75, 70, 83, 74, 71, 73, 78, 72, 76, 4, 75, 78, 7, 7, 79,
        42, 74, 80, 76, 80, 85, 77, 78, 82, 79, 42, 84, 45, 83, 89, 80, 73, 87, 82, 85, 100, 81,
        79, 88, 84, 69, 90, 87, 86, 94, 83, 89, 102, 85, 82, 96, 88, 60, 91, 86, 90, 92, 93, 58,
        147, 91, 91, 155, 94, 93, 95, 87, 94, 157, 96, 95, 103, 89, 98, 241, 48, 45, 99, 97, 100,
        239, 98, 84, 101, 99, 102, 237, 100, 88, 245, 101, 160, 246, 96, 105, 104, 106, 107, 105,
        108, 105, 110, 106, 109, 111, 108, 112, 107, 111, 114, 108, 110, 115, 109, 113, 118, 112,
        117, 110, 153, 111, 118, 151, 117, 119, 113, 116, 121, 112, 115, 122, 116, 120, 125, 119,
        123, 117, 122, 126, 118, 121, 150, 120, 124, 129, 123, 127, 119, 126, 130, 121, 125, 169,
        124, 128, 133, 127, 131, 123, 130, 133, 125, 129, 149, 128, 132, 136, 131, 134, 127, 129,
        137, 132, 135, 141, 134, 140, 131, 137, 142, 133, 136, 161, 139, 140, 138, 49, 135, 138,
        143, 134, 142, 145, 136, 141, 148, 49, 140, 144, 143, 145, 146, 141, 144, 147, 50, 144, 58,
        92, 145, 155, 142, 158, 156, 130, 161, 165, 122, 152, 173, 115, 153, 152, 150, 151, 154,
        114, 151, 152, 174, 175, 147, 156, 93, 148, 155, 157, 156, 159, 95, 148, 162, 159, 157,
        158, 160, 159, 164, 103, 137, 149, 162, 158, 161, 163, 162, 166, 164, 160, 163, 200, 149,
        168, 166, 163, 165, 167, 166, 171, 191, 165, 169, 170, 126, 173, 168, 168, 176, 172, 167,
        172, 184, 170, 179, 171, 150, 169, 175, 154, 177, 154, 173, 176, 170, 175, 178, 174, 178,
        176, 177, 180, 172, 180, 183, 178, 179, 181, 180, 182, 185, 181, 187, 179, 184, 186, 171,
        183, 190, 181, 186, 188, 183, 185, 192, 182, 188, 189, 185, 187, 194, 187, 196, 184, 191,
        193, 167, 190, 199, 186, 193, 195, 190, 192, 201, 188, 195, 197, 192, 194, 207, 189, 197,
        198, 194, 196, 206, 196, 203, 191, 200, 202, 164, 199, 246, 193, 202, 208, 199, 201, 235,
        198, 205, 204, 203, 210, 203, 206, 214, 197, 209, 205, 195, 208, 209, 201, 207, 226, 206,
        207, 219, 204, 213, 211, 210, 212, 211, 215, 210, 214, 216, 205, 218, 213, 212, 216, 217,
        213, 220, 215, 215, 222, 214, 219, 221, 209, 225, 218, 216, 221, 223, 218, 227, 220, 217,
        223, 224, 220, 229, 222, 222, 231, 219, 226, 228, 208, 234, 225, 221, 228, 230, 225, 236,
        227, 223, 230, 232, 227, 238, 229, 224, 232, 233, 229, 240, 231, 231, 242, 226, 235, 237,
        202, 245, 234, 228, 237, 239, 101, 234, 236, 230, 239, 241, 99, 236, 238, 232, 241, 243,
        97, 238, 240, 233, 243, 244, 48, 240, 242, 46, 242, 235, 246, 102, 103, 200, 245,
    ];
    (xadj, adjacency)
}
