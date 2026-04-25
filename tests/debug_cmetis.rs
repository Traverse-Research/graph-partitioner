/// Debug test to call C METIS directly and compare intermediate results.
///
/// This calls METIS_PartGraphRecursive directly via FFI to understand what
/// the inner recursive bisection produces.
extern crate metis_sys;

use metis_sys::*;
use std::ptr;

/// Call METIS_PartGraphRecursive on the irregular_6v graph.
#[test]
fn debug_c_recursive_irregular_6v() {
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![1, 3, 0, 2, 3, 1, 4, 0, 1, 4, 5, 2, 3, 5, 3, 4];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;
    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 6];

    // Options: only set seed=-1 (default)
    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[9] = -1; // Seed = -1 (default)

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(), // vertex_weights
            ptr::null_mut(), // vertex_sizes
            ptr::null_mut(), // edge_weights
            &mut nparts,
            ptr::null_mut(), // target_part_weights
            ptr::null_mut(), // ubvec
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive (seed=-1): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );

    // Now try with seed=42
    options[9] = 42;
    let mut part2 = vec![0 as idx_t; 6];
    let mut edgecut2: idx_t = 0;
    let ret2 = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut2,
            part2.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive (seed=42): ret={}, edgecut={}, part={:?}",
        ret2, edgecut2, part2
    );
}

/// Call PartGraphRecursive on the coarsened graph to simulate what InitKWayPartitioning does.
/// The coarsened graph has reordered adjacency from CreateCoarseGraph.
#[test]
fn debug_c_recursive_coarsened_6v() {
    // The SINGLY-coarsened adjacency (from k-way CreateCoarseGraph with self-loop removal):
    // v0: [3, 1], v1: [3, 0, 2], v2: [4, 1], v3: [5, 0, 1, 4], v4: [5, 2, 3], v5: [4, 3]
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![3, 1, 3, 0, 2, 4, 1, 5, 0, 1, 4, 5, 2, 3, 4, 3];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    // Match the inner PMETIS parameters:
    options[8] = 4; // NCUTS = 4 (ctrl->niparts)
    options[9] = -1; // Seed = -1

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 6];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive on coarsened graph: ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
}

/// Replicate EXACTLY what InitKWayPartitioning does inside PartGraphKway.
/// This means calling PartGraphRecursive with:
/// - ncuts=4 (ctrl->nIparts for CoarsenTo=60=30*2)
/// - niter=10 (default)
/// - seed=-1 (inner PMETIS always gets seed=-1)
/// - ubvec computed from kway imbalance_tols: pow(1.0300499, 1/ln(2))
/// - The ORIGINAL graph (since kway can't coarsen 6v with CoarsenTo=60)
#[test]
fn debug_c_recursive_as_init_kway() {
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![1, 3, 0, 2, 3, 1, 4, 0, 1, 4, 5, 2, 3, 5, 3, 4];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    // Match InitKWayPartitioning exactly:
    options[7] = 10; // NITER = 10 (ctrl->niter default)
    options[8] = 4; // NCUTS = 4 (ctrl->nIparts)
                    // Seed stays -1 (default from METIS_SetDefaultOptions)
                    // OBJTYPE stays default (CUT)

    // Compute ubvec as InitKWayPartitioning does:
    // ubvec[i] = pow(ctrl->imbalance_tols[i], 1.0/log(nparts))
    // ctrl->imbalance_tols[0] = 1.0 + 0.001 * 30 + 0.0000499 = 1.0300499 (kway ufactor=30)
    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (2.0_f64).ln()) as f32;
    eprintln!("DEBUG ubvec_val = {}", ubvec_val);
    let mut ubvec = vec![ubvec_val];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 6];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive (as InitKWayPartitioning): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
}

/// DEFINITIVE test: Call PartGraphRecursive on the COARSENED graph
/// with EXACTLY the parameters that InitKWayPartitioning uses.
/// This should match what our init_kway_partition produces.
#[test]
fn debug_c_recursive_coarsened_with_ubvec() {
    // The SINGLY-coarsened adjacency (from k-way CreateCoarseGraph with self-loop removal):
    // v0: [3, 1], v1: [3, 0, 2], v2: [4, 1], v3: [5, 0, 1, 4], v4: [5, 2, 3], v5: [4, 3]
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![3, 1, 3, 0, 2, 4, 1, 5, 0, 1, 4, 5, 2, 3, 4, 3];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    // Match InitKWayPartitioning exactly:
    options[7] = 10; // NITER = 10 (ctrl->niter default)
    options[8] = 4; // NCUTS = 4 (ctrl->nIparts)
                    // Seed stays -1 (default)

    // Compute ubvec as InitKWayPartitioning does:
    // ubvec[i] = pow(ctrl->imbalance_tols[i], 1.0/log(nparts))
    // ctrl->imbalance_tols[0] = 1.0 + 0.001 * 30 + 0.0000499 = 1.0300499 (kway ufactor=30)
    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (2.0_f64).ln()) as f32;
    eprintln!("DEBUG ubvec_val = {}", ubvec_val);
    let mut ubvec = vec![ubvec_val];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 6];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive (coarsened + ubvec): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
}

/// Direct comparison of C rand() vs our Rust RNG (MSVC LCG)
#[test]
fn debug_compare_rng() {
    // Link to C's srand/rand (from MSVC runtime)
    extern "C" {
        fn srand(seed: u32);
        fn rand() -> i32;
    }

    // Test with seed 4321 (the default inner PMETIS seed)
    unsafe {
        srand(4321);
        eprint!("C rand() from seed 4321: ");
        for i in 0..20 {
            let v = rand();
            eprint!("{}", v);
            if i < 19 {
                eprint!(", ");
            }
        }
        eprintln!();
    }

    // Our Rust RNG (MSVC LCG: state = state * 214013 + 2531011; output = (state >> 16) & 0x7fff)
    let mut state: u32 = 4321;
    eprint!("Rust LCG from seed 4321:  ");
    for i in 0..20 {
        state = state.wrapping_mul(214013).wrapping_add(2531011);
        let v = ((state >> 16) & 0x7fff) as i32;
        eprint!("{}", v);
        if i < 19 {
            eprint!(", ");
        }
    }
    eprintln!();

    // Also compare after 12 calls (the point where GrowBisection starts)
    unsafe {
        srand(4321);
        for _ in 0..12 {
            rand();
        }
        eprint!("C rand()%6 after 12 calls:  ");
        for i in 0..10 {
            let v = rand() % 6;
            eprint!("{}", v);
            if i < 9 {
                eprint!(", ");
            }
        }
        eprintln!();
    }

    let mut state2: u32 = 4321;
    for _ in 0..12 {
        state2 = state2.wrapping_mul(214013).wrapping_add(2531011);
    }
    eprint!("Rust LCG%6 after 12 calls: ");
    for i in 0..10 {
        state2 = state2.wrapping_mul(214013).wrapping_add(2531011);
        let v = ((state2 >> 16) & 0x7fff) as i32 % 6;
        eprint!("{}", v);
        if i < 9 {
            eprint!(", ");
        }
    }
    eprintln!();
}

/// Call PartGraphRecursive on the coarsened graph with ncuts=1 to isolate first iteration.
#[test]
fn debug_c_recursive_coarsened_ncuts1() {
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![3, 1, 3, 0, 2, 4, 1, 5, 0, 1, 4, 5, 2, 3, 4, 3];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (2.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];

    for ncuts_val in 1..=4 {
        let mut options = vec![0 as idx_t; 40];
        unsafe {
            METIS_SetDefaultOptions(options.as_mut_ptr());
        }
        options[7] = 10; // NITER = 10
        options[8] = ncuts_val; // NCUTS = ncuts_val

        let mut edgecut: idx_t = 0;
        let mut part = vec![0 as idx_t; 6];

        let ret = unsafe {
            METIS_PartGraphRecursive(
                &mut num_vertices,
                &mut ncon,
                xadj.as_mut_ptr(),
                adjacency.as_mut_ptr(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                &mut nparts,
                ptr::null_mut(),
                ubvec.as_mut_ptr(),
                options.as_mut_ptr(),
                &mut edgecut,
                part.as_mut_ptr(),
            )
        };

        eprintln!(
            "DEBUG C METIS PartGraphRecursive (ncuts={}): ret={}, edgecut={}, part={:?}",
            ncuts_val, ret, edgecut, part
        );
    }
}

/// Call PartGraphRecursive on the 53-vertex coarsened graph from 10x10 grid.
/// This is the exact graph that init_kway_partition sees.
#[test]
fn debug_c_recursive_coarsened_53v() {
    let mut xadj: Vec<idx_t> = vec![
        0, 3, 7, 12, 16, 21, 24, 27, 33, 37, 42, 47, 51, 55, 59, 65, 69, 75, 78, 83, 87, 91, 96,
        101, 105, 108, 113, 116, 120, 126, 132, 137, 142, 146, 150, 155, 161, 167, 171, 175, 179,
        184, 189, 194, 198, 203, 206, 210, 215, 218, 221, 225, 229, 232,
    ];
    let mut adjacency: Vec<idx_t> = vec![
        7, 6, 1, 8, 0, 7, 2, 15, 1, 3, 8, 9, 10, 2, 9, 4, 16, 3, 5, 10, 11, 12, 4, 11, 13, 0, 7,
        14, 0, 6, 13, 1, 8, 14, 1, 7, 2, 21, 3, 2, 10, 15, 22, 3, 9, 4, 16, 16, 5, 4, 12, 23, 5,
        11, 16, 14, 6, 17, 7, 19, 7, 13, 18, 8, 15, 20, 2, 14, 9, 23, 4, 10, 22, 11, 12, 18, 13,
        24, 28, 14, 17, 19, 24, 28, 14, 18, 20, 29, 15, 19, 21, 29, 9, 20, 22, 25, 23, 10, 21, 25,
        16, 12, 16, 22, 26, 18, 17, 27, 31, 22, 21, 30, 26, 32, 23, 25, 34, 24, 33, 28, 35, 18, 27,
        34, 19, 29, 36, 20, 28, 35, 21, 30, 42, 25, 29, 31, 36, 43, 25, 30, 32, 37, 38, 26, 31, 37,
        45, 27, 34, 39, 35, 27, 33, 39, 28, 41, 28, 34, 40, 29, 36, 47, 29, 35, 30, 41, 42, 44, 32,
        31, 38, 48, 32, 37, 44, 40, 34, 33, 46, 50, 35, 39, 41, 46, 50, 35, 40, 36, 47, 51, 36, 41,
        42, 47, 51, 31, 42, 44, 52, 37, 43, 38, 48, 49, 33, 46, 40, 39, 45, 49, 51, 36, 41, 42, 50,
        52, 38, 44, 50, 46, 45, 47, 40, 49, 41, 52, 42, 47, 43, 48, 44, 51,
    ];
    let mut edge_weights: Vec<idx_t> = vec![
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2,
        2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1,
        1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ];
    let mut vertex_weights: Vec<idx_t> = vec![
        2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
    ];
    let mut num_vertices: idx_t = 53;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4 (ctrl->niparts)

    // Compute ubvec as init_kway_partition does:
    // kway ctrl->imbalance_tols[0] = 1.0 + 0.001 * 30 + 0.0000499 = 1.0300499
    // init_kway_partition: ubvec[i] = pow(ctrl->imbalance_tols[i], 1.0/log(nparts))
    // Then C's SetupCtrl for the inner PMETIS adds ANOTHER +0.0000499.
    // So we pass the ubvec WITHOUT the extra fudge — C will add it internally.
    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (2.0_f64).ln()) as f32;
    // Do NOT add 0.0000499 here — C's SetupCtrl adds it internally
    let mut ubvec = vec![ubvec_val];
    eprintln!(
        "DEBUG 53v ubvec_val = {} (C will add +0.0000499 internally)",
        ubvec_val
    );

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 53];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            vertex_weights.as_mut_ptr(),
            ptr::null_mut(), // vertex_sizes
            edge_weights.as_mut_ptr(),
            &mut nparts,
            ptr::null_mut(), // target_part_weights
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive on 53v (ncuts=4): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );

    // Also try ncuts=1 to see first iteration only
    options[8] = 1;
    let mut part1 = vec![0 as idx_t; 53];
    let mut edgecut1: idx_t = 0;
    let ret1 = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            vertex_weights.as_mut_ptr(),
            ptr::null_mut(),
            edge_weights.as_mut_ptr(),
            &mut nparts,
            ptr::null_mut(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut1,
            part1.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive on 53v (ncuts=1): ret={}, edgecut={}, part={:?}",
        ret1, edgecut1, part1
    );
}

/// Call PartGraphRecursive with nparts=3 on grid_3x5 to match what init_kway_partition does.
/// For grid_3x5 (15v) with kway nparts=3:
///   CoarsenTo = max(15/(40*1), 30*3) = 90, so no kway coarsening happens.
///   init_kway_partition calls METIS_PartGraphRecursive with nparts=3 on the original graph.
///   Parameters: NCUTS=4, NITER=10, seed=-1, ubvec=pow(1.0300499, 1/ln(3))
#[test]
fn debug_c_recursive_grid15_nparts3() {
    // Build 3x5 grid
    let rows = 3;
    let cols = 5;
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjacency_vec = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if r > 0 {
                adjacency_vec.push((v - cols) as i32);
            }
            if c > 0 {
                adjacency_vec.push((v - 1) as i32);
            }
            if c + 1 < cols {
                adjacency_vec.push((v + 1) as i32);
            }
            if r + 1 < rows {
                adjacency_vec.push((v + cols) as i32);
            }
            xadj[v + 1] = adjacency_vec.len() as i32;
        }
    }

    let mut num_vertices: idx_t = n as idx_t;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 3;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4

    // Compute ubvec as init_kway_partition does:
    // kway ctrl->imbalance_tols[0] = 1.0 + 0.001*30 + 0.0000499 = 1.0300499
    // ubvec = pow(1.0300499, 1/ln(3))
    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (3.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];
    eprintln!("DEBUG grid15 nparts=3: ubvec_val = {}", ubvec_val);

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; n];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency_vec.as_mut_ptr(),
            ptr::null_mut(), // vertex_weights (defaults to all 1s)
            ptr::null_mut(), // vertex_sizes
            ptr::null_mut(), // edge_weights (defaults to all 1s)
            &mut nparts,
            ptr::null_mut(), // target_part_weights (defaults to uniform 1/3)
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphRecursive (grid15, nparts=3): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );

    // Also call PartGraphKway directly for comparison
    let mut kway_options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(kway_options.as_mut_ptr());
    }
    kway_options[9] = 42; // Seed = 42

    let mut kway_edgecut: idx_t = 0;
    let mut kway_part = vec![0 as idx_t; n];
    let mut kway_nparts: idx_t = 3;

    let ret2 = unsafe {
        METIS_PartGraphKway(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency_vec.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut kway_nparts,
            ptr::null_mut(),
            ptr::null_mut(),
            kway_options.as_mut_ptr(),
            &mut kway_edgecut,
            kway_part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphKway (grid15, nparts=3, seed=42): ret={}, edgecut={}, part={:?}",
        ret2, kway_edgecut, kway_part
    );
}

/// Isolate the FIRST bisection of the nparts=3 recursive bisection.
/// Calls C's PartGraphRecursive with nparts=2 and target_part_weights=[1/3, 2/3].
/// This should produce the same partition as the first step of nparts=3 recursive bisection.
#[test]
fn debug_c_recursive_grid15_first_bisection() {
    // Build 3x5 grid
    let rows = 3;
    let cols = 5;
    let n = rows * cols;
    let mut xadj = vec![0i32; n + 1];
    let mut adjacency_vec = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            let v = r * cols + c;
            if r > 0 {
                adjacency_vec.push((v - cols) as i32);
            }
            if c > 0 {
                adjacency_vec.push((v - 1) as i32);
            }
            if c + 1 < cols {
                adjacency_vec.push((v + 1) as i32);
            }
            if r + 1 < rows {
                adjacency_vec.push((v + cols) as i32);
            }
            xadj[v + 1] = adjacency_vec.len() as i32;
        }
    }

    let mut num_vertices: idx_t = n as idx_t;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4

    // ubvec same as for nparts=3 call
    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (3.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];

    // target_part_weights = [1/3, 2/3] to match the first bisection step of nparts=3
    let mut target_part_weights = vec![1.0f32 / 3.0, 2.0f32 / 3.0];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; n];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency_vec.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            target_part_weights.as_mut_ptr(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!("DEBUG C first bisection (nparts=2, target_part_weights=[1/3, 2/3]): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part);

    // Part 0 should have ~5 vertices (1/3 of 15)
    let count0 = part.iter().filter(|&&p| p == 0).count();
    let count1 = part.iter().filter(|&&p| p == 1).count();
    eprintln!("  Part 0: {} vertices, Part 1: {} vertices", count0, count1);
}

/// Call PartGraphRecursive on 10v subgraph (right side after first bisection)
/// with the SAME seed (to match the RNG state).
/// This tests the second bisection in isolation.
#[test]
fn debug_c_recursive_grid15_second_bisection() {
    // The 10v right subgraph after splitting grid_3x5 at {4, 8, 9, 13, 14}
    // Original vertices: [0, 1, 2, 3, 5, 6, 7, 10, 11, 12]
    // Mapped to subgraph vertices 0..9
    //
    // Original 3x5 grid:
    //  0- 1- 2- 3- 4      →  0- 1- 2- 3
    //  |  |  |  |  |         |  |  |  |
    //  5- 6- 7- 8- 9      →  4- 5- 6
    //  |  |  |  |  |         |  |  |
    // 10-11-12-13-14      →  7- 8- 9
    //
    // After removing {4,8,9,13,14} and their cross-partition edges:
    // sub 0 (orig 0): [1, 5] → [sub 1, sub 4]
    // sub 1 (orig 1): [0, 2, 6] → [sub 0, sub 2, sub 5]
    // sub 2 (orig 2): [1, 3, 7] → [sub 1, sub 3, sub 6]
    // sub 3 (orig 3): [2] → [sub 2]  (orig neighbor 4 removed, orig neighbor 8 removed)
    // sub 4 (orig 5): [0, 6, 10] → [sub 0, sub 5, sub 7]
    // sub 5 (orig 6): [1, 5, 7, 11] → [sub 1, sub 4, sub 6, sub 8]
    // sub 6 (orig 7): [2, 6, 12] → [sub 2, sub 5, sub 9]  (orig neighbor 8 removed)
    // sub 7 (orig 10): [5, 11] → [sub 4, sub 8]
    // sub 8 (orig 11): [6, 10, 12] → [sub 5, sub 7, sub 9]
    // sub 9 (orig 12): [7, 11] → [sub 6, sub 8]  (orig neighbor 13 removed)
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 8, 9, 12, 16, 19, 21, 24, 26];
    let mut adjacency: Vec<idx_t> = vec![
        1, 4, // sub 0
        0, 2, 5, // sub 1
        1, 3, 6, // sub 2
        2, // sub 3
        0, 5, 7, // sub 4
        1, 4, 6, 8, // sub 5
        2, 5, 9, // sub 6
        4, 8, // sub 7
        5, 7, 9, // sub 8
        6, 8, // sub 9
    ];

    let mut num_vertices: idx_t = 10;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4
                    // Seed stays -1 (fresh start)

    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (3.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 10];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(), // target_part_weights: uniform [1/2, 1/2]
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C second bisection (10v subgraph, fresh seed): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
}

/// Verify: Call C's PartGraphRecursive on 10v subgraph with seed=1128597453
/// (which is our RNG state after 700 calls = after first bisection).
/// If this matches our second bisection, then C consumed a different number of RNG calls
/// in its first bisection.
#[test]
fn debug_c_recursive_10v_with_state_700() {
    // Same 10v subgraph as debug_c_recursive_grid15_second_bisection
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 8, 9, 12, 16, 19, 21, 24, 26];
    let mut adjacency: Vec<idx_t> = vec![
        1, 4, // sub 0
        0, 2, 5, // sub 1
        1, 3, 6, // sub 2
        2, // sub 3
        0, 5, 7, // sub 4
        1, 4, 6, 8, // sub 5
        2, 5, 9, // sub 6
        4, 8, // sub 7
        5, 7, 9, // sub 8
        6, 8, // sub 9
    ];

    let mut num_vertices: idx_t = 10;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4
    options[9] = 1128597453; // Seed = our RNG state after first bisection

    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (3.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 10];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C 10v with seed=1128597453 (state after 700 calls): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
    eprintln!("  Our second bisection result: where=[1, 0, 0, 0, 1, 0, 0, 1, 1, 1]");
    eprintln!("  If these match, C consumed different RNG calls in the first bisection.");
}

/// CRITICAL TEST: Call C on the 10v subgraph as it actually appears after
/// create_coarse_graph (edges reordered by hash table self-loop removal).
/// The 10v subgraph from SplitGraphPart gets its adjacency from the COARSENED
/// graph (where edges were reordered by CreateCoarseGraph).
#[test]
fn debug_c_recursive_10v_reordered() {
    // After outer kway coarsening + split, the 10v subgraph edges are reordered.
    // Our debug shows the graph entering init_kway has:
    // v0: [5, 1] → after removing {4,8,9,13,14}: [5] remains but 5→sub4, and 1→sub1
    // Actually, let me use the adjacency from our SplitGraphPart output.
    //
    // The graph entering init_kway_partition has the outer-coarsened adjacency.
    // After the FIRST bisection (where=[1,1,1,1,0,1,1,1,0,0,1,1,1,0,0]):
    // Right subgraph (where==1): vertices {0,1,2,3,5,6,7,10,11,12}
    //
    // In the outer-coarsened graph, vertex adjacency is:
    // v0: [5, 1]      → remove {4,8,9,13,14} neighbors: keep [5→sub4, 1→sub1] = [4, 1]
    // v1: [6, 0, 2]   → keep [6→sub5, 0→sub0, 2→sub2] = [5, 0, 2]
    // v2: [7, 1, 3]   → keep [7→sub6, 1→sub1, 3→sub3] = [6, 1, 3]
    // v3: [8, 2, 4]   → remove 4(Part0), 8(Part0): keep [2→sub2] = [2]
    // v5: [10, 0, 6]  → keep [10→sub7, 0→sub0, 6→sub5] = [7, 0, 5]
    // v6: [11, 1, 5, 7] → keep [11→sub8, 1→sub1, 5→sub4, 7→sub6] = [8, 1, 4, 6]
    // v7: [12, 2, 6]  → remove 8(Part0): keep [12→sub9, 2→sub2, 6→sub5] = [9, 2, 5]
    // v10: [5, 11]     → keep [5→sub4, 11→sub8] = [4, 8]
    // v11: [6, 10, 12] → keep [6→sub5, 10→sub7, 12→sub9] = [5, 7, 9]
    // v12: [7, 11]     → remove 13(Part0): keep [7→sub6, 11→sub8] = [6, 8]
    //
    // Wait, but the outer coarsened graph has DIFFERENT adjacency order from the original!
    // The outer-coarsened adjacency for vertex 0 is [5, 1], not [1, 5].
    // So the 10v subgraph has edges in the outer-coarsened order.
    //
    // But actually, our SplitGraphPart iterates through the original order and preserves it.
    // So the subgraph adjacency follows the parent graph's order (outer-coarsened).
    //
    // Let me use the actual adjacency from our SplitGraphPart:
    // sub 0 (orig 0): neighbors from parent [5, 1] → filter same-partition → both in → [4, 1]
    // sub 1 (orig 1): neighbors from parent [6, 0, 2] → all in → [5, 0, 2]
    // sub 2 (orig 2): neighbors from parent [7, 1, 3] → all in → [6, 1, 3]
    // sub 3 (orig 3): neighbors from parent [8, 2, 4] → 8=Part0, 4=Part0 → [2]
    // sub 4 (orig 5): neighbors from parent [10, 0, 6] → all in → [7, 0, 5]
    // sub 5 (orig 6): neighbors from parent [11, 1, 5, 7] → all in → [8, 1, 4, 6]
    // sub 6 (orig 7): neighbors from parent [12, 2, 6, 8] → 8=Part0 → [9, 2, 5]
    // Wait, the outer-coarsened adjacency for v7 is [12, 2, 6, 8]. Let me check.
    //
    // Outer coarsened adjacency (from debug):
    // [5, 1, 6, 0, 2, 7, 1, 3, 8, 2, 4, 9, 3, 10, 0, 6, 11, 1, 5, 7, 12, 2, 6, 8, 13, 3, 7, 9, 14, 4, 8, 11, 5, 12, 6, 10, 13, 7, 11, 14, 8, 12, 13, 9]
    // xadj=[0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44]
    //
    // v7 (xadj[7..8] = [20, 24]): adjacency[20..24] = [12, 2, 6, 8]
    // But 8 is in Part0! So filter: [12→sub9, 2→sub2, 6→sub5] = [9, 2, 5]
    //
    // v9 (xadj[9..10] = [28, 31]): adjacency[28..31] = [14, 4, 8]
    // But 4=Part0, 8=Part0, 14=Part0! So filter: nothing → [] (isolated vertex?!)
    // Wait, v9 is NOT in the right subgraph! partition[9]=0 (Part0). Skip.
    //
    // Let me recheck. where=[1,1,1,1,0,1,1,1,0,0,1,1,1,0,0]
    // Part0 (where=0): {4, 8, 9, 13, 14}
    // Right subgraph (where=1): {0, 1, 2, 3, 5, 6, 7, 10, 11, 12}
    //
    // v10 (xadj[10..11] = [31, 33]): adjacency[31..33] = [5, 11]
    // Both in right → [sub4, sub8] = [4, 8]
    //
    // v11 (xadj[11..12] = [33, 36]): adjacency[33..36] = [12, 6, 10]
    // Wait, that's [12→sub9, 6→sub5, 10→sub7] = [9, 5, 7]
    // Hmm, that doesn't match my earlier calculation. Let me recheck.
    //
    // Actually, I need the adjacency from the graph that enters init_kway_partition.
    // From debug:
    // adjacency=[5, 1, 6, 0, 2, 7, 1, 3, 8, 2, 4, 9, 3, 10, 0, 6, 11, 1, 5, 7, 12, 2, 6, 8, 13, 3, 7, 9, 14, 4, 8, 11, 5, 12, 6, 10, 13, 7, 11, 14, 8, 12, 13, 9]
    //
    // v0:  adjacency[0..2]   = [5, 1]
    // v1:  adjacency[2..5]   = [6, 0, 2]
    // v2:  adjacency[5..8]   = [7, 1, 3]
    // v3:  adjacency[8..11]  = [8, 2, 4]     → remove 8(P0), 4(P0) → keep [2]
    // v4:  adjacency[11..13] = [9, 3]         → v4 is Part0, skip
    // v5:  adjacency[13..16] = [10, 0, 6]
    // v6:  adjacency[16..20] = [11, 1, 5, 7]
    // v7:  adjacency[20..24] = [12, 2, 6, 8]  → remove 8(P0) → keep [12, 2, 6]
    // v8:  adjacency[24..28] = [13, 3, 7, 9]  → v8 is Part0, skip
    // v9:  adjacency[28..31] = [14, 4, 8]     → v9 is Part0, skip
    // v10: adjacency[31..33] = [5, 11]
    // v11: adjacency[33..36] = [12, 6, 10]
    //
    // Wait, let me recheck v11. xadj[11]=33, xadj[12]=36.
    // adjacency[33..36] = [12, 6, 10]
    // But the outer coarsened adjacency has v11 neighbors as... let me recount.
    // xadj=[0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44]
    //        v0 v1 v2 v3 v4 v5  v6  v7  v8  v9 v10 v11 v12 v13 v14
    // v11: xadj[11..12] = [33, 36], adjacency[33..36] = ...
    // adjacency index 33: let me count. adjacency has 44 entries.
    // [5, 1, | 6, 0, 2, | 7, 1, 3, | 8, 2, 4, | 9, 3, | 10, 0, 6, | 11, 1, 5, 7, | 12, 2, 6, 8, | 13, 3, 7, 9, | 14, 4, 8, | 5, 11, | 12, 6, 10, | 13, 7, 11, | 14, 8, 12, | 13, 9]
    //  0  1    2  3  4    5  6  7    8  9 10   11 12   13 14 15   16 17 18 19   20 21 22 23   24 25 26 27   28 29 30   31 32   33 34 35   36 37 38   39 40 41   42 43
    //
    // v10: xadj[10]=31, xadj[11]=33. adjacency[31..33] = [5, 11] → Hmm, but in the original grid, v10's neighbors are [5, 11]. Wait, the original grid:
    // v10 (row 2, col 0): up(5), right(11) → [5, 11]
    // And the outer-coarsened graph has v10: [5, 11] (both in right subgraph)
    // So sub7 (orig 10): neighbors = [4, 8] (5→sub4, 11→sub8)
    //
    // v11: xadj[11]=33, xadj[12]=36. adjacency[33..36] = [12, 6, 10]
    // But v11 in the original grid: up(6), left(10), right(12) → [6, 10, 12]
    // After outer coarsening: the order is [12, 6, 10]!
    // This is because outer CreateCoarseGraph reordered: [6, 10, 12] → self-loop removal → first becomes last → [12, 6, 10]
    //
    // Wait, that's not how it works for vertices that aren't self-loops. Let me think again.
    // For v11 in outer CreateCoarseGraph:
    // v11 is self-matched (matching[11]=11)
    // cv = coarse_map[11] = 11
    // Self-loop sentinel at position 0
    // Process edges of v=11: adjacency[33..36] of ORIGINAL graph = [6, 10, 12]
    // Edge 0: k=coarse_map[6]=6, hash(6)=6, htable[6]=cnum_edges, cadjacency[...]=6, cnum_edges++
    // Edge 1: k=coarse_map[10]=10, hash(10)=10, htable[10]=cnum_edges, cadjacency[...]=10, cnum_edges++
    // Edge 2: k=coarse_map[12]=12, hash(12)=12, htable[12]=cnum_edges, cadjacency[...]=12, cnum_edges++
    // Adjacency: [11(self), 6, 10, 12]
    // Self-loop removal: cadjacency[start] = cadjacency[cnum_edges-1] = 12
    // Result: [12, 6, 10]
    //
    // So the outer coarsening reverses the first and last edges. For v11: [6, 10, 12] → [12, 6, 10].
    //
    // Now in the subgraph, v11 becomes sub8, with neighbors mapped: 12→sub9, 6→sub5, 10→sub7
    // So sub8's adjacency is [9, 5, 7]
    //
    // But in my original manual calculation, I had sub8 (orig 11) with neighbors [5, 7, 9].
    // These are the SAME vertices but in a DIFFERENT ORDER! [9, 5, 7] vs [5, 7, 9].

    // This is important - the 10v subgraph's adjacency ORDER depends on the outer coarsening.
    //
    // Let me construct the correct 10v subgraph:
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 8, 9, 12, 16, 19, 21, 24, 26];
    let mut adjacency: Vec<idx_t> = vec![
        4, 1, // sub 0 (orig 0): parent [5, 1] → [sub4, sub1]
        5, 0, 2, // sub 1 (orig 1): parent [6, 0, 2] → [sub5, sub0, sub2]
        6, 1, 3, // sub 2 (orig 2): parent [7, 1, 3] → [sub6, sub1, sub3]
        2, // sub 3 (orig 3): parent [8, 2, 4] → remove 8,4 → [sub2]
        7, 0, 5, // sub 4 (orig 5): parent [10, 0, 6] → [sub7, sub0, sub5]
        8, 1, 4, 6, // sub 5 (orig 6): parent [11, 1, 5, 7] → [sub8, sub1, sub4, sub6]
        9, 2, 5, // sub 6 (orig 7): parent [12, 2, 6, 8] → remove 8 → [sub9, sub2, sub5]
        4, 8, // sub 7 (orig 10): parent [5, 11] → [sub4, sub8]
        9, 5, 7, // sub 8 (orig 11): parent [12, 6, 10] → [sub9, sub5, sub7]
        6, 8, // sub 9 (orig 12): parent [7, 11] → remove 13 → [sub6, sub8]
    ];
    // Wait, v12 in the parent graph: xadj[12]=36, xadj[13]=39
    // adjacency[36..39] = [13, 7, 11]
    // Remove 13(Part0): [7→sub6, 11→sub8] = [6, 8]
    // But that's [6, 8]. Hmm, let me check the self-loop reorder.
    // For v12 in outer coarsening: original neighbors are [7, 11, 13]
    // After outer CreateCoarseGraph:
    // Self-loop at position 0, then process edges [7, 11, 13]
    // Result: [self, 7, 11, 13] → remove self → [13, 7, 11]
    // So outer-coarsened v12 has [13, 7, 11]
    // In split: remove 13 → [7→sub6, 11→sub8] = [6, 8]

    let mut num_vertices: idx_t = 10;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[7] = 10; // NITER = 10
    options[8] = 4; // NCUTS = 4
    options[9] = 1128597453; // Seed = our RNG state after first bisection (700 calls)

    let kway_ubfactor: f32 = 1.0 + 0.001 * 30.0 + 0.0000499;
    let ubvec_val: f32 = (kway_ubfactor as f64).powf(1.0 / (3.0_f64).ln()) as f32;
    let mut ubvec = vec![ubvec_val];

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 10];

    let ret = unsafe {
        METIS_PartGraphRecursive(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ubvec.as_mut_ptr(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C 10v reordered with seed=1128597453: ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
    eprintln!("  Our second bisection result: where=[1, 0, 0, 0, 1, 0, 0, 1, 1, 1]");
}

/// Call PartGraphKway on the irregular_6v graph, same as compare test
#[test]
fn debug_c_kway_irregular_6v() {
    let mut xadj: Vec<idx_t> = vec![0, 2, 5, 7, 11, 14, 16];
    let mut adjacency: Vec<idx_t> = vec![1, 3, 0, 2, 3, 1, 4, 0, 1, 4, 5, 2, 3, 5, 3, 4];
    let mut num_vertices: idx_t = 6;
    let mut ncon: idx_t = 1;
    let mut nparts: idx_t = 2;

    let mut options = vec![0 as idx_t; 40];
    unsafe {
        METIS_SetDefaultOptions(options.as_mut_ptr());
    }
    options[9] = 42; // Seed = 42

    let mut edgecut: idx_t = 0;
    let mut part = vec![0 as idx_t; 6];

    let ret = unsafe {
        METIS_PartGraphKway(
            &mut num_vertices,
            &mut ncon,
            xadj.as_mut_ptr(),
            adjacency.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
            &mut nparts,
            ptr::null_mut(),
            ptr::null_mut(),
            options.as_mut_ptr(),
            &mut edgecut,
            part.as_mut_ptr(),
        )
    };

    eprintln!(
        "DEBUG C METIS PartGraphKway (seed=42, nparts=2): ret={}, edgecut={}, part={:?}",
        ret, edgecut, part
    );
}
