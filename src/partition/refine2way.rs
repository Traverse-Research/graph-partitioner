use crate::types::{Idx, Real};
use crate::ctrl::Ctrl;
use crate::graph::GraphData;
use crate::partition::fm::fm_2way_cut_refine;
use crate::partition::balance::balance_2way;
use crate::partition::initpart::compute_2way_partition_params;

/// Refine2Way: refine the 2-way partition through the entire coarsening chain.
///
/// Matches the C METIS Refine2Way algorithm:
///   - `orggraph` is a raw pointer to the original (finest) graph
///   - `graph` is the coarsest graph where initial partitioning was performed
///   - Starts at the coarsest level, iterates: Balance -> FM refine -> project to finer
///   - Stops when we reach the original graph
pub fn refine_2way(ctrl: &mut Ctrl, orggraph: *mut GraphData, graph: &mut GraphData, tpwgts: &[Real]) {
    let niter = ctrl.niter;

    // Compute 2-way partition parameters at the coarsest level
    compute_2way_partition_params(ctrl, graph);

    // Use raw pointer to walk the chain from coarsest to finest
    let mut cur: *mut GraphData = graph as *mut GraphData;

    loop {
        let g = unsafe { &mut *cur };

        // Balance and refine at the current level
        balance_2way(ctrl, g, tpwgts);
        fm_2way_cut_refine(ctrl, g, tpwgts, niter);

        // If we have reached the original (finest) graph, we are done
        if cur == orggraph {
            break;
        }

        // Move to finer level
        let finer_ptr = g.finer;
        assert!(!finer_ptr.is_null(), "refine_2way: finer pointer is null before reaching orggraph");

        // Project the partition from the current (coarser) level to the finer level
        project_2way_partition(unsafe { &mut *finer_ptr });

        cur = finer_ptr;
    }
}

/// Project2WayPartition: project the 2-way partition from the coarser graph to the
/// finer graph.
///
/// Matches the C METIS Project2WayPartition algorithm:
///   - `graph` is the finer level; `graph.coarser` is the coarser level
///   - Uses `graph.cmap` to map fine vertices to coarse vertices
///   - Reuses `cmap[i]` to cache `cbndptr[cmap[i]]` for the interior-node optimization
///   - Computes id/ed for each vertex and builds the boundary list
///   - Copies mincut and pwgts from the coarser graph
///   - Frees (drops) the coarser graph
fn project_2way_partition(graph: &mut GraphData) {
    // Extract all needed data from the coarser graph before mutating this graph.
    // This avoids borrow checker conflicts (coarser is owned by graph).
    let (cwhere, cbndptr, cmincut, ncon, cpwgts) = {
        let cgraph = graph.coarser.as_ref()
            .expect("project_2way_partition: coarser graph is None");
        let ncon = cgraph.ncon as usize;
        (
            cgraph.where_.clone(),
            cgraph.bndptr.clone(),
            cgraph.mincut,
            ncon,
            cgraph.pwgts[..2 * ncon].to_vec(),
        )
    };

    // Allocate 2-way partition memory for this (finer) level
    graph.alloc_2way();

    let nvtxs = graph.nvtxs as usize;

    // Project partition: where[i] = cwhere[cmap[i]]
    // Also cache cbndptr for the interior-node optimization: cmap[i] = cbndptr[cmap[i]]
    for i in 0..nvtxs {
        let j = graph.cmap[i] as usize;
        graph.where_[i] = cwhere[j];
        graph.cmap[i] = cbndptr[j]; // Reuse cmap to store boundary info
    }

    // Compute id/ed for each vertex and build the boundary list
    let mut nbnd: Idx = 0;

    for i in 0..nvtxs {
        let istart = graph.xadj[i] as usize;
        let iend = graph.xadj[i + 1] as usize;
        let mut tid: Idx = 0;
        let mut ted: Idx = 0;

        if graph.cmap[i] == -1 {
            // Interior node: the coarser node was interior (cbndptr == -1),
            // so all neighbors must be in the same partition.
            for j in istart..iend {
                tid += graph.adjwgt[j];
            }
        } else {
            // Potentially interface node: check each neighbor's partition
            let me = graph.where_[i];
            for j in istart..iend {
                if me == graph.where_[graph.adjncy[j] as usize] {
                    tid += graph.adjwgt[j];
                } else {
                    ted += graph.adjwgt[j];
                }
            }
        }

        graph.id[i] = tid;
        graph.ed[i] = ted;

        if ted > 0 || istart == iend {
            // BNDInsert(nbnd, bndind, bndptr, i)
            graph.bndind[nbnd as usize] = i as Idx;
            graph.bndptr[i] = nbnd;
            nbnd += 1;
        }
    }

    graph.mincut = cmincut;
    graph.nbnd = nbnd;

    // Copy pwgts from coarser
    graph.pwgts[..2 * ncon].copy_from_slice(&cpwgts);

    // Free coarser graph (FreeGraph(&graph->coarser))
    graph.coarser = None;
}
