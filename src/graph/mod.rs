pub mod setup;
pub mod coarsen;
pub mod contract;

use crate::types::{Idx, Real};

/// Neighbor info for k-way cut refinement.
#[derive(Clone, Copy, Default, Debug)]
pub struct CnbrInfo {
    pub pid: Idx,
    pub ed: Idx,
}

/// Per-vertex k-way cut refinement info.
#[derive(Clone, Default, Debug)]
pub struct CKRInfo {
    pub id: Idx,
    pub ed: Idx,
    pub nnbrs: Idx,
    pub inbr: i32, // index into cnbrpool, or -1
}

/// Internal graph data structure (CSR format with partitioning state).
pub struct GraphData {
    pub nvtxs: Idx,
    pub nedges: Idx, // total number of directed edges (= xadj[nvtxs])
    pub ncon: Idx,

    // CSR adjacency
    pub xadj: Vec<Idx>,
    pub adjncy: Vec<Idx>,
    pub vwgt: Vec<Idx>,
    pub vsize: Vec<Idx>,
    pub adjwgt: Vec<Idx>,

    // 2-way partition state
    pub where_: Vec<Idx>,
    pub pwgts: Vec<Idx>,
    pub bndptr: Vec<Idx>, // bndptr[v] = index in bndind, or -1
    pub bndind: Vec<Idx>, // boundary vertex list
    pub nbnd: Idx,
    pub id: Vec<Idx>, // internal degree
    pub ed: Vec<Idx>, // external degree
    pub mincut: Idx,

    // K-way refinement info
    pub ckrinfo: Vec<CKRInfo>,

    // Coarsening map: cmap[fine_vertex] = coarse_vertex
    pub cmap: Vec<Idx>,
    // match_[v] = matching partner, or v for self-match
    pub match_: Vec<Idx>,

    // Multi-level chain
    pub coarser: Option<Box<GraphData>>,
    pub finer: *mut GraphData, // raw back-pointer

    // Graph-level metadata
    pub tvwgt: Vec<Idx>,     // total vertex weight per constraint
    pub invtvwgt: Vec<Real>, // 1.0 / tvwgt
    pub label: Vec<Idx>,     // original vertex labels (for recursive bisection)
}

impl GraphData {
    pub fn new() -> Self {
        GraphData {
            nvtxs: 0,
            nedges: 0,
            ncon: 1,
            xadj: Vec::new(),
            adjncy: Vec::new(),
            vwgt: Vec::new(),
            vsize: Vec::new(),
            adjwgt: Vec::new(),
            where_: Vec::new(),
            pwgts: Vec::new(),
            bndptr: Vec::new(),
            bndind: Vec::new(),
            nbnd: 0,
            id: Vec::new(),
            ed: Vec::new(),
            mincut: 0,
            ckrinfo: Vec::new(),
            cmap: Vec::new(),
            match_: Vec::new(),
            coarser: None,
            finer: std::ptr::null_mut(),
            tvwgt: Vec::new(),
            invtvwgt: Vec::new(),
            label: Vec::new(),
        }
    }

    /// Insert vertex into boundary list.
    pub fn bnd_insert(&mut self, v: usize) {
        let pos = self.nbnd;
        self.bndind[pos as usize] = v as Idx;
        self.bndptr[v] = pos;
        self.nbnd += 1;
    }

    /// Delete vertex from boundary list.
    pub fn bnd_delete(&mut self, v: usize) {
        let pos = self.bndptr[v];
        if pos == -1 {
            return;
        }
        self.nbnd -= 1;
        let last = self.nbnd;
        if pos != last {
            let moved_vtx = self.bndind[last as usize] as usize;
            self.bndind[pos as usize] = moved_vtx as Idx;
            self.bndptr[moved_vtx] = pos;
        }
        self.bndptr[v] = -1;
    }

    /// Allocate 2-way partition memory.
    pub fn alloc_2way(&mut self) {
        let n = self.nvtxs as usize;
        let ncon = self.ncon as usize;
        self.where_ = vec![0; n];
        self.pwgts = vec![0; 2 * ncon];
        self.bndptr = vec![-1; n];
        self.bndind = vec![0; n];
        self.nbnd = 0;
        self.id = vec![0; n];
        self.ed = vec![0; n];
    }
}
