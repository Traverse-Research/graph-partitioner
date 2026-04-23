pub mod setup;
pub mod coarsen;
pub mod contract;

use crate::types::{Idx, Real};

/// Neighbor info for k-way cut refinement.
#[derive(Clone, Copy, Default, Debug)]
pub struct NeighborPartInfo {
    pub part_id: Idx,
    pub external_degree: Idx,
}

/// Per-vertex k-way cut refinement info.
#[derive(Clone, Default, Debug)]
pub struct KwayCutInfo {
    pub internal_degree: Idx,
    pub external_degree: Idx,
    pub num_neighbors: Idx,
    pub neighbor_offset: i32, // index into neighbor_pool, or -1
}

/// Internal graph data structure (CSR format with partitioning state).
pub struct GraphData {
    pub num_vertices: Idx,
    pub num_edges: Idx, // total number of directed edges (= xadj[num_vertices])
    pub num_constraints: Idx,

    // CSR adjacency
    pub xadj: Vec<Idx>,
    pub adjacency: Vec<Idx>,
    pub vertex_weights: Vec<Idx>,
    pub vertex_sizes: Vec<Idx>,
    pub edge_weights: Vec<Idx>,

    // 2-way partition state
    pub partition: Vec<Idx>,
    pub part_weights: Vec<Idx>,
    pub boundary_map: Vec<Idx>, // boundary_map[v] = index in boundary_list, or -1
    pub boundary_list: Vec<Idx>, // boundary vertex list
    pub num_boundary: Idx,
    pub internal_degree: Vec<Idx>, // internal degree
    pub external_degree: Vec<Idx>, // external degree
    pub edge_cut: Idx,

    // K-way refinement info
    pub kway_refinement_info: Vec<KwayCutInfo>,

    // Coarsening map: coarse_map[fine_vertex] = coarse_vertex
    pub coarse_map: Vec<Idx>,
    // matching[v] = matching partner, or v for self-match
    pub matching: Vec<Idx>,

    // (Multi-level chain managed externally via Vec<GraphData>)

    // Graph-level metadata
    pub total_vertex_weight: Vec<Idx>,     // total vertex weight per constraint
    pub inv_total_vertex_weight: Vec<Real>, // 1.0 / total_vertex_weight
    pub label: Vec<Idx>,     // original vertex labels (for recursive bisection)
}

impl GraphData {
    pub fn new() -> Self {
        GraphData {
            num_vertices: 0,
            num_edges: 0,
            num_constraints: 1,
            xadj: Vec::new(),
            adjacency: Vec::new(),
            vertex_weights: Vec::new(),
            vertex_sizes: Vec::new(),
            edge_weights: Vec::new(),
            partition: Vec::new(),
            part_weights: Vec::new(),
            boundary_map: Vec::new(),
            boundary_list: Vec::new(),
            num_boundary: 0,
            internal_degree: Vec::new(),
            external_degree: Vec::new(),
            edge_cut: 0,
            kway_refinement_info: Vec::new(),
            coarse_map: Vec::new(),
            matching: Vec::new(),
            total_vertex_weight: Vec::new(),
            inv_total_vertex_weight: Vec::new(),
            label: Vec::new(),
        }
    }

    /// Insert vertex into boundary list.
    pub fn add_to_boundary(&mut self, v: usize) {
        let pos = self.num_boundary;
        self.boundary_list[pos as usize] = v as Idx;
        self.boundary_map[v] = pos;
        self.num_boundary += 1;
    }

    /// Delete vertex from boundary list.
    pub fn remove_from_boundary(&mut self, v: usize) {
        let pos = self.boundary_map[v];
        if pos == -1 {
            return;
        }
        self.num_boundary -= 1;
        let last = self.num_boundary;
        if pos != last {
            let moved_vtx = self.boundary_list[last as usize] as usize;
            self.boundary_list[pos as usize] = moved_vtx as Idx;
            self.boundary_map[moved_vtx] = pos;
        }
        self.boundary_map[v] = -1;
    }

    /// Allocate 2-way partition memory.
    pub fn alloc_2way(&mut self) {
        let n = self.num_vertices as usize;
        let ncon = self.num_constraints as usize;
        self.partition = vec![0; n];
        self.part_weights = vec![0; 2 * ncon];
        self.boundary_map = vec![-1; n];
        self.boundary_list = vec![0; n];
        self.num_boundary = 0;
        self.internal_degree = vec![0; n];
        self.external_degree = vec![0; n];
    }
}
