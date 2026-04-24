pub mod option;
mod types;
mod rng;
mod ctrl;
mod graph;
mod partition;
mod mesh;
mod util;
mod balance_kway;
mod contig;
mod minconn;

pub use types::{Idx, Real, Error, NewGraphError, NewMeshError, Result};
pub use mesh::create_graph_dual;
pub use option::Options;
use std::result::Result as StdResult;

#[derive(Debug, PartialEq)]
pub struct Graph<'a> {
    num_constraints: Idx,
    num_parts: Idx,
    xadj: &'a [Idx],
    adjacency: &'a [Idx],
    vertex_weights: Option<&'a [Idx]>,
    vertex_sizes: Option<&'a [Idx]>,
    edge_weights: Option<&'a [Idx]>,
    target_part_weights: Option<&'a [Real]>,
    ubvec: Option<&'a [Real]>,
    options: Options,
}

impl<'a> Graph<'a> {
    pub fn new(
        num_constraints: Idx,
        num_parts: Idx,
        xadj: &'a [Idx],
        adjacency: &'a [Idx],
    ) -> StdResult<Graph<'a>, NewGraphError> {
        if num_constraints < 1 {
            return Err(NewGraphError::NoConstraints);
        }
        if num_parts < 1 {
            return Err(NewGraphError::NoParts);
        }
        if xadj.len() < 1 {
            return Err(NewGraphError::msg("xadj must have at least one element"));
        }

        let num_vertices = xadj.len() - 1;
        if Idx::try_from(num_vertices).is_err() {
            return Err(NewGraphError::TooLarge);
        }
        if Idx::try_from(adjacency.len()).is_err() {
            return Err(NewGraphError::TooLarge);
        }

        // Validate xadj is non-decreasing and in bounds
        for i in 0..num_vertices {
            if xadj[i] > xadj[i + 1] {
                return Err(NewGraphError::msg("xadj must be non-decreasing"));
            }
        }
        if xadj[0] != 0 {
            return Err(NewGraphError::msg("xadj[0] must be 0"));
        }
        let last = xadj[num_vertices] as usize;
        if last != adjacency.len() {
            return Err(NewGraphError::msg(
                "xadj[n] must equal adjacency.len()",
            ));
        }

        // Validate adjacency values
        for &adj in adjacency {
            if adj < 0 || adj >= num_vertices as Idx {
                return Err(NewGraphError::msg("adjacency values out of range"));
            }
        }

        Ok(Graph {
            num_constraints,
            num_parts,
            xadj,
            adjacency,
            vertex_weights: None,
            vertex_sizes: None,
            edge_weights: None,
            target_part_weights: None,
            ubvec: None,
            options: Options::default(),
        })
    }

    pub fn new_unchecked(
        num_constraints: Idx,
        num_parts: Idx,
        xadj: &'a [Idx],
        adjacency: &'a [Idx],
    ) -> Graph<'a> {
        Graph {
            num_constraints,
            num_parts,
            xadj,
            adjacency,
            vertex_weights: None,
            vertex_sizes: None,
            edge_weights: None,
            target_part_weights: None,
            ubvec: None,
            options: Options::default(),
        }
    }

    pub fn set_vertex_weights(mut self, vertex_weights: &'a [Idx]) -> Self {
        self.vertex_weights = Some(vertex_weights);
        self
    }

    pub fn set_vertex_sizes(mut self, vertex_sizes: &'a [Idx]) -> Self {
        self.vertex_sizes = Some(vertex_sizes);
        self
    }

    pub fn set_edge_weights(mut self, edge_weights: &'a [Idx]) -> Self {
        self.edge_weights = Some(edge_weights);
        self
    }

    pub fn set_target_part_weights(mut self, target_part_weights: &'a [Real]) -> Self {
        self.target_part_weights = Some(target_part_weights);
        self
    }

    pub fn set_ubvec(mut self, ubvec: &'a [Real]) -> Self {
        self.ubvec = Some(ubvec);
        self
    }

    // --- Option builder methods ---

    pub fn seed(mut self, seed: Idx) -> Self {
        self.options.seed = Some(seed);
        self
    }

    pub fn obj_type(mut self, obj_type: option::ObjType) -> Self {
        self.options.obj_type = Some(obj_type);
        self
    }

    pub fn coarsen_type(mut self, ctype: option::CType) -> Self {
        self.options.coarsen_type = Some(ctype);
        self
    }

    pub fn init_part_type(mut self, iptype: option::IpType) -> Self {
        self.options.init_part_type = Some(iptype);
        self
    }

    pub fn refine_type(mut self, rtype: option::RType) -> Self {
        self.options.refine_type = Some(rtype);
        self
    }

    pub fn num_iter(mut self, niter: Idx) -> Self {
        self.options.num_iter = Some(niter);
        self
    }

    pub fn num_cuts(mut self, ncuts: Idx) -> Self {
        self.options.num_cuts = Some(ncuts);
        self
    }

    pub fn num_separators(mut self, nseps: Idx) -> Self {
        self.options.num_separators = Some(nseps);
        self
    }

    pub fn imbalance_factor(mut self, ufactor: Idx) -> Self {
        self.options.imbalance_factor = Some(ufactor);
        self
    }

    pub fn minimize_connectivity(mut self, minconn: bool) -> Self {
        self.options.minimize_connectivity = Some(minconn);
        self
    }

    pub fn force_contiguous(mut self, contig: bool) -> Self {
        self.options.force_contiguous = Some(contig);
        self
    }

    pub fn compress(mut self, compress: bool) -> Self {
        self.options.compress = Some(compress);
        self
    }

    pub fn cc_order(mut self, cc_order: bool) -> Self {
        self.options.cc_order = Some(cc_order);
        self
    }

    pub fn prune_factor(mut self, pfactor: Idx) -> Self {
        self.options.prune_factor = Some(pfactor);
        self
    }

    pub fn debug_level(mut self, dbglvl: option::DbgLvl) -> Self {
        self.options.debug_level = Some(dbglvl);
        self
    }

    pub fn disable_2hop(mut self, no2hop: bool) -> Self {
        self.options.disable_2hop = Some(no2hop);
        self
    }

    pub fn part_kway(self, part: &mut [Idx]) -> Result<Idx> {
        let num_vertices = self.xadj.len() - 1;
        if self.num_parts == 1 {
            for p in part.iter_mut().take(num_vertices) {
                *p = 0;
            }
            return Ok(0);
        }
        partition::kway::part_kway(self, part)
    }

    pub fn part_recursive(self, part: &mut [Idx]) -> Result<Idx> {
        let num_vertices = self.xadj.len() - 1;
        if self.num_parts == 1 {
            for p in part.iter_mut().take(num_vertices) {
                *p = 0;
            }
            return Ok(0);
        }
        partition::recursive::part_recursive(self, part)
    }
}

#[derive(Debug, PartialEq)]
pub struct Mesh<'a> {
    nn: Idx,
    num_parts: Idx,
    min_common_nodes: Idx,
    element_offsets: &'a [Idx],
    element_indices: &'a [Idx],
    vertex_weights: Option<&'a [Idx]>,
    vertex_sizes: Option<&'a [Idx]>,
    target_part_weights: Option<&'a [Real]>,
    options: Options,
}

impl<'a> Mesh<'a> {
    pub fn new(
        num_parts: Idx,
        element_offsets: &'a [Idx],
        element_indices: &'a [Idx],
    ) -> StdResult<Mesh<'a>, NewMeshError> {
        if num_parts < 1 {
            return Err(NewMeshError::NoParts);
        }
        if element_offsets.len() < 1 {
            return Err(NewMeshError::msg("element_offsets must have at least one element"));
        }

        let ne = element_offsets.len() - 1;
        if Idx::try_from(ne).is_err() {
            return Err(NewMeshError::TooLarge);
        }
        if Idx::try_from(element_indices.len()).is_err() {
            return Err(NewMeshError::TooLarge);
        }

        // Validate element_offsets
        if element_offsets[0] != 0 {
            return Err(NewMeshError::msg("element_offsets[0] must be 0"));
        }
        for i in 0..ne {
            if element_offsets[i] > element_offsets[i + 1] {
                return Err(NewMeshError::msg("element_offsets must be non-decreasing"));
            }
        }
        let last = element_offsets[ne] as usize;
        if last != element_indices.len() {
            return Err(NewMeshError::msg("element_offsets[ne] must equal element_indices.len()"));
        }

        // Compute nn (max node index + 1)
        let nn = element_indices.iter().copied().max().map(|m| m + 1).unwrap_or(0);
        if nn < 0 {
            return Err(NewMeshError::msg("negative node indices in element_indices"));
        }

        Ok(Mesh {
            nn,
            num_parts,
            min_common_nodes: 1,
            element_offsets,
            element_indices,
            vertex_weights: None,
            vertex_sizes: None,
            target_part_weights: None,
            options: Options::default(),
        })
    }

    pub fn new_unchecked(
        nn: Idx,
        num_parts: Idx,
        element_offsets: &'a [Idx],
        element_indices: &'a [Idx],
    ) -> Mesh<'a> {
        Mesh {
            nn,
            num_parts,
            min_common_nodes: 1,
            element_offsets,
            element_indices,
            vertex_weights: None,
            vertex_sizes: None,
            target_part_weights: None,
            options: Options::default(),
        }
    }

    pub fn set_vertex_weights(mut self, vertex_weights: &'a [Idx]) -> Self {
        self.vertex_weights = Some(vertex_weights);
        self
    }

    pub fn set_vertex_sizes(mut self, vertex_sizes: &'a [Idx]) -> Self {
        self.vertex_sizes = Some(vertex_sizes);
        self
    }

    pub fn set_target_part_weights(mut self, target_part_weights: &'a [Real]) -> Self {
        self.target_part_weights = Some(target_part_weights);
        self
    }

    // --- Option builder methods ---

    pub fn seed(mut self, seed: Idx) -> Self {
        self.options.seed = Some(seed);
        self
    }

    pub fn obj_type(mut self, obj_type: option::ObjType) -> Self {
        self.options.obj_type = Some(obj_type);
        self
    }

    pub fn num_cuts(mut self, ncuts: Idx) -> Self {
        self.options.num_cuts = Some(ncuts);
        self
    }

    pub fn num_iter(mut self, niter: Idx) -> Self {
        self.options.num_iter = Some(niter);
        self
    }

    pub fn imbalance_factor(mut self, ufactor: Idx) -> Self {
        self.options.imbalance_factor = Some(ufactor);
        self
    }

    pub fn part_dual(mut self, epart: &mut [Idx], npart: &mut [Idx]) -> Result<Idx> {
        let ne = self.element_offsets.len() - 1;
        if self.num_parts == 1 {
            for p in epart.iter_mut().take(ne) {
                *p = 0;
            }
            for p in npart.iter_mut().take(self.nn as usize) {
                *p = 0;
            }
            return Ok(0);
        }
        self.options.numbering = Some(option::Numbering::C);
        mesh::meshpart::part_dual(self, epart, npart)
    }
}
