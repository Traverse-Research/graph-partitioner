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

pub use types::{Idx, Real, NOPTIONS, Error, NewGraphError, NewMeshError, Result};
use option::Opt;
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
    options: [Idx; NOPTIONS],
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
            options: [-1; NOPTIONS],
        })
    }

    pub unsafe fn new_unchecked(
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
            options: [-1; NOPTIONS],
        }
    }

    pub fn set_vertex_weights(mut self, vertex_weights: &'a [Idx]) -> Graph<'a> {
        self.vertex_weights = Some(vertex_weights);
        self
    }

    pub fn set_vertex_sizes(mut self, vertex_sizes: &'a [Idx]) -> Graph<'a> {
        self.vertex_sizes = Some(vertex_sizes);
        self
    }

    pub fn set_edge_weights(mut self, edge_weights: &'a [Idx]) -> Graph<'a> {
        self.edge_weights = Some(edge_weights);
        self
    }

    pub fn set_target_part_weights(mut self, target_part_weights: &'a [Real]) -> Graph<'a> {
        self.target_part_weights = Some(target_part_weights);
        self
    }

    pub fn set_ubvec(mut self, ubvec: &'a [Real]) -> Graph<'a> {
        self.ubvec = Some(ubvec);
        self
    }

    pub fn set_options(mut self, options: &[Idx; NOPTIONS]) -> Graph<'a> {
        self.options = *options;
        self
    }

    pub fn set_option<O: Opt>(mut self, option: O) -> Graph<'a> {
        self.options[O::INDEX] = option.value();
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
    options: [Idx; NOPTIONS],
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
            options: [-1; NOPTIONS],
        })
    }

    pub unsafe fn new_unchecked(
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
            options: [-1; NOPTIONS],
        }
    }

    pub fn set_vertex_weights(mut self, vertex_weights: &'a [Idx]) -> Mesh<'a> {
        self.vertex_weights = Some(vertex_weights);
        self
    }

    pub fn set_vertex_sizes(mut self, vertex_sizes: &'a [Idx]) -> Mesh<'a> {
        self.vertex_sizes = Some(vertex_sizes);
        self
    }

    pub fn set_target_part_weights(mut self, target_part_weights: &'a [Real]) -> Mesh<'a> {
        self.target_part_weights = Some(target_part_weights);
        self
    }

    pub fn set_options(mut self, options: &[Idx; NOPTIONS]) -> Mesh<'a> {
        self.options = *options;
        self
    }

    pub fn set_option<O: Opt>(mut self, option: O) -> Mesh<'a> {
        self.options[O::INDEX] = option.value();
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
        self.options[option::Numbering::INDEX] = option::Numbering::C.value();
        mesh::meshpart::part_dual(self, epart, npart)
    }
}
