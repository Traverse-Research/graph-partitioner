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
    ncon: Idx,
    nparts: Idx,
    xadj: &'a [Idx],
    adjncy: &'a [Idx],
    vwgt: Option<&'a [Idx]>,
    vsize: Option<&'a [Idx]>,
    adjwgt: Option<&'a [Idx]>,
    tpwgts: Option<&'a [Real]>,
    ubvec: Option<&'a [Real]>,
    options: [Idx; NOPTIONS],
}

impl<'a> Graph<'a> {
    pub fn new(
        ncon: Idx,
        nparts: Idx,
        xadj: &'a [Idx],
        adjncy: &'a [Idx],
    ) -> StdResult<Graph<'a>, NewGraphError> {
        if ncon < 1 {
            return Err(NewGraphError::NoConstraints);
        }
        if nparts < 1 {
            return Err(NewGraphError::NoParts);
        }
        if xadj.len() < 1 {
            return Err(NewGraphError::msg("xadj must have at least one element"));
        }

        let nvtxs = xadj.len() - 1;
        if Idx::try_from(nvtxs).is_err() {
            return Err(NewGraphError::TooLarge);
        }
        if Idx::try_from(adjncy.len()).is_err() {
            return Err(NewGraphError::TooLarge);
        }

        // Validate xadj is non-decreasing and in bounds
        for i in 0..nvtxs {
            if xadj[i] > xadj[i + 1] {
                return Err(NewGraphError::msg("xadj must be non-decreasing"));
            }
        }
        if xadj[0] != 0 {
            return Err(NewGraphError::msg("xadj[0] must be 0"));
        }
        let last = xadj[nvtxs] as usize;
        if last != adjncy.len() {
            return Err(NewGraphError::msg(
                "xadj[n] must equal adjncy.len()",
            ));
        }

        // Validate adjncy values
        for &adj in adjncy {
            if adj < 0 || adj >= nvtxs as Idx {
                return Err(NewGraphError::msg("adjncy values out of range"));
            }
        }

        Ok(Graph {
            ncon,
            nparts,
            xadj,
            adjncy,
            vwgt: None,
            vsize: None,
            adjwgt: None,
            tpwgts: None,
            ubvec: None,
            options: [-1; NOPTIONS],
        })
    }

    pub unsafe fn new_unchecked(
        ncon: Idx,
        nparts: Idx,
        xadj: &'a [Idx],
        adjncy: &'a [Idx],
    ) -> Graph<'a> {
        Graph {
            ncon,
            nparts,
            xadj,
            adjncy,
            vwgt: None,
            vsize: None,
            adjwgt: None,
            tpwgts: None,
            ubvec: None,
            options: [-1; NOPTIONS],
        }
    }

    pub fn set_vwgt(mut self, vwgt: &'a [Idx]) -> Graph<'a> {
        self.vwgt = Some(vwgt);
        self
    }

    pub fn set_vsize(mut self, vsize: &'a [Idx]) -> Graph<'a> {
        self.vsize = Some(vsize);
        self
    }

    pub fn set_adjwgt(mut self, adjwgt: &'a [Idx]) -> Graph<'a> {
        self.adjwgt = Some(adjwgt);
        self
    }

    pub fn set_tpwgts(mut self, tpwgts: &'a [Real]) -> Graph<'a> {
        self.tpwgts = Some(tpwgts);
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
        let nvtxs = self.xadj.len() - 1;
        if self.nparts == 1 {
            for p in part.iter_mut().take(nvtxs) {
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
    nparts: Idx,
    ncommon: Idx,
    eptr: &'a [Idx],
    eind: &'a [Idx],
    vwgt: Option<&'a [Idx]>,
    vsize: Option<&'a [Idx]>,
    tpwgts: Option<&'a [Real]>,
    options: [Idx; NOPTIONS],
}

impl<'a> Mesh<'a> {
    pub fn new(
        nparts: Idx,
        eptr: &'a [Idx],
        eind: &'a [Idx],
    ) -> StdResult<Mesh<'a>, NewMeshError> {
        if nparts < 1 {
            return Err(NewMeshError::NoParts);
        }
        if eptr.len() < 1 {
            return Err(NewMeshError::msg("eptr must have at least one element"));
        }

        let ne = eptr.len() - 1;
        if Idx::try_from(ne).is_err() {
            return Err(NewMeshError::TooLarge);
        }
        if Idx::try_from(eind.len()).is_err() {
            return Err(NewMeshError::TooLarge);
        }

        // Validate eptr
        if eptr[0] != 0 {
            return Err(NewMeshError::msg("eptr[0] must be 0"));
        }
        for i in 0..ne {
            if eptr[i] > eptr[i + 1] {
                return Err(NewMeshError::msg("eptr must be non-decreasing"));
            }
        }
        let last = eptr[ne] as usize;
        if last != eind.len() {
            return Err(NewMeshError::msg("eptr[ne] must equal eind.len()"));
        }

        // Compute nn (max node index + 1)
        let nn = eind.iter().copied().max().map(|m| m + 1).unwrap_or(0);
        if nn < 0 {
            return Err(NewMeshError::msg("negative node indices in eind"));
        }

        Ok(Mesh {
            nn,
            nparts,
            ncommon: 1,
            eptr,
            eind,
            vwgt: None,
            vsize: None,
            tpwgts: None,
            options: [-1; NOPTIONS],
        })
    }

    pub unsafe fn new_unchecked(
        nn: Idx,
        nparts: Idx,
        eptr: &'a [Idx],
        eind: &'a [Idx],
    ) -> Mesh<'a> {
        Mesh {
            nn,
            nparts,
            ncommon: 1,
            eptr,
            eind,
            vwgt: None,
            vsize: None,
            tpwgts: None,
            options: [-1; NOPTIONS],
        }
    }

    pub fn set_vwgt(mut self, vwgt: &'a [Idx]) -> Mesh<'a> {
        self.vwgt = Some(vwgt);
        self
    }

    pub fn set_vsize(mut self, vsize: &'a [Idx]) -> Mesh<'a> {
        self.vsize = Some(vsize);
        self
    }

    pub fn set_tpwgts(mut self, tpwgts: &'a [Real]) -> Mesh<'a> {
        self.tpwgts = Some(tpwgts);
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
        let ne = self.eptr.len() - 1;
        if self.nparts == 1 {
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
