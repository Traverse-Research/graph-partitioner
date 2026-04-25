use crate::types::Idx;

// --- Enum options ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PType {
    Rb,
    Kway,
}

impl PType {
    #[expect(dead_code)]
    pub(crate) fn value(self) -> Idx {
        match self {
            PType::Rb => 0,
            PType::Kway => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjType {
    Cut,
    Vol,
}

impl ObjType {
    pub(crate) fn value(self) -> Idx {
        match self {
            ObjType::Cut => 0,
            ObjType::Vol => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CType {
    Rm,
    Shem,
}

impl CType {
    pub(crate) fn value(self) -> Idx {
        match self {
            CType::Rm => 0,
            CType::Shem => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IpType {
    Grow,
    Random,
    Edge,
    Node,
}

impl IpType {
    pub(crate) fn value(self) -> Idx {
        match self {
            IpType::Grow => 0,
            IpType::Random => 1,
            IpType::Edge => 2,
            IpType::Node => 3,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RType {
    Fm,
    Greedy,
    Sep2Sided,
    Sep1Sided,
}

impl RType {
    pub(crate) fn value(self) -> Idx {
        match self {
            RType::Fm => 0,
            RType::Greedy => 1,
            RType::Sep2Sided => 2,
            RType::Sep1Sided => 3,
        }
    }
}

// --- Bitfield option ---

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DbgLvl {
    pub info: bool,
    pub time: bool,
    pub coarsen: bool,
    pub refine: bool,
    pub ipart: bool,
    pub move_info: bool,
    pub sep_info: bool,
    pub conn_info: bool,
    pub contig_info: bool,
}

impl DbgLvl {
    pub(crate) fn value(self) -> Idx {
        let mut v: Idx = 0;
        if self.info {
            v |= 1;
        }
        if self.time {
            v |= 2;
        }
        if self.coarsen {
            v |= 4;
        }
        if self.refine {
            v |= 8;
        }
        if self.ipart {
            v |= 16;
        }
        if self.move_info {
            v |= 32;
        }
        if self.sep_info {
            v |= 64;
        }
        if self.conn_info {
            v |= 128;
        }
        if self.contig_info {
            v |= 256;
        }
        v
    }
}

// --- Internal-only Numbering option ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[expect(dead_code)]
pub(crate) enum Numbering {
    C,
    Fortran,
}

impl Numbering {
    pub(crate) fn value(self) -> Idx {
        match self {
            Numbering::C => 0,
            Numbering::Fortran => 1,
        }
    }
}

// --- Options struct ---

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Options {
    pub obj_type: Option<ObjType>,
    pub coarsen_type: Option<CType>,
    pub init_part_type: Option<IpType>,
    pub refine_type: Option<RType>,
    pub debug_level: Option<DbgLvl>,
    pub num_iter: Option<Idx>,
    pub num_cuts: Option<Idx>,
    pub seed: Option<Idx>,
    pub minimize_connectivity: Option<bool>,
    pub force_contiguous: Option<bool>,
    pub compress: Option<bool>,
    pub cc_order: Option<bool>,
    pub prune_factor: Option<Idx>,
    pub num_separators: Option<Idx>,
    pub imbalance_factor: Option<Idx>,
    pub disable_2hop: Option<bool>,
    pub(crate) numbering: Option<Numbering>,
}
