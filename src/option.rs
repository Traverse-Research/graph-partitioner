use crate::types::Idx;

mod private {
    pub trait Sealed {}
}

pub trait Opt: private::Sealed {
    const INDEX: usize;
    fn value(self) -> Idx;
}

// --- Enum options ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PType {
    Rb,
    Kway,
}

impl private::Sealed for PType {}
impl Opt for PType {
    const INDEX: usize = 0;
    fn value(self) -> Idx {
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

impl private::Sealed for ObjType {}
impl Opt for ObjType {
    const INDEX: usize = 1;
    fn value(self) -> Idx {
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

impl private::Sealed for CType {}
impl Opt for CType {
    const INDEX: usize = 2;
    fn value(self) -> Idx {
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

impl private::Sealed for IpType {}
impl Opt for IpType {
    const INDEX: usize = 3;
    fn value(self) -> Idx {
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

impl private::Sealed for RType {}
impl Opt for RType {
    const INDEX: usize = 4;
    fn value(self) -> Idx {
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

impl private::Sealed for DbgLvl {}
impl Opt for DbgLvl {
    const INDEX: usize = 5;
    fn value(self) -> Idx {
        let mut v: Idx = 0;
        if self.info { v |= 1; }
        if self.time { v |= 2; }
        if self.coarsen { v |= 4; }
        if self.refine { v |= 8; }
        if self.ipart { v |= 16; }
        if self.move_info { v |= 32; }
        if self.sep_info { v |= 64; }
        if self.conn_info { v |= 128; }
        if self.contig_info { v |= 256; }
        v
    }
}

// --- Integer-valued tuple struct options ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NIter(pub Idx);
impl private::Sealed for NIter {}
impl Opt for NIter {
    const INDEX: usize = 7;
    fn value(self) -> Idx { self.0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NCuts(pub Idx);
impl private::Sealed for NCuts {}
impl Opt for NCuts {
    const INDEX: usize = 8;
    fn value(self) -> Idx { self.0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed(pub Idx);
impl private::Sealed for Seed {}
impl Opt for Seed {
    const INDEX: usize = 9;
    fn value(self) -> Idx { self.0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinConn(pub bool);
impl private::Sealed for MinConn {}
impl Opt for MinConn {
    const INDEX: usize = 11;
    fn value(self) -> Idx { self.0 as Idx }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Contig(pub bool);
impl private::Sealed for Contig {}
impl Opt for Contig {
    const INDEX: usize = 12;
    fn value(self) -> Idx { self.0 as Idx }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Compress(pub bool);
impl private::Sealed for Compress {}
impl Opt for Compress {
    const INDEX: usize = 13;
    fn value(self) -> Idx { self.0 as Idx }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CCOrder(pub bool);
impl private::Sealed for CCOrder {}
impl Opt for CCOrder {
    const INDEX: usize = 14;
    fn value(self) -> Idx { self.0 as Idx }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PFactor(pub Idx);
impl private::Sealed for PFactor {}
impl Opt for PFactor {
    const INDEX: usize = 15;
    fn value(self) -> Idx { self.0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NSeps(pub Idx);
impl private::Sealed for NSeps {}
impl Opt for NSeps {
    const INDEX: usize = 16;
    fn value(self) -> Idx { self.0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UFactor(pub Idx);
impl private::Sealed for UFactor {}
impl Opt for UFactor {
    const INDEX: usize = 17;
    fn value(self) -> Idx { self.0 }
}

// --- Internal-only Numbering option ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum Numbering {
    C,
    Fortran,
}

impl private::Sealed for Numbering {}
impl Opt for Numbering {
    const INDEX: usize = 18;
    fn value(self) -> Idx {
        match self {
            Numbering::C => 0,
            Numbering::Fortran => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct No2Hop(pub bool);
impl private::Sealed for No2Hop {}
impl Opt for No2Hop {
    const INDEX: usize = 20;
    fn value(self) -> Idx { self.0 as Idx }
}
