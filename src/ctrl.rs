use crate::types::{Idx, Real, NOPTIONS};
use crate::rng::Rng;
use crate::option::{self, Opt};
use crate::graph::CnbrInfo;

/// Internal control structure holding parsed options and runtime state.
#[allow(dead_code)]
pub struct Ctrl {
    pub optype: Idx,
    pub objtype: Idx,
    pub ctype: Idx,
    pub iptype: Idx,
    pub rtype: Idx,

    pub ncon: Idx,
    pub nparts: Idx,
    pub ncuts: Idx,
    pub nseps: Idx,
    pub niter: Idx,
    pub niparts: Idx,
    pub seed: Idx,
    pub minconn: bool,
    pub contig: bool,
    pub compress: bool,
    pub ccorder: bool,
    pub pfactor: Idx,
    pub ufactor: Idx,
    pub no2hop: bool,
    pub dbglvl: Idx,
    pub numflag: Idx,

    pub coarsen_to: Idx,

    pub ubfactors: Vec<Real>,
    pub tpwgts: Vec<Real>,
    pub maxvwgt: Vec<Idx>,

    // 2-way balance multipliers: pijbm[i*ncon+j] = invtvwgt[j] / tpwgts[i*ncon+j]
    pub pijbm: Vec<Real>,

    pub rng: Rng,

    // K-way neighbor pool
    pub cnbrpool: Vec<CnbrInfo>,
    pub cnbrpool_pos: usize,
}

impl Ctrl {
    pub fn new(options: &[Idx; NOPTIONS], ncon: Idx, nparts: Idx, is_kway: bool) -> Self {
        let get = |idx: usize, default: Idx| -> Idx {
            if options[idx] == -1 { default } else { options[idx] }
        };

        let optype: Idx = if is_kway { 0 } else { 1 };
        let objtype = get(<option::ObjType as Opt>::INDEX, 0);
        let ctype = get(<option::CType as Opt>::INDEX, 1);
        let iptype = get(<option::IpType as Opt>::INDEX, 0);
        let rtype = get(<option::RType as Opt>::INDEX, if is_kway { 1 } else { 0 });
        let ncuts = get(<option::NCuts as Opt>::INDEX, 1);
        let nseps = get(<option::NSeps as Opt>::INDEX, 1);
        let niter = get(<option::NIter as Opt>::INDEX, 10);
        let seed = get(<option::Seed as Opt>::INDEX, -1);
        let minconn = get(<option::MinConn as Opt>::INDEX, 0) != 0;
        let contig = get(<option::Contig as Opt>::INDEX, 0) != 0;
        let compress = get(<option::Compress as Opt>::INDEX, 0) != 0;
        let ccorder = get(<option::CCOrder as Opt>::INDEX, 0) != 0;
        let pfactor = get(<option::PFactor as Opt>::INDEX, 0);
        let ufactor_default = if is_kway { 30 } else { 1 };
        let ufactor = get(<option::UFactor as Opt>::INDEX, ufactor_default);
        let no2hop = get(<option::No2Hop as Opt>::INDEX, 0) != 0;
        let dbglvl = get(<option::DbgLvl as Opt>::INDEX, 0);
        let numflag = get(<option::Numbering as Opt>::INDEX, 0);

        let rng = Rng::new(seed);

        let ubfactors: Vec<Real> = (0..ncon as usize)
            .map(|_| 1.0 + 0.001 * ufactor as Real + 0.0000499)
            .collect();

        let tpwgts = vec![1.0 / nparts as Real; (nparts * ncon) as usize];
        let maxvwgt = vec![0; ncon as usize];

        Ctrl {
            optype,
            objtype,
            ctype,
            iptype,
            rtype,
            ncon,
            nparts,
            ncuts,
            nseps,
            niter,
            niparts: -1,
            seed,
            minconn,
            contig,
            compress,
            ccorder,
            pfactor,
            ufactor,
            no2hop,
            dbglvl,
            numflag,
            coarsen_to: 0,
            ubfactors,
            tpwgts,
            maxvwgt,
            pijbm: Vec::new(),
            rng,
            cnbrpool: Vec::new(),
            cnbrpool_pos: 0,
        }
    }

    pub fn set_tpwgts(&mut self, tpwgts: &[Real]) {
        self.tpwgts = tpwgts.to_vec();
    }

    pub fn set_ubvec(&mut self, ubvec: &[Real]) {
        self.ubfactors = ubvec.to_vec();
    }

    /// Allocate cnbrpool for k-way refinement.
    pub fn cnbrpool_init(&mut self, capacity: usize) {
        self.cnbrpool = vec![CnbrInfo::default(); capacity];
        self.cnbrpool_pos = 0;
    }

    pub fn cnbrpool_reset(&mut self) {
        self.cnbrpool_pos = 0;
    }

    /// Get next chunk of nbrs from pool.
    pub fn cnbrpool_get_next(&mut self, nnbrs: usize) -> i32 {
        let pos = self.cnbrpool_pos;
        self.cnbrpool_pos += nnbrs;
        // Grow pool if needed
        if self.cnbrpool_pos > self.cnbrpool.len() {
            self.cnbrpool.resize(self.cnbrpool_pos * 2, CnbrInfo::default());
        }
        pos as i32
    }

    /// Setup 2-way balance multipliers.
    pub fn setup_2way_bal_multipliers(&mut self, invtvwgt: &[Real], tpwgts2: &[Real]) {
        let ncon = self.ncon as usize;
        self.pijbm = vec![0.0; 2 * ncon];
        for i in 0..2 {
            for j in 0..ncon {
                if tpwgts2[i * ncon + j] > 0.0 {
                    self.pijbm[i * ncon + j] = invtvwgt[j] / tpwgts2[i * ncon + j];
                } else {
                    self.pijbm[i * ncon + j] = 0.0;
                }
            }
        }
    }

    /// Setup k-way balance multipliers.
    /// pijbm[i*ncon+j] = invtvwgt[j] / tpwgts[i*ncon+j]
    pub fn setup_kway_bal_multipliers(&mut self, invtvwgt: &[Real]) {
        let ncon = self.ncon as usize;
        let nparts = self.nparts as usize;
        self.pijbm = vec![0.0; nparts * ncon];
        for i in 0..nparts {
            for j in 0..ncon {
                let idx = i * ncon + j;
                if self.tpwgts[idx] > 0.0 {
                    self.pijbm[idx] = invtvwgt[j] / self.tpwgts[idx];
                } else {
                    self.pijbm[idx] = 0.0;
                }
            }
        }
    }
}
