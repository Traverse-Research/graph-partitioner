use crate::types::{Idx, Real};
use crate::rng::Rng;
use crate::option::Options;
use crate::graph::NeighborPartInfo;

/// Internal control structure holding parsed options and runtime state.
#[allow(dead_code)]
pub struct Control {
    pub op_type: Idx,
    pub obj_type: Idx,
    pub coarsen_type: Idx,
    pub init_part_type: Idx,
    pub refine_type: Idx,

    pub num_constraints: Idx,
    pub num_parts: Idx,
    pub num_cuts: Idx,
    pub num_separators: Idx,
    pub num_iter: Idx,
    pub num_init_parts: Idx,
    pub seed: Idx,
    pub minimize_connectivity: bool,
    pub force_contiguous: bool,
    pub compress: bool,
    pub cc_order: bool,
    pub prune_factor: Idx,
    pub imbalance_factor: Idx,
    pub disable_2hop: bool,
    pub debug_level: Idx,
    pub base_numbering: Idx,

    pub coarsen_to: Idx,

    pub imbalance_tols: Vec<Real>,
    pub target_part_weights: Vec<Real>,
    pub max_vertex_weight: Vec<Idx>,

    // 2-way balance multipliers: partition_ij_balance_multipliers[i*ncon+j] = inv_total_vertex_weight[j] / target_part_weights[i*ncon+j]
    pub partition_ij_balance_multipliers: Vec<Real>,

    pub rng: Rng,

    // K-way neighbor pool
    pub neighbor_pool: Vec<NeighborPartInfo>,
    pub neighbor_pool_pos: usize,
}

impl Control {
    pub fn new(options: &Options, ncon: Idx, nparts: Idx, is_kway: bool) -> Self {
        let optype: Idx = if is_kway { 0 } else { 1 };
        let objtype = options.obj_type.map_or(0, |v| v.value());
        let ctype = options.coarsen_type.map_or(1, |v| v.value());
        let iptype = options.init_part_type.map_or(0, |v| v.value());
        let rtype = options.refine_type.map_or(if is_kway { 1 } else { 0 }, |v| v.value());
        let ncuts = options.num_cuts.unwrap_or(1);
        let nseps = options.num_separators.unwrap_or(1);
        let niter = options.num_iter.unwrap_or(10);
        let seed = options.seed.unwrap_or(-1);
        let minconn = options.minimize_connectivity.unwrap_or(false);
        let contig = options.force_contiguous.unwrap_or(false);
        let compress = options.compress.unwrap_or(false);
        let cc_order = options.cc_order.unwrap_or(false);
        let prune_factor = options.prune_factor.unwrap_or(0);
        let ufactor_default = if is_kway { 30 } else { 1 };
        let ufactor = options.imbalance_factor.unwrap_or(ufactor_default);
        let disable_2hop = options.disable_2hop.unwrap_or(false);
        let debug_level = options.debug_level.map_or(0, |v| v.value());
        let base_numbering = options.numbering.map_or(0, |v| v.value());

        let rng = Rng::new(seed);

        let imbalance_tols: Vec<Real> = (0..ncon as usize)
            .map(|_| 1.0 + 0.001 * ufactor as Real + 0.0000499)
            .collect();

        let target_part_weights = vec![1.0 / nparts as Real; (nparts * ncon) as usize];
        let max_vertex_weights = vec![0; ncon as usize];

        Control {
            op_type: optype,
            obj_type: objtype,
            coarsen_type: ctype,
            init_part_type: iptype,
            refine_type: rtype,
            num_constraints: ncon,
            num_parts: nparts,
            num_cuts: ncuts,
            num_separators: nseps,
            num_iter: niter,
            num_init_parts: -1,
            seed,
            minimize_connectivity: minconn,
            force_contiguous: contig,
            compress,
            cc_order,
            prune_factor,
            imbalance_factor: ufactor,
            disable_2hop,
            debug_level,
            base_numbering,
            coarsen_to: 0,
            imbalance_tols,
            target_part_weights: target_part_weights,
            max_vertex_weight: max_vertex_weights,
            partition_ij_balance_multipliers: Vec::new(),
            rng,
            neighbor_pool: Vec::new(),
            neighbor_pool_pos: 0,
        }
    }

    pub fn set_target_part_weights(&mut self, target_part_weights: &[Real]) {
        self.target_part_weights = target_part_weights.to_vec();
    }

    pub fn set_ubvec(&mut self, ubvec: &[Real]) {
        self.imbalance_tols = ubvec.to_vec();
    }

    /// Allocate neighbor_pool for k-way refinement.
    pub fn init_neighbor_pool(&mut self, capacity: usize) {
        self.neighbor_pool = vec![NeighborPartInfo::default(); capacity];
        self.neighbor_pool_pos = 0;
    }

    pub fn reset_neighbor_pool(&mut self) {
        self.neighbor_pool_pos = 0;
    }

    /// Get next chunk of nbrs from pool.
    /// Clamps nnbrs to min(nparts, nnbrs) matching C METIS cnbrpoolGetNext.
    pub fn alloc_neighbor_info(&mut self, nnbrs: usize) -> i32 {
        let clamped = nnbrs.min(self.num_parts as usize);
        let pos = self.neighbor_pool_pos;
        self.neighbor_pool_pos += clamped;
        // Grow pool if needed
        if self.neighbor_pool_pos > self.neighbor_pool.len() {
            self.neighbor_pool.resize(self.neighbor_pool_pos * 2, NeighborPartInfo::default());
        }
        pos as i32
    }

    /// Setup 2-way balance multipliers.
    pub fn setup_2way_balance_multipliers(&mut self, inv_total_vertex_weight: &[Real], target_part_weights2: &[Real]) {
        let ncon = self.num_constraints as usize;
        self.partition_ij_balance_multipliers = vec![0.0; 2 * ncon];
        for i in 0..2 {
            for j in 0..ncon {
                if target_part_weights2[i * ncon + j] > 0.0 {
                    self.partition_ij_balance_multipliers[i * ncon + j] = inv_total_vertex_weight[j] / target_part_weights2[i * ncon + j];
                } else {
                    self.partition_ij_balance_multipliers[i * ncon + j] = 0.0;
                }
            }
        }
    }

    /// Setup k-way balance multipliers.
    /// partition_ij_balance_multipliers[i*ncon+j] = inv_total_vertex_weight[j] / target_part_weights[i*ncon+j]
    pub fn setup_kway_balance_multipliers(&mut self, inv_total_vertex_weight: &[Real]) {
        let ncon = self.num_constraints as usize;
        let nparts = self.num_parts as usize;
        self.partition_ij_balance_multipliers = vec![0.0; nparts * ncon];
        for i in 0..nparts {
            for j in 0..ncon {
                let idx = i * ncon + j;
                if self.target_part_weights[idx] > 0.0 {
                    self.partition_ij_balance_multipliers[idx] = inv_total_vertex_weight[j] / self.target_part_weights[idx];
                } else {
                    self.partition_ij_balance_multipliers[idx] = 0.0;
                }
            }
        }
    }
}
