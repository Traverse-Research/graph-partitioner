use crate::types::{Idx, Real};

#[allow(dead_code)]
/// Compute load imbalance for each constraint.
pub fn compute_load_imbalance(
    ncon: usize,
    nparts: usize,
    pwgts: &[Idx],
    tpwgts: &[Real],
    tvwgt: &[Idx],
) -> Vec<Real> {
    let mut lb = vec![0.0 as Real; ncon];

    for j in 0..ncon {
        let mut max_imbal: Real = 0.0;
        for p in 0..nparts {
            let actual = pwgts[p * ncon + j] as Real;
            let target = tpwgts[p * ncon + j] * tvwgt[j] as Real;
            if target > 0.0 {
                let imbal = actual / target;
                if imbal > max_imbal {
                    max_imbal = imbal;
                }
            }
        }
        lb[j] = max_imbal;
    }

    lb
}

#[allow(dead_code)]
/// Check if the partition is balanced within tolerance.
pub fn is_balanced(
    ncon: usize,
    nparts: usize,
    pwgts: &[Idx],
    tpwgts: &[Real],
    tvwgt: &[Idx],
    ubfactors: &[Real],
) -> bool {
    let lb = compute_load_imbalance(ncon, nparts, pwgts, tpwgts, tvwgt);
    for j in 0..ncon {
        if lb[j] > ubfactors[j] {
            return false;
        }
    }
    true
}
