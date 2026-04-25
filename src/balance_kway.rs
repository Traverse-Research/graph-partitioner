use crate::types::{Idx, Real};

/// Compute load imbalance for each constraint.
pub fn compute_load_imbalance(
    ncon: usize,
    nparts: usize,
    part_weights: &[Idx],
    target_part_weights: &[Real],
    total_vertex_weight: &[Idx],
) -> Vec<Real> {
    let mut lb = vec![0.0 as Real; ncon];

    for j in 0..ncon {
        let mut max_imbal: Real = 0.0;
        for p in 0..nparts {
            let actual = part_weights[p * ncon + j] as Real;
            let target = target_part_weights[p * ncon + j] * total_vertex_weight[j] as Real;
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

#[expect(dead_code)]
/// Check if the partition is balanced within tolerance.
pub fn is_balanced(
    ncon: usize,
    nparts: usize,
    part_weights: &[Idx],
    target_part_weights: &[Real],
    total_vertex_weight: &[Idx],
    imbalance_tols: &[Real],
) -> bool {
    let lb = compute_load_imbalance(
        ncon,
        nparts,
        part_weights,
        target_part_weights,
        total_vertex_weight,
    );
    for j in 0..ncon {
        if lb[j] > imbalance_tols[j] {
            return false;
        }
    }
    true
}
