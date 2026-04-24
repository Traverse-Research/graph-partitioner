use crate::types::Idx;

/// MSVC-compatible LCG random number generator.
/// Replicates C stdlib rand()/srand() on MSVC:
///   state = state * 214013 + 2531011
///   output = (state >> 16) & 0x7fff
pub struct Rng {
    pub state: u32,
    pub count: u64,
}

impl Rng {
    pub fn new(seed: Idx) -> Self {
        let s = if seed == -1 { 4321 } else { seed as u32 };
        Rng { state: s, count: 0 }
    }

    pub fn call_count(&self) -> u64 {
        self.count
    }

    /// Returns a random number in [0, 32767] (matches MSVC rand()).
    pub fn rand(&mut self) -> Idx {
        self.count += 1;
        self.state = self.state.wrapping_mul(214013).wrapping_add(2531011);
        ((self.state >> 16) & 0x7fff) as Idx
    }

    /// Returns a random number in [0, max). Matches irandInRange.
    pub fn rand_in_range(&mut self, max: Idx) -> Idx {
        self.rand() % max
    }

    /// Permute an array in-place, matching GKlib irandArrayPermute exactly.
    ///
    /// If flag == true, initializes arr[offset+i] = i for i in 0..n.
    /// For n < 10: n iterations, each picks two random positions in [0,n), swaps them.
    /// For n >= 10: nshuffles iterations, each picks v,u in [0,n-3), does 4 cross-swaps.
    pub fn rand_array_permute(&mut self, n: usize, arr: &mut [Idx], offset: usize, flag: bool) {
        if flag {
            for i in 0..n {
                arr[offset + i] = i as Idx;
            }
        }

        if n < 10 {
            for _ in 0..n {
                let v = self.rand_in_range(n as Idx) as usize;
                let u = self.rand_in_range(n as Idx) as usize;
                arr.swap(offset + v, offset + u);
            }
        }
        // n >= 10 case is handled by rand_array_permute_n
    }

    /// Permute with explicit nshuffles parameter (for n >= 10).
    /// This matches the n >= 10 path of irandArrayPermute.
    pub fn rand_array_permute_with_nshuffles(
        &mut self,
        n: usize,
        arr: &mut [Idx],
        offset: usize,
        nshuffles: usize,
        flag: bool,
    ) {
        if flag {
            for i in 0..n {
                arr[offset + i] = i as Idx;
            }
        }

        if n < 10 {
            for _ in 0..n {
                let v = self.rand_in_range(n as Idx) as usize;
                let u = self.rand_in_range(n as Idx) as usize;
                arr.swap(offset + v, offset + u);
            }
        } else {
            for _ in 0..nshuffles {
                let v = self.rand_in_range((n - 3) as Idx) as usize;
                let u = self.rand_in_range((n - 3) as Idx) as usize;
                // Cross-swap pattern:
                // swap(p[v+0], p[u+2])
                // swap(p[v+1], p[u+3])
                // swap(p[v+2], p[u+0])
                // swap(p[v+3], p[u+1])
                arr.swap(offset + v, offset + u + 2);
                arr.swap(offset + v + 1, offset + u + 3);
                arr.swap(offset + v + 2, offset + u);
                arr.swap(offset + v + 3, offset + u + 1);
            }
        }
    }
}
