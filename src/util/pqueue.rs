use crate::types::Idx;

/// A priority queue (max-heap) keyed by f64, with vertex IDs.
/// Supports insert, delete, update, and get-max operations.
pub struct PQueue {
    heap: Vec<(f64, Idx)>,
    locator: Vec<i32>, // Maps vertex ID -> index in heap (-1 if not present)
    max_nodes: usize,
}

impl PQueue {
    pub fn new(max_nodes: usize) -> Self {
        PQueue {
            heap: Vec::with_capacity(max_nodes),
            locator: vec![-1; max_nodes],
            max_nodes,
        }
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    pub fn contains(&self, node: Idx) -> bool {
        let n = node as usize;
        n < self.locator.len() && self.locator[n] >= 0
    }

    pub fn insert(&mut self, node: Idx, key: f64) {
        let n = node as usize;
        debug_assert!(n < self.max_nodes);
        debug_assert!(self.locator[n] == -1);

        let pos = self.heap.len();
        self.heap.push((key, node));
        self.locator[n] = pos as i32;
        self.sift_up(pos);
    }

    pub fn delete(&mut self, node: Idx) -> f64 {
        let n = node as usize;
        let pos = self.locator[n] as usize;
        let key = self.heap[pos].0;
        self.locator[n] = -1;

        let last = self.heap.len() - 1;
        if pos != last {
            self.heap[pos] = self.heap[last];
            self.locator[self.heap[pos].1 as usize] = pos as i32;
            self.heap.pop();
            if pos < self.heap.len() {
                self.sift_up(pos);
                self.sift_down(pos);
            }
        } else {
            self.heap.pop();
        }

        key
    }

    pub fn get_top(&mut self) -> Option<(Idx, f64)> {
        if self.heap.is_empty() {
            return None;
        }
        let (key, node) = self.heap[0];
        self.delete(node);
        Some((node, key))
    }

    pub fn peek_top(&self) -> Option<(Idx, f64)> {
        self.heap.first().map(|&(key, node)| (node, key))
    }

    /// Update node's key in-place, matching C METIS rpqUpdate exactly.
    /// If node is not in the queue, inserts it.
    pub fn update(&mut self, node: Idx, new_key: f64) {
        let n = node as usize;
        if n >= self.locator.len() || self.locator[n] < 0 {
            self.insert(node, new_key);
            return;
        }

        let pos = self.locator[n] as usize;
        self.heap[pos].0 = new_key;

        // C METIS: if (i > 0 && newkey > queue->heap[(i-1)/2].key) sift_up else sift_down
        if pos > 0 && new_key > self.heap[(pos - 1) / 2].0 {
            self.sift_up(pos);
        } else {
            self.sift_down(pos);
        }
    }

    pub fn get_key(&self, node: Idx) -> f64 {
        let pos = self.locator[node as usize] as usize;
        self.heap[pos].0
    }

    pub fn reset(&mut self) {
        for &(_, node) in &self.heap {
            self.locator[node as usize] = -1;
        }
        self.heap.clear();
    }

    fn sift_up(&mut self, mut pos: usize) {
        while pos > 0 {
            let parent = (pos - 1) / 2;
            if self.heap[pos].0 > self.heap[parent].0 {
                self.locator[self.heap[pos].1 as usize] = parent as i32;
                self.locator[self.heap[parent].1 as usize] = pos as i32;
                self.heap.swap(pos, parent);
                pos = parent;
            } else {
                break;
            }
        }
    }

    fn sift_down(&mut self, mut pos: usize) {
        let n = self.heap.len();
        loop {
            let left = 2 * pos + 1;
            let right = 2 * pos + 2;
            let mut largest = pos;

            if left < n && self.heap[left].0 > self.heap[largest].0 {
                largest = left;
            }
            if right < n && self.heap[right].0 > self.heap[largest].0 {
                largest = right;
            }

            if largest != pos {
                self.locator[self.heap[pos].1 as usize] = largest as i32;
                self.locator[self.heap[largest].1 as usize] = pos as i32;
                self.heap.swap(pos, largest);
                pos = largest;
            } else {
                break;
            }
        }
    }
}
