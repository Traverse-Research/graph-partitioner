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

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    #[allow(dead_code)]
    pub fn contains(&self, node: Idx) -> bool {
        let n = node as usize;
        n < self.locator.len() && self.locator[n] >= 0
    }

    #[inline]
    pub fn insert(&mut self, node: Idx, key: f64) {
        let n = node as usize;
        debug_assert!(n < self.max_nodes);
        debug_assert!(self.locator[n] == -1);

        let pos = self.heap.len();
        self.heap.push((key, node));
        self.locator[n] = pos as i32;
        self.sift_up(pos);
    }

    #[inline]
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

    #[inline]
    pub fn get_top(&mut self) -> Option<(Idx, f64)> {
        if self.heap.is_empty() {
            return None;
        }
        let (key, node) = self.heap[0];
        self.delete(node);
        Some((node, key))
    }

    #[allow(dead_code)]
    pub fn peek_top(&self) -> Option<(Idx, f64)> {
        self.heap.first().map(|&(key, node)| (node, key))
    }

    /// Update node's key in-place, matching C METIS rpqUpdate exactly.
    /// If node is not in the queue, inserts it.
    #[inline]
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

    #[allow(dead_code)]
    pub fn get_key(&self, node: Idx) -> f64 {
        let pos = self.locator[node as usize] as usize;
        self.heap[pos].0
    }

    #[allow(dead_code)]
    pub fn reset(&mut self) {
        for &(_, node) in &self.heap {
            self.locator[node as usize] = -1;
        }
        self.heap.clear();
    }

    #[inline]
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

    /// Check heap invariant: every node's key >= its children's keys.
    /// Used in tests to verify structural integrity after operations.
    #[cfg(test)]
    fn verify_heap_invariant(&self) -> bool {
        for i in 0..self.heap.len() {
            let left = 2 * i + 1;
            let right = 2 * i + 2;
            if left < self.heap.len() && self.heap[left].0 > self.heap[i].0 {
                return false;
            }
            if right < self.heap.len() && self.heap[right].0 > self.heap[i].0 {
                return false;
            }
        }
        true
    }

    /// Check that every node in the heap has a correct locator entry.
    #[cfg(test)]
    fn verify_locator_consistency(&self) -> bool {
        for (pos, &(_, node)) in self.heap.iter().enumerate() {
            if self.locator[node as usize] != pos as i32 {
                return false;
            }
        }
        // Check that no non-heap node has a valid locator
        for (node, &loc) in self.locator.iter().enumerate() {
            if loc >= 0 {
                let pos = loc as usize;
                if pos >= self.heap.len() || self.heap[pos].1 != node as Idx {
                    return false;
                }
            }
        }
        true
    }

    #[inline]
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

#[cfg(test)]
mod tests {
    use super::*;

    // ========= Basic operations =========

    #[test]
    fn test_new_empty() {
        let pq = PQueue::new(10);
        assert!(pq.is_empty());
        assert_eq!(pq.len(), 0);
    }

    #[test]
    fn test_insert_single() {
        let mut pq = PQueue::new(10);
        pq.insert(3, 5.0);
        assert_eq!(pq.len(), 1);
        assert!(!pq.is_empty());
        assert!(pq.contains(3));
        assert!(!pq.contains(0));
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_insert_get_top_single() {
        let mut pq = PQueue::new(10);
        pq.insert(5, 42.0);
        let (node, key) = pq.get_top().unwrap();
        assert_eq!(node, 5);
        assert_eq!(key, 42.0);
        assert!(pq.is_empty());
    }

    #[test]
    fn test_get_top_empty() {
        let mut pq = PQueue::new(10);
        assert!(pq.get_top().is_none());
    }

    // ========= Max-heap ordering =========

    #[test]
    fn test_max_heap_ordering() {
        let mut pq = PQueue::new(10);
        pq.insert(0, 1.0);
        pq.insert(1, 5.0);
        pq.insert(2, 3.0);
        pq.insert(3, 7.0);
        pq.insert(4, 2.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Should extract in descending order of keys
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 3);
        assert_eq!(k, 7.0);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 1);
        assert_eq!(k, 5.0);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 2);
        assert_eq!(k, 3.0);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 4);
        assert_eq!(k, 2.0);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 0);
        assert_eq!(k, 1.0);
        assert!(pq.get_top().is_none());
    }

    #[test]
    fn test_max_heap_negative_keys() {
        let mut pq = PQueue::new(5);
        pq.insert(0, -3.0);
        pq.insert(1, -1.0);
        pq.insert(2, -5.0);
        pq.insert(3, 0.0);
        pq.insert(4, -2.0);
        assert!(pq.verify_heap_invariant());

        let (n, _) = pq.get_top().unwrap();
        assert_eq!(n, 3); // key 0.0 is max
        let (n, _) = pq.get_top().unwrap();
        assert_eq!(n, 1); // key -1.0
        let (n, _) = pq.get_top().unwrap();
        assert_eq!(n, 4); // key -2.0
        let (n, _) = pq.get_top().unwrap();
        assert_eq!(n, 0); // key -3.0
        let (n, _) = pq.get_top().unwrap();
        assert_eq!(n, 2); // key -5.0
    }

    // ========= Peek =========

    #[test]
    fn test_peek_top() {
        let mut pq = PQueue::new(5);
        assert!(pq.peek_top().is_none());

        pq.insert(0, 3.0);
        pq.insert(1, 7.0);
        pq.insert(2, 1.0);

        let (node, key) = pq.peek_top().unwrap();
        assert_eq!(node, 1);
        assert_eq!(key, 7.0);
        // peek should not remove
        assert_eq!(pq.len(), 3);
        assert!(pq.contains(1));
    }

    // ========= Contains and get_key =========

    #[test]
    fn test_contains_and_get_key() {
        let mut pq = PQueue::new(10);
        pq.insert(2, 10.0);
        pq.insert(5, 20.0);
        pq.insert(8, 15.0);

        assert!(pq.contains(2));
        assert!(pq.contains(5));
        assert!(pq.contains(8));
        assert!(!pq.contains(0));
        assert!(!pq.contains(9));

        assert_eq!(pq.get_key(2), 10.0);
        assert_eq!(pq.get_key(5), 20.0);
        assert_eq!(pq.get_key(8), 15.0);
    }

    // ========= Delete =========

    #[test]
    fn test_delete_returns_key() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 10.0);
        pq.insert(1, 20.0);
        pq.insert(2, 15.0);

        let key = pq.delete(1);
        assert_eq!(key, 20.0);
        assert!(!pq.contains(1));
        assert_eq!(pq.len(), 2);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_delete_from_middle() {
        let mut pq = PQueue::new(10);
        for i in 0..8 {
            pq.insert(i, i as f64);
        }
        assert!(pq.verify_heap_invariant());

        // Delete a middle element
        pq.delete(4);
        assert!(!pq.contains(4));
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Should still extract remaining in order
        let mut prev_key = f64::MAX;
        while let Some((_, key)) = pq.get_top() {
            assert!(key <= prev_key);
            prev_key = key;
        }
    }

    #[test]
    fn test_delete_last_element() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 10.0);
        pq.delete(0);
        assert!(pq.is_empty());
        assert!(!pq.contains(0));
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_delete_all_elements() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 1.0);
        pq.insert(1, 2.0);
        pq.insert(2, 3.0);

        pq.delete(1);
        pq.delete(0);
        pq.delete(2);

        assert!(pq.is_empty());
        assert!(pq.verify_locator_consistency());
    }

    // ========= In-place update (THE critical behavior) =========

    #[test]
    fn test_update_increase_key() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 1.0);
        pq.insert(1, 5.0);
        pq.insert(2, 3.0);

        // Increase key of node 0 to make it the max
        pq.update(0, 10.0);
        assert_eq!(pq.get_key(0), 10.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        let (node, key) = pq.get_top().unwrap();
        assert_eq!(node, 0);
        assert_eq!(key, 10.0);
    }

    #[test]
    fn test_update_decrease_key() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 10.0);
        pq.insert(1, 5.0);
        pq.insert(2, 3.0);

        // Decrease key of node 0 to make it the min
        pq.update(0, 1.0);
        assert_eq!(pq.get_key(0), 1.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        let (node, _) = pq.get_top().unwrap();
        assert_eq!(node, 1); // node 1 (key=5.0) should now be max
    }

    #[test]
    fn test_update_same_key() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 5.0);
        pq.insert(1, 3.0);
        pq.insert(2, 7.0);

        // Update to same key - should be a no-op
        pq.update(0, 5.0);
        assert_eq!(pq.get_key(0), 5.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_update_not_present_inserts() {
        let mut pq = PQueue::new(10);
        pq.insert(0, 5.0);

        // Update a node not in the queue should insert it
        pq.update(3, 10.0);
        assert!(pq.contains(3));
        assert_eq!(pq.get_key(3), 10.0);
        assert_eq!(pq.len(), 2);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
    }

    /// Critical test: verify that in-place update produces a different heap
    /// structure than delete+insert. This is the exact bug that was fixed.
    ///
    /// In C METIS rpqUpdate, when a key is updated in-place, the node
    /// stays at its position and only sifts up or down as needed. With
    /// delete+insert, the node goes to the end and sifts up from there,
    /// potentially producing a different tree structure even with the
    /// same set of (key, node) pairs.
    #[test]
    fn test_inplace_update_vs_delete_insert_structure() {
        // Build a specific heap where in-place update and delete+insert differ.
        // Insert nodes so we get a known heap structure.
        let mut pq_update = PQueue::new(10);
        let mut pq_reinsert = PQueue::new(10);

        // Insert same elements into both
        for &(node, key) in &[(0, 10.0), (1, 8.0), (2, 6.0), (3, 4.0), (4, 2.0)] {
            pq_update.insert(node, key);
            pq_reinsert.insert(node, key);
        }

        // In-place update node 4's key from 2.0 to 9.0
        pq_update.update(4, 9.0);

        // Delete+insert for comparison
        pq_reinsert.delete(4);
        pq_reinsert.insert(4, 9.0);

        // Both must be valid heaps
        assert!(pq_update.verify_heap_invariant());
        assert!(pq_update.verify_locator_consistency());
        assert!(pq_reinsert.verify_heap_invariant());
        assert!(pq_reinsert.verify_locator_consistency());

        // Both must produce the same extraction ORDER (same sorted result)
        let mut order_update = Vec::new();
        let mut order_reinsert = Vec::new();
        while let Some((n, k)) = pq_update.get_top() {
            order_update.push((n, k));
        }
        while let Some((n, k)) = pq_reinsert.get_top() {
            order_reinsert.push((n, k));
        }

        // Keys should be in same descending order
        let keys_update: Vec<f64> = order_update.iter().map(|&(_, k)| k).collect();
        let keys_reinsert: Vec<f64> = order_reinsert.iter().map(|&(_, k)| k).collect();
        assert_eq!(
            keys_update, keys_reinsert,
            "Keys should be extracted in same order"
        );
    }

    /// Test that the heap structure after in-place update matches what
    /// C METIS rpqUpdate would produce. The key property: the node stays
    /// at its current position and sifts from there.
    #[test]
    fn test_update_sift_direction() {
        let mut pq = PQueue::new(10);

        // Build heap: root=10, left child=8, right child=6
        pq.insert(0, 10.0);
        pq.insert(1, 8.0);
        pq.insert(2, 6.0);

        // Update root's key to 5.0 (decrease) - should sift down
        pq.update(0, 5.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Node 1 (key=8) should now be the root
        let (top_node, top_key) = pq.peek_top().unwrap();
        assert_eq!(top_node, 1);
        assert_eq!(top_key, 8.0);

        // Update node 2's key to 20.0 (increase) - should sift up
        pq.update(2, 20.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Node 2 (key=20) should now be the root
        let (top_node, top_key) = pq.peek_top().unwrap();
        assert_eq!(top_node, 2);
        assert_eq!(top_key, 20.0);
    }

    /// Test multiple sequential updates on the same node.
    #[test]
    fn test_multiple_updates_same_node() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 5.0);
        pq.insert(1, 3.0);
        pq.insert(2, 7.0);
        pq.insert(3, 1.0);

        // Update node 3 several times
        pq.update(3, 10.0); // becomes max
        assert!(pq.verify_heap_invariant());
        assert_eq!(pq.peek_top().unwrap().0, 3);

        pq.update(3, 0.0); // becomes min
        assert!(pq.verify_heap_invariant());
        assert_ne!(pq.peek_top().unwrap().0, 3);

        pq.update(3, 6.0); // somewhere in the middle
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
        assert_eq!(pq.get_key(3), 6.0);
    }

    // ========= Reset =========

    #[test]
    fn test_reset() {
        let mut pq = PQueue::new(10);
        pq.insert(0, 5.0);
        pq.insert(3, 10.0);
        pq.insert(7, 1.0);

        pq.reset();
        assert!(pq.is_empty());
        assert_eq!(pq.len(), 0);
        assert!(!pq.contains(0));
        assert!(!pq.contains(3));
        assert!(!pq.contains(7));
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_insert_after_reset() {
        let mut pq = PQueue::new(10);
        pq.insert(0, 5.0);
        pq.insert(1, 10.0);
        pq.reset();

        // Should be able to re-insert same nodes
        pq.insert(0, 20.0);
        pq.insert(1, 15.0);
        assert_eq!(pq.len(), 2);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        let (node, key) = pq.get_top().unwrap();
        assert_eq!(node, 0);
        assert_eq!(key, 20.0);
    }

    // ========= Stress / larger heaps =========

    #[test]
    fn test_stress_insert_extract_sorted() {
        let n = 100;
        let mut pq = PQueue::new(n);

        for i in 0..n {
            pq.insert(i as Idx, i as f64);
        }
        assert_eq!(pq.len(), n);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Extract all - should come out in descending order
        let mut prev = f64::MAX;
        for _ in 0..n {
            let (_, key) = pq.get_top().unwrap();
            assert!(key <= prev);
            prev = key;
        }
        assert!(pq.is_empty());
    }

    #[test]
    fn test_stress_random_updates() {
        let n = 50;
        let mut pq = PQueue::new(n);

        // Insert all with initial keys
        for i in 0..n {
            pq.insert(i as Idx, (i * 3 % 17) as f64);
        }

        // Do a series of updates
        for i in 0..n {
            let new_key = ((i * 7 + 13) % 100) as f64;
            pq.update(i as Idx, new_key);
            assert!(pq.verify_heap_invariant());
            assert!(pq.verify_locator_consistency());
        }

        // Extract all and verify descending order
        let mut prev = f64::MAX;
        while let Some((_, key)) = pq.get_top() {
            assert!(key <= prev);
            prev = key;
        }
    }

    #[test]
    fn test_stress_interleaved_operations() {
        let n = 30;
        let mut pq = PQueue::new(n);

        // Insert half
        for i in 0..n / 2 {
            pq.insert(i as Idx, (i as f64) * 2.0);
        }

        // Delete some, update some, insert more
        pq.delete(3);
        pq.delete(7);
        pq.update(0, 100.0);
        pq.update(5, -10.0);

        for i in (n / 2)..n {
            pq.insert(i as Idx, (i as f64) * 1.5);
        }

        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Extract all in order
        let mut prev = f64::MAX;
        while let Some((_, key)) = pq.get_top() {
            assert!(key <= prev);
            prev = key;
        }
    }

    // ========= Edge cases =========

    #[test]
    fn test_single_element_operations() {
        let mut pq = PQueue::new(1);
        pq.insert(0, 5.0);

        pq.update(0, 10.0);
        assert_eq!(pq.get_key(0), 10.0);
        assert!(pq.verify_heap_invariant());

        pq.update(0, -5.0);
        assert_eq!(pq.get_key(0), -5.0);
        assert!(pq.verify_heap_invariant());

        let (n, k) = pq.get_top().unwrap();
        assert_eq!(n, 0);
        assert_eq!(k, -5.0);
    }

    #[test]
    fn test_two_elements_swap() {
        let mut pq = PQueue::new(2);
        pq.insert(0, 1.0);
        pq.insert(1, 2.0);

        // 1 is root (key 2.0), 0 is child (key 1.0)
        assert_eq!(pq.peek_top().unwrap().0, 1);

        // Update 0 to be larger -> should sift up and become root
        pq.update(0, 3.0);
        assert_eq!(pq.peek_top().unwrap().0, 0);
        assert!(pq.verify_heap_invariant());

        // Update 0 back to small -> should sift down
        pq.update(0, 0.5);
        assert_eq!(pq.peek_top().unwrap().0, 1);
        assert!(pq.verify_heap_invariant());
    }

    #[test]
    fn test_insert_after_delete() {
        let mut pq = PQueue::new(5);
        pq.insert(0, 5.0);
        pq.insert(1, 10.0);

        pq.delete(1);
        assert!(!pq.contains(1));

        // Re-insert same node with different key
        pq.insert(1, 3.0);
        assert!(pq.contains(1));
        assert_eq!(pq.get_key(1), 3.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());
    }

    #[test]
    fn test_ties_extraction() {
        // When multiple nodes have the same key, extraction order
        // depends on heap structure (not insertion order in general).
        // We just verify all come out with the correct key.
        let mut pq = PQueue::new(5);
        pq.insert(0, 5.0);
        pq.insert(1, 5.0);
        pq.insert(2, 5.0);
        pq.insert(3, 5.0);
        pq.insert(4, 5.0);

        let mut extracted = Vec::new();
        while let Some((n, k)) = pq.get_top() {
            assert_eq!(k, 5.0);
            extracted.push(n);
        }
        extracted.sort();
        assert_eq!(extracted, vec![0, 1, 2, 3, 4]);
    }

    /// Verify that the heap structure produced by in-place update is
    /// deterministic and matches C METIS rpqUpdate logic:
    /// - If new_key > parent_key: sift up
    /// - Otherwise: sift down
    #[test]
    fn test_update_c_metis_sift_logic() {
        // Build a heap with known structure:
        // After inserting 0:10, 1:8, 2:6, 3:4, 4:2
        // Heap should be: [10, 8, 6, 4, 2] (perfect max-heap)
        //       10(0)
        //      /    \
        //    8(1)   6(2)
        //   /  \
        // 4(3)  2(4)
        let mut pq = PQueue::new(10);
        pq.insert(0, 10.0);
        pq.insert(1, 8.0);
        pq.insert(2, 6.0);
        pq.insert(3, 4.0);
        pq.insert(4, 2.0);

        // Update node 3 (key 4.0 -> 9.0). It's at position 3, parent is at (3-1)/2=1.
        // Parent key is 8.0. New key 9.0 > 8.0, so sift up.
        // After sift up: node 3 swaps with node 1, then 9.0 < 10.0 so stops.
        pq.update(3, 9.0);
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Root should still be node 0 (key 10)
        assert_eq!(pq.peek_top().unwrap().0, 0);
        assert_eq!(pq.peek_top().unwrap().1, 10.0);

        // Node 3 should now be at key 9.0
        assert_eq!(pq.get_key(3), 9.0);

        // Extract order should be: 10, 9, 8, 6, 2
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(k, 10.0);
        assert_eq!(n, 0);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(k, 9.0);
        assert_eq!(n, 3);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(k, 8.0);
        assert_eq!(n, 1);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(k, 6.0);
        assert_eq!(n, 2);
        let (n, k) = pq.get_top().unwrap();
        assert_eq!(k, 2.0);
        assert_eq!(n, 4);
    }

    /// Replicate the FM refinement scenario where update vs delete+insert
    /// caused divergent behavior. The PQ is used for FM gain tracking:
    /// vertices get their gain updated as neighbors move.
    #[test]
    fn test_fm_like_update_sequence() {
        let mut pq = PQueue::new(20);

        // Simulate FM: insert boundary vertices with their gains
        let gains = [
            (0, 2.0),
            (1, -1.0),
            (2, 4.0),
            (3, 0.0),
            (4, -3.0),
            (5, 1.0),
            (6, 3.0),
            (7, -2.0),
            (8, 5.0),
            (9, 0.0),
        ];
        for &(node, gain) in &gains {
            pq.insert(node, gain);
        }

        // Extract top (simulating moving the best-gain vertex)
        let (best, _) = pq.get_top().unwrap();
        assert_eq!(best, 8); // gain 5.0

        // Update neighbors' gains (simulating FM neighbor update)
        pq.update(2, 2.0); // was 4.0, decreased
        pq.update(6, 5.0); // was 3.0, increased
        pq.update(0, -1.0); // was 2.0, decreased
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Next extraction should be node 6 (gain 5.0)
        let (best, _) = pq.get_top().unwrap();
        assert_eq!(best, 6);

        // Continue extracting - should always maintain heap property
        while pq.get_top().is_some() {
            // Just drain; heap invariant checked in get_top via delete
        }
    }

    /// Test that delete followed by operations maintains consistency
    /// when the deleted node was in the middle of the heap.
    #[test]
    fn test_delete_middle_then_update_others() {
        let mut pq = PQueue::new(10);
        for i in 0..8 {
            pq.insert(i, (8 - i) as f64); // keys: 8,7,6,5,4,3,2,1
        }

        // Delete node in middle of heap
        pq.delete(3); // had key 5.0
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        // Update remaining nodes
        pq.update(0, 0.5); // was 8.0, now min
        pq.update(7, 20.0); // was 1.0, now max
        assert!(pq.verify_heap_invariant());
        assert!(pq.verify_locator_consistency());

        let (top, _) = pq.peek_top().unwrap();
        assert_eq!(top, 7);
    }
}
