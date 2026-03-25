use std::collections::VecDeque;

/// Compares two bitstrings lexicographically.
///
/// Returns true if the first bitstring is less than the second.
fn bitstrings_less_than(x: &[i32], y: &[i32]) -> bool {
    x.iter().zip(y.iter()).any(|(a, b)| a < b)
}

/// Checks if two bitstrings are equal.
fn bitstrings_equal(x: &[i32], y: &[i32]) -> bool {
    x.iter().zip(y.iter()).all(|(a, b)| a == b)
}

/// A struct to represent and manipulate an ordered rooted tree
/// in doubly linked adjacency list representation.
///
/// This tree structure is used in Gray code algorithms for rectangulations.
#[derive(Debug, Clone)]
pub struct Tree {
    num_vertices: usize,
    root: usize,
    children: Vec<VecDeque<usize>>,
    parent: Vec<usize>,
}

impl Tree {
    /// Constructs a new tree from a bitstring vector.
    ///
    /// # Panics
    ///
    /// Panics if the bitstring length is not odd.
    ///
    /// # Arguments
    ///
    /// * `xv` - A vector of integers representing the tree encoding.
    pub fn new(xv: Vec<i32>) -> Self {
        assert!(xv.len() % 2 == 1);

        let num_vertices = (xv.len() - 1) / 2 + 1;
        let mut children = vec![VecDeque::new(); num_vertices];
        let mut parent = vec![0; num_vertices];
        let root = 0;
        let mut u = 0; // the current vertex
        let mut n = 1; // the number of vertices created so far
        for x in xv.iter().take(xv.len() - 1) {
            if *x == 1 {
                children[u].push_back(n);
                parent[n] = u;
                u = n;
                n += 1;
            } else {
                u = parent[u];
            }
        }
        assert_eq!(n, num_vertices);

        Self {
            num_vertices,
            root,
            children,
            parent,
        }
    }

    /// Computes the degree of a vertex.
    fn deg(&self, u: usize) -> usize {
        assert!(u < self.num_vertices);
        if u == self.root {
            self.children[u].len()
        } else {
            self.children[u].len() + 1
        }
    }

    /// Returns the number of children of a vertex.
    fn num_children(&self, u: usize) -> usize {
        assert!(u < self.num_vertices);
        self.children[u].len()
    }

    /// Returns the i-th child of a vertex.
    fn ith_child(&self, u: usize, i: usize) -> usize {
        assert!(u < self.num_vertices && i < self.num_children(u));
        *self.children[u].get(i).unwrap()
    }

    /// Checks if the tree is a tau preimage.
    fn is_tau_preimage(&self) -> bool {
        if self.num_vertices < 3 {
            return false;
        }
        let u = self.ith_child(self.root, 0);
        if self.num_children(u) == 0 {
            return false;
        }
        let v = self.ith_child(u, 0);
        if self.num_children(v) != 0 {
            return false;
        }
        true
    }

    /// Checks if the tree is a tau image.
    fn is_tau_image(&self) -> bool {
        if self.num_vertices < 3
            || self.num_children(self.root) < 2
            || self.num_children(self.ith_child(self.root, 0)) > 0
        {
            return false;
        }
        true
    }

    /// Applies the tau transformation to the tree.
    fn tau(&mut self) {
        assert!(self.is_tau_preimage());
        let u = self.ith_child(self.root, 0);
        let v = self.ith_child(u, 0);
        self.move_leaf(v, self.root, 0);
    }

    /// Applies the inverse tau transformation to the tree.
    fn tau_inverse(&mut self) {
        assert!(self.is_tau_image());
        let v = self.ith_child(self.root, 0);
        let u = self.ith_child(self.root, 1);
        self.move_leaf(v, u, 0);
    }

    /// Moves a leaf node to a new parent at a specified position.
    ///
    /// # Arguments
    ///
    /// * `leaf` - The leaf node to move.
    /// * `new_parent` - The new parent node.
    /// * `pos` - The position in the new parent's children list.
    pub fn move_leaf(&mut self, leaf: usize, new_parent: usize, pos: usize) {
        assert!(leaf < self.num_vertices);
        assert!(new_parent < self.num_vertices);
        assert!(pos <= self.children[new_parent].len());
        assert!(self.num_children(leaf) == 0);

        let old_parent = self.parent[leaf];
        if let Some(index) = self.children[old_parent].iter().position(|&x| x == leaf) {
            self.children[old_parent].remove(index);
        }
        self.children[new_parent].insert(pos, leaf);
        self.parent[leaf] = new_parent;
    }

    /// Rotates the tree by moving the root to its first child.
    ///
    /// This is a tree rotation operation that changes the root
    /// while maintaining the tree structure.
    pub fn rotate(&mut self) {
        assert!(self.num_vertices >= 2);
        let u = self.ith_child(self.root, 0);
        self.parent[self.root] = u;
        if let Some(first) = self.children[self.root].pop_front() {
            self.children[u].push_back(first);
        }
        self.children[u].push_back(self.root);
        self.root = u;
    }

    /// Rotates the tree until the specified vertex becomes the root.
    ///
    /// # Arguments
    ///
    /// * `u` - The target vertex to become the new root.
    pub fn rotate_to_vertex(&mut self, u: usize) {
        while self.root != u {
            self.rotate();
        }
    }

    /// Rotates children by one position (default).
    pub fn rotate_children_default(&mut self) {
        self.rotate_children(1);
    }

    /// Rotates the children of the root by k positions.
    ///
    /// # Arguments
    ///
    /// * `k` - The number of positions to rotate.
    pub fn rotate_children(&mut self, k: usize) {
        let mut queue: VecDeque<_> = self.children[self.root].iter().cloned().collect();
        for _ in 0..k {
            if let Some(front) = queue.pop_front() {
                queue.push_back(front);
            }
        }
        self.children[self.root] = queue.into_iter().collect();
    }

    /// Flips the tree structure if possible.
    ///
    /// Returns true if the flip was successful, false otherwise.
    pub fn flip_tree(&mut self) -> bool {
        if self.is_tau_preimage() && Self::is_flip_tree_tau(self) {
            self.tau();
            true
        } else if self.is_tau_image() {
            self.tau_inverse();
            if Self::is_flip_tree_tau(self) {
                return true;
            }
            self.tau(); // undo tau^{-1}
            false
        } else {
            false
        }
    }
}

impl Tree {
    /// Roots the tree canonically based on its center.
    fn root_canonically(&mut self) {
        let (c1, c2): (usize, usize) = self.compute_center(); // center vertices

        if c2 != usize::MAX {
            // centers are different
            let num_bits = 2 * (self.num_vertices - 1);
            let mut x1: Vec<i32> = vec![0; num_bits];
            let mut x2: Vec<i32> = vec![0; num_bits];
            self.rotate_to_vertex(c1);
            while self.ith_child(self.root, 0) != c2 {
                self.rotate_children_default();
            }
            self.to_bitstring(&mut x1);

            self.rotate();
            self.rotate_children(self.num_children(self.root) - 1);
            assert!(self.root == c2 && self.ith_child(self.root, 0) == c1);
            self.to_bitstring(&mut x2);

            if bitstrings_less_than(&x1[..num_bits], &x2[..num_bits]) {
                self.rotate();
                self.rotate_children(self.num_children(self.root) - 1);
                assert!(self.root == c1 && self.ith_child(self.root, 0) == c2);
            }
        } else {
            // centers are the same
            let num_bits = 2 * (self.num_vertices - 1);
            let mut x: Vec<i32> = vec![0; num_bits];
            self.rotate_to_vertex(c1);
            self.to_bitstring(&mut x);

            let mut subtree_count: Vec<i32> = vec![0; num_bits];
            let (mut c, mut depth) = (0, 0);
            for i in 0..num_bits {
                if x[i] == 1 {
                    depth += 1;
                } else {
                    // x[i] == 0
                    depth -= 1;
                }
                subtree_count[i] = c;
                if depth == 0 {
                    c += 1;
                }
            }
            assert_eq!(c, self.num_children(self.root) as i32);

            let k = self.min_string_rotation(&x, num_bits);
            self.rotate_children(subtree_count[k] as usize);
        }
    }
}

impl Tree {
    /// Computes the center vertex (or vertices) of the tree.
    ///
    /// Returns a tuple of (center1, center2) where center2 is usize::MAX
    /// if there is only one center.
    fn compute_center(&self) -> (usize, usize) {
        let mut degs = vec![0; self.num_vertices];
        let mut leaves = VecDeque::new();

        for (i, deg) in degs.iter_mut().enumerate() {
            *deg = self.deg(i);
            if *deg == 1 {
                leaves.push_back(i);
            }
        }

        let mut num_vertices_remaining = self.num_vertices;

        while num_vertices_remaining > 2 {
            for _ in 0..leaves.len() {
                let u = leaves.pop_front().unwrap();
                for &child in &self.children[u] {
                    degs[child] -= 1;
                    if degs[child] == 1 {
                        leaves.push_back(child);
                    }
                }
                if u != self.root {
                    let parent = self.parent[u];
                    degs[parent] -= 1;
                    if degs[parent] == 1 {
                        leaves.push_back(parent);
                    }
                }
            }
            num_vertices_remaining -= leaves.len();
        }

        assert!(!leaves.is_empty() && leaves.len() <= 2);

        let c1: usize;
        let c2: usize;

        if leaves.len() == 1 {
            c1 = leaves[0];
            c2 = usize::MAX;
        } else {
            c1 = leaves[0];
            c2 = leaves[1];
        }
        (c1, c2)
    }
}

impl Tree {
    /// Checks if the tree is a flip tree using tau transformation.
    fn is_flip_tree_tau(&mut self) -> bool {
        if self.is_star() {
            return false;
        }

        let r = self.root;
        let u = self.ith_child(self.root, 0);

        let num_bits = 2 * (self.num_vertices - 1);
        let mut this_bitstring = vec![0; num_bits];
        let mut canon_bitstring = vec![0; num_bits];

        let v = self.ith_child(self.root, 0);
        if (self.num_children(v) == 1) && (self.num_children(self.ith_child(v, 0)) == 0) {
            // tree has the form 1100...
            self.to_bitstring(&mut this_bitstring);
            self.root_canonically();
            let mut v = self.ith_child(self.root, 0);
            while (self.num_children(v) != 1) || (self.num_children(self.ith_child(v, 0)) != 0) {
                self.rotate();
                v = self.ith_child(self.root, 0);
            }
        } else {
            if self.has_thin_leaf() {
                return false;
            }
            let v = self.ith_child(self.root, 0);
            let c = self.count_pending_edges(v);
            if (c < self.num_children(v)) || (c < 2) || (self.is_light_dumbbell()) {
                return false;
            }
            // tree has the form 1(10)^k0... with k>=2
            self.to_bitstring(&mut this_bitstring);
            self.root_canonically();
            let mut v = self.ith_child(self.root, 0);
            let mut c = self.count_pending_edges(v);
            while (c < self.num_children(v)) || (c < 2) {
                self.rotate();
                self.rotate_children(c);
                v = self.ith_child(self.root, 0);
                c = self.count_pending_edges(v);
            }
        }

        self.to_bitstring(&mut canon_bitstring);

        self.rotate_to_vertex(r);
        while self.ith_child(self.root, 0) != u {
            self.rotate_children_default();
        }

        bitstrings_equal(&this_bitstring[..num_bits], &canon_bitstring[..num_bits])
    }

    /// Checks if the tree is a star (root connected to all other vertices).
    fn is_star(&self) -> bool {
        self.num_vertices <= 3
            || self.deg(self.root) == self.num_vertices - 1
            || self.deg(self.ith_child(self.root, 0)) == self.num_vertices - 1
    }

    /// Checks if the tree is a light dumbbell.
    fn is_light_dumbbell(&self) -> bool {
        if self.num_vertices < 5 {
            return false;
        }
        let u = self.ith_child(self.root, 0);
        let k = self.num_children(u);
        let l = self.num_children(self.root) - 1;
        !(k + l + 1 < self.num_vertices - 1 || k <= l)
    }
}

impl Tree {
    /// Checks if a vertex is a thin leaf.
    ///
    /// A thin leaf is a leaf whose parent has degree 2.
    ///
    /// # Arguments
    ///
    /// * `u` - The vertex index to check.
    pub fn is_thin_leaf(&self, u: usize) -> bool {
        if self.deg(u) > 1 {
            return false;
        }
        (u == self.root && self.deg(self.children[u][0]) == 2)
            || (u != self.root && self.deg(self.parent[u]) == 2)
    }

    /// Checks if the tree has any thin leaves.
    pub fn has_thin_leaf(&self) -> bool {
        for i in 0..self.num_vertices {
            if self.is_thin_leaf(i) {
                return true;
            }
        }
        false
    }

    /// Counts the number of pending edges from a vertex.
    ///
    /// # Arguments
    ///
    /// * `u` - The vertex index.
    pub fn count_pending_edges(&self, u: usize) -> usize {
        let mut c = 0;
        for i in 0..self.children[u].len() {
            let v = self.children[u][i];
            if self.children[v].is_empty() {
                c += 1;
            } else {
                return c;
            }
        }
        c
    }
}

impl Tree {
    /// Converts the tree to a bitstring representation.
    ///
    /// # Arguments
    ///
    /// * `x` - The output bitstring vector.
    pub fn to_bitstring(&self, x: &mut [i32]) {
        let mut pos = 0;
        self.to_bitstring_rec(x, self.root, &mut pos);
    }

    /// Recursive helper for converting tree to bitstring.
    fn to_bitstring_rec(&self, x: &mut [i32], u: usize, pos: &mut usize) {
        if self.num_children(u) == 0 {
        } else {
            for &child in &self.children[u] {
                x[*pos] = 1;
                *pos += 1;
                self.to_bitstring_rec(x, child, pos);
                x[*pos] = 0;
                *pos += 1;
            }
        }
    }

    /// Finds the minimum string rotation using Booth's algorithm.
    ///
    /// # Arguments
    ///
    /// * `x` - The input bitstring.
    /// * `length` - The length of the bitstring.
    ///
    /// Returns the starting index of the minimum rotation.
    fn min_string_rotation(&self, x: &[i32], length: usize) -> usize {
        let mut xx: Vec<i32> = x.to_vec();
        xx.extend_from_slice(x);
        let mut fail: Vec<usize> = vec![usize::MAX; 2 * length];
        let mut k = 0;
        for j in 1..2 * length {
            let xj = xx[j];
            let mut i = fail[j - k - 1];
            while i != usize::MAX && xj != xx[k + i + 1] {
                if xj < xx[k + i + 1] {
                    k = j - i - 1;
                }
                i = fail[i];
            }
            if xj != xx[k + i + 1] {
                if xj < xx[k] {
                    k = j;
                }
                fail[j - k] = usize::MAX;
            } else {
                fail[j - k] = i + 1;
            }
        }
        k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitstrings_less_than() {
        assert!(bitstrings_less_than(&[0, 1, 0], &[0, 1, 1]));
        assert!(!bitstrings_less_than(&[0, 1, 1], &[0, 1, 0]));
    }

    #[test]
    fn test_bitstrings_equal() {
        assert!(bitstrings_equal(&[0, 1, 0], &[0, 1, 0]));
        assert!(!bitstrings_equal(&[0, 1, 0], &[0, 1, 1]));
    }
}
