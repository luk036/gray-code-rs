use std::collections::{LinkedList, VecDeque};
use std::fmt;
use std::cmp;

#[derive(Debug, Clone, PartialEq)]
struct Vertex {
    bits: Vec<i32>,
}

impl Vertex {
    fn new(x: Vec<i32>) -> Self {
        assert!(x.len() % 2 == 1);
        assert!(x.len() >= 3);
        Vertex { bits: x }
    }

    fn len(&self) -> usize {
        self.bits.len()
    }

    fn rev_inv(&mut self, left: Option<usize>, right: Option<usize>) {
        let left = left.unwrap_or(0);
        let right = right.unwrap_or(self.bits.len() - 2);
        
        for i in left..=right {
            self.bits[i] = 1 - self.bits[i];
        }
        self.bits[left..=right].reverse();
    }

    fn first_touchdown(&self, a: usize) -> Option<usize> {
        let mut height = 0;
        for i in a..self.bits.len() - 1 {
            height += 2 * self.bits[i] - 1;
            if height == 0 {
                return Some(i);
            }
        }
        None
    }

    fn first_dive(&self) -> Option<usize> {
        let mut height = 0;
        for i in 0..self.bits.len() - 1 {
            height += 2 * self.bits[i] - 1;
            if height == -1 {
                return Some(i);
            }
        }
        None
    }

    fn steps_height(&self) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let mut usteps_neg = Vec::new();
        let mut usteps_pos = Vec::new();
        let mut dsteps_neg = Vec::new();
        let mut dsteps_pos = Vec::new();
        
        let mut height = 0;
        let mut min_height = 0;
        let mut max_height = 0;
        
        for i in 0..self.bits.len() - 1 {
            if self.bits[i] == 0 && height <= 0 {
                if height == min_height {
                    usteps_neg.push(Vec::new());
                    dsteps_neg.push(Vec::new());
                }
                dsteps_neg[-height as usize].push(i);
            }
            
            if self.bits[i] == 1 && height >= 0 {
                if height == max_height {
                    usteps_pos.push(Vec::new());
                    dsteps_pos.push(Vec::new());
                }
                usteps_pos[height as usize].push(i);
            }
            
            height += 2 * self.bits[i] - 1;
            min_height = cmp::min(height, min_height);
            max_height = cmp::max(height, max_height);
            
            if self.bits[i] == 0 && height >= 0 {
                dsteps_pos[height as usize].push(i);
                assert_eq!(dsteps_pos[height as usize].len(), usteps_pos[height as usize].len());
            }
            
            if self.bits[i] == 1 && height <= 0 {
                usteps_neg[-height as usize].push(i);
                assert_eq!(usteps_neg[-height as usize].len(), dsteps_neg[-height as usize].len());
            }
        }
        
        assert_eq!(usteps_neg.len(), dsteps_neg.len());
        (usteps_neg, usteps_pos, dsteps_neg, dsteps_pos)
    }

    fn count_flaws(&self) -> usize {
        let mut c = 0;
        let mut height = 0;
        for i in 0..self.bits.len() - 1 {
            if height <= 0 && self.bits[i] == 0 {
                c += 1;
            }
            height += 2 * self.bits[i] - 1;
        }
        c
    }

    fn count_ones(&self) -> usize {
        self.bits[..self.bits.len() - 1].iter().filter(|&&x| x == 1).count()
    }

    fn is_first_vertex(&self) -> bool {
        self.count_flaws() == 0 && self.count_ones() == self.bits.len() / 2
    }

    fn is_last_vertex(&self) -> bool {
        self.count_flaws() == 1 && self.count_ones() == self.bits.len() / 2
    }

    fn to_first_vertex(&mut self) -> usize {
        if self.is_last_vertex() {
            let b = self.first_dive().unwrap();
            self.bits.copy_within(..b, 1);
            self.bits[0] = 1;
            self.bits[b + 1] = 0;
            2 * b + 2
        } else {
            let (usteps_neg, usteps_pos, dsteps_neg, dsteps_pos) = self.steps_height();
            let min_zero = usteps_neg.is_empty();
            let unique_min = if min_zero {
                usteps_pos[0].len() == 1
            } else {
                usteps_neg.last().unwrap().len() == 1
            };
            let middle_level = 2 * self.count_ones() + 1 == self.bits.len();
            
            let to = if (!unique_min && middle_level) || (unique_min && !middle_level) {
                if min_zero {
                    usteps_pos[0][0]
                } else {
                    usteps_neg.last().unwrap()[0]
                } - 1
            } else {
                if min_zero {
                    usteps_pos[0].last().unwrap()
                } else {
                    usteps_neg.last().unwrap().last().unwrap()
                } - 1
            };
            
            self.bits.copy_within(..=to, 1);
            self.bits[0] = 1;
            
            for d in 0..dsteps_neg.len() - if unique_min && middle_level { 1 } else { 0 } {
                self.bits[dsteps_neg[d][0] + 1] = 1;
            }
            
            for d in 0..usteps_neg.len() - if unique_min && !middle_level { 1 } else { 0 } {
                self.bits[usteps_neg[d].last().unwrap()] = 0;
            }
            
            if !middle_level {
                for d in if min_zero && unique_min { 1 } else { 0 }..=1 {
                    self.bits[usteps_pos[d].last().unwrap()] = 0;
                }
            }
            
            2 * (to + 1) + if middle_level { 0 } else { 1 }
        }
    }

    fn to_last_vertex(&mut self) -> usize {
        let mut d = 0;
        if !self.is_first_vertex() {
            d = -(self.to_first_vertex() as isize) as usize;
        }
        assert!(self.is_first_vertex());

        let b = self.first_touchdown(0).unwrap();
        self.bits.copy_within(1..b, 0);
        self.bits[b - 1] = 0;
        self.bits[b] = 1;
        d + 2 * (b - 1) + 2
    }

    fn compute_flip_seq_0(&mut self, flip: bool) -> Vec<usize> {
        assert!(self.is_first_vertex());
        
        if !flip {
            let b = self.first_touchdown(0).unwrap();
            let length = 2 * (b - 1) + 2;
            let mut seq = vec![0; length];
            let mut next_step = vec![0; b + 1];
            self.aux_pointers(0, b, &mut next_step);
            
            let mut idx = 0;
            seq[idx] = b;
            idx += 1;
            seq[idx] = 0;
            idx += 1;
            self.compute_flip_seq_0_rec(&mut seq, &mut idx, 1, b - 1, &next_step);
            seq
        } else {
            assert!(self.bits[0] == 1);
            
            if self.bits[1] == 1 {
                assert!(self.bits[2] == 0);
                vec![2, 0]
            } else {
                self.bits[1] = 1;
                self.bits[2] = 0;
                
                let b = self.first_touchdown(0).unwrap();
                let length = 2 * (b - 1) + 2;
                let mut seq = vec![0; length];
                let mut next_step = vec![0; b + 1];
                self.aux_pointers(0, b, &mut next_step);
                
                let mut idx = 0;
                seq[idx] = b;
                idx += 1;
                seq[idx] = 0;
                idx += 1;
                self.compute_flip_seq_0_rec(&mut seq, &mut idx, 1, b - 1, &next_step);
                
                self.bits[1] = 0;
                self.bits[2] = 1;
                
                assert!(seq[0] == b && seq[1] == 0 && seq[2] == 2 && seq[3] == 1 && seq[4] == 0 && seq[5] == 2);
                seq[0] = b;
                seq[1] = 0;
                seq[2] = 1;
                seq[3] = 2;
                seq[4] = 0;
                seq[5] = 1;
                seq
            }
        }
    }

    fn compute_flip_seq_0_rec(
        &self,
        seq: &mut [usize],
        idx: &mut usize,
        left: usize,
        right: usize,
        next_step: &[usize],
    ) {
        let length = right - left + 1;
        if length <= 0 {
            return;
        }
        
        assert!(self.bits[left] == 1 && self.bits[right] == 0 && length % 2 == 0);
        
        let m = next_step[left];
        assert!(m <= right && self.bits[m] == 0);
        seq[*idx] = m;
        *idx += 1;
        seq[*idx] = left;
        *idx += 1;
        self.compute_flip_seq_0_rec(seq, idx, left + 1, m - 1, next_step);
        seq[*idx] = left - 1;
        *idx += 1;
        seq[*idx] = m;
        *idx += 1;
        self.compute_flip_seq_0_rec(seq, idx, m + 1, right, next_step);
    }

    fn compute_flip_seq_1(&self) -> Vec<usize> {
        assert!(self.is_last_vertex());
        let b = self.first_dive().unwrap();
        let length = 2 * ((self.bits.len() - 2) - (b + 2) + 1) + 2;
        let mut seq = vec![0; length];
        let mut next_step = vec![0; self.bits.len() - 1];
        self.aux_pointers(b + 2, self.bits.len() - 2, &mut next_step);
        
        let mut idx = 0;
        seq[idx] = b + 1;
        idx += 1;
        self.compute_flip_seq_1_rec(&mut seq, &mut idx, b + 2, self.bits.len() - 2, &next_step);
        seq[idx] = b;
        idx += 1;
        seq
    }

    fn compute_flip_seq_1_rec(
        &self,
        seq: &mut [usize],
        idx: &mut usize,
        left: usize,
        right: usize,
        next_step: &[usize],
    ) {
        let length = right - left + 1;
        if length <= 0 {
            return;
        }
        
        assert!(self.bits[left] == 1 && self.bits[right] == 0 && length % 2 == 0);
        
        let m = next_step[left];
        seq[*idx] = m;
        *idx += 1;
        seq[*idx] = left;
        *idx += 1;
        self.compute_flip_seq_1_rec(seq, idx, left + 1, m - 1, next_step);
        seq[*idx] = left - 1;
        *idx += 1;
        seq[*idx] = m;
        *idx += 1;
        self.compute_flip_seq_1_rec(seq, idx, m + 1, right, next_step);
    }

    fn aux_pointers(&self, a: usize, b: usize, next_step: &mut [usize]) {
        assert!(a == b + 1 || (self.bits[a] == 1 && self.bits[b] == 0));
        let mut left_ustep_height = vec![0; b - a + 1];
        let mut height = 0;
        
        for i in a..=b {
            if self.bits[i] == 0 {
                assert!(height >= 1);
                let left = left_ustep_height[(height - 1) as usize];
                assert!(left < i);
                next_step[left] = i;
                next_step[i] = left;
            } else {
                assert!(height >= 0);
                left_ustep_height[height as usize] = i;
            }
            height += 2 * self.bits[i] - 1;
        }
        assert!(height == 0);
    }
}

impl fmt::Display for Vertex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for bit in &self.bits {
            write!(f, "{}", bit)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
struct Tree {
    root: usize,
    num_vertices: usize,
    children: Vec<LinkedList<usize>>,
    parent: Vec<usize>,
}

impl Tree {
    fn new(x: &Vertex) -> Self {
        let xv = &x.bits;
        assert!(xv.len() % 2 == 1);
        
        let root = 0;
        let num_vertices = (xv.len() - 1) / 2 + 1;
        let mut children = vec![LinkedList::new(); num_vertices];
        let mut parent = vec![0; num_vertices];
        
        let mut u = root;
        let mut n = 1;
        let mut height = 0;
        
        for i in 0..xv.len() - 1 {
            if xv[i] == 1 {
                children[u].push_back(n);
                parent[n] = u;
                u = n;
                n += 1;
            } else {
                u = parent[u];
            }
            height += 2 * xv[i] - 1;
            assert!(height >= 0);
        }
        assert!(n == num_vertices);
        
        Tree {
            root,
            num_vertices,
            children,
            parent,
        }
    }

    fn deg(&self, u: usize) -> usize {
        assert!(u < self.num_vertices);
        self.children[u].len() + if u == self.root { 0 } else { 1 }
    }

    fn num_children(&self, u: usize) -> usize {
        assert!(u < self.num_vertices);
        self.children[u].len()
    }

    fn ith_child(&self, u: usize, i: usize) -> usize {
        assert!(u < self.num_vertices);
        assert!(i < self.num_children(u));
        *self.children[u].iter().nth(i).unwrap()
    }

    fn is_tau_preimage(&self) -> bool {
        if self.num_vertices < 3 {
            return false;
        }
        let u = self.ith_child(self.root, 0);
        if self.num_children(u) == 0 {
            return false;
        }
        let v = self.ith_child(u, 0);
        self.num_children(v) == 0
    }

    fn is_tau_image(&self) -> bool {
        self.num_vertices >= 3 &&
        self.num_children(self.root) >= 2 &&
        self.num_children(self.ith_child(self.root, 0)) == 0
    }

    fn tau(&mut self) {
        assert!(self.is_tau_preimage());
        let u = self.ith_child(self.root, 0);
        let v = self.ith_child(u, 0);
        self.move_leaf(v, self.root, 0);
    }

    fn tau_inverse(&mut self) {
        assert!(self.is_tau_image());
        let v = self.ith_child(self.root, 0);
        let u = self.ith_child(self.root, 1);
        self.move_leaf(v, u, 0);
    }

    fn move_leaf(&mut self, leaf: usize, new_parent: usize, pos: usize) {
        assert!(leaf < self.num_vertices);
        assert!(new_parent < self.num_vertices);
        assert!(pos <= self.num_children(new_parent));
        assert!(self.num_children(leaf) == 0);
        
        let old_parent = self.parent[leaf];
        let index = self.children[old_parent]
            .iter()
            .position(|&x| x == leaf)
            .unwrap();
        let removed = self.children[old_parent].remove(index);
        assert_eq!(removed, leaf);
        
        let mut split = self.children[new_parent].split_off(pos);
        self.children[new_parent].push_back(leaf);
        self.children[new_parent].append(&mut split);
        
        self.parent[leaf] = new_parent;
    }

    fn rotate(&mut self) {
        assert!(self.num_vertices >= 2);
        let u = self.ith_child(self.root, 0);
        self.parent[self.root] = u;
        let moved_child = self.children[self.root].pop_front().unwrap();
        self.children[u].push_back(moved_child);
        self.children[u].back_mut().unwrap();
        self.root = u;
    }

    fn rotate_to_vertex(&mut self, u: usize) {
        while self.root != u {
            self.rotate();
        }
    }

    fn rotate_children(&mut self, k: usize) {
        let mut temp = VecDeque::from_iter(self.children[self.root].iter().cloned());
        temp.rotate_left(k);
        self.children[self.root] = LinkedList::from_iter(temp);
    }

    fn flip_tree(&mut self) -> bool {
        if self.is_tau_preimage() && self.is_flip_tree_tau() {
            self.tau();
            true
        } else if self.is_tau_image() {
            self.tau_inverse();
            if self.is_flip_tree_tau() {
                true
            } else {
                self.tau();
                false
            }
        } else {
            false
        }
    }

    fn root_canonically(&mut self) {
        let (c1, c2) = self.compute_center();
        if let Some(c2) = c2 {
            self.rotate_to_vertex(c1);
            while self.ith_child(self.root, 0) != c2 {
                self.rotate_children(1);
            }
            let x1 = self.to_bitstring();
            
            self.rotate();
            self.rotate_children(self.num_children(self.root) - 1);
            assert!(self.root == c2 && self.ith_child(self.root, 0) == c1);
            let x2 = self.to_bitstring();
            
            if bitstrings_less_than(&x1, &x2) {
                self.rotate();
                self.rotate_children(self.num_children(self.root) - 1);
                assert!(self.root == c1 && self.ith_child(self.root, 0) == c2);
            }
        } else {
            self.rotate_to_vertex(c1);
            let x = self.to_bitstring();
            let mut subtree_count = vec![0; x.len()];
            let mut c = 0;
            let mut depth = 0;
            
            for i in 0..x.len() {
                if x[i] == 1 {
                    depth += 1;
                } else {
                    depth -= 1;
                }
                subtree_count[i] = c;
                if depth == 0 {
                    c += 1;
                }
            }
            
            let k = self.min_string_rotation(&x);
            self.rotate_children(subtree_count[k]);
        }
    }

    fn compute_center(&self) -> (usize, Option<usize>) {
        let mut degs: Vec<usize> = (0..self.num_vertices).map(|u| self.deg(u)).collect();
        let mut leaves: Vec<usize> = (0..self.num_vertices).filter(|&u| degs[u] == 1).collect();
        let mut num_leaves = leaves.len();
        let mut num_vertices_remaining = self.num_vertices;
        let mut num_new_leaves = 0;
        
        while num_vertices_remaining > 2 {
            for i in 0..num_leaves {
                let u = leaves[i];
                for &child in &self.children[u] {
                    degs[child] -= 1;
                    if degs[child] == 1 {
                        leaves[num_new_leaves] = child;
                        num_new_leaves += 1;
                    }
                }
                if u != self.root {
                    degs[self.parent[u]] -= 1;
                    if degs[self.parent[u]] == 1 {
                        leaves[num_new_leaves] = self.parent[u];
                        num_new_leaves += 1;
                    }
                }
            }
            num_vertices_remaining -= num_leaves;
            num_leaves = num_new_leaves;
            num_new_leaves = 0;
        }
        
        assert!(num_leaves >= 1 && num_leaves <= 2);
        if num_leaves == 1 {
            (leaves[0], None)
        } else {
            (leaves[0], Some(leaves[1]))
        }
    }

    fn is_flip_tree_tau(&mut self) -> bool {
        if self.is_star() {
            return false;
        }
        
        let r = self.root;
        let u = self.ith_child(self.root, 0);
        let num_bits = 2 * (self.num_vertices - 1);
        
        let v = self.ith_child(self.root, 0);
        if self.num_children(v) == 1 && self.num_children(self.ith_child(v, 0)) == 0 {
            let this_bitstring = self.to_bitstring();
            self.root_canonically();
            let mut v = self.ith_child(self.root, 0);
            while !(self.num_children(v) == 1 && self.num_children(self.ith_child(v, 0)) == 0) {
                self.rotate();
                v = self.ith_child(self.root, 0);
            }
            let canon_bitstring = self.to_bitstring();
            
            self.rotate_to_vertex(r);
            while self.ith_child(self.root, 0) != u {
                self.rotate_children(1);
            }
            
            bitstrings_equal(&this_bitstring, &canon_bitstring)
        } else {
            if self.has_thin_leaf() {
                return false;
            }
            let mut v = self.ith_child(self.root, 0);
            let mut c = self.count_pending_edges(v);
            if c < self.num_children(v) || c < 2 || self.is_light_dumbbell() {
                return false;
            }
            let this_bitstring = self.to_bitstring();
            self.root_canonically();
            v = self.ith_child(self.root, 0);
            c = self.count_pending_edges(v);
            while c < self.num_children(v) || c < 2 {
                self.rotate();
                self.rotate_children(c);
                v = self.ith_child(self.root, 0);
                c = self.count_pending_edges(v);
            }
            let canon_bitstring = self.to_bitstring();
            
            self.rotate_to_vertex(r);
            while self.ith_child(self.root, 0) != u {
                self.rotate_children(1);
            }
            
            bitstrings_equal(&this_bitstring, &canon_bitstring)
        }
    }

    fn is_star(&self) -> bool {
        self.num_vertices <= 3 ||
        self.deg(self.root) == self.num_vertices - 1 ||
        self.deg(self.ith_child(self.root, 0)) == self.num_vertices - 1
    }

    fn is_light_dumbbell(&self) -> bool {
        if self.num_vertices < 5 {
            return false;
        }
        let u = self.ith_child(self.root, 0);
        let k = self.num_children(u);
        let l = self.num_children(self.root) - 1;
        k + l + 1 >= self.num_vertices - 1 && k > l
    }

    fn is_thin_leaf(&self, u: usize) -> bool {
        if self.deg(u) > 1 {
            return false;
        }
        (u == self.root && self.deg(self.ith_child(u, 0)) == 2) ||
        (u != self.root && self.deg(self.parent[u]) == 2)
    }

    fn has_thin_leaf(&self) -> bool {
        (0..self.num_vertices).any(|u| self.is_thin_leaf(u))
    }

    fn count_pending_edges(&self, u: usize) -> usize {
        let mut c = 0;
        for i in 0..self.num_children(u) {
            let v = self.ith_child(u, i);
            if self.num_children(v) == 0 {
                c += 1;
            } else {
                break;
            }
        }
        c
    }

    fn to_bitstring(&self) -> Vec<i32> {
        let mut x = Vec::new();
        self.to_bitstring_rec(&mut x, self.root);
        x
    }

    fn to_bitstring_rec(&self, x: &mut Vec<i32>, u: usize) {
        if self.num_children(u) == 0 {
            return;
        }
        for &child in &self.children[u] {
            x.push(1);
            self.to_bitstring_rec(x, child);
            x.push(0);
        }
    }

    fn min_string_rotation(&self, x: &[i32]) -> usize {
        let xx: Vec<i32> = x.iter().chain(x.iter()).cloned().collect();
        let mut fail = vec![-1; 2 * x.len()];
        let mut k = 0;
        
        for j in 1..2 * x.len() {
            let xj = xx[j];
            let mut i = fail[j - k - 1];
            
            while i != -1 && xj != xx[k + i as usize + 1] {
                if xj < xx[k + i as usize + 1] {
                    k = j - i as usize - 1;
                }
                i = fail[i as usize];
            }
            
            if xj != xx[k + i as usize + 1] {
                if xj < xx[k] {
                    k = j;
                }
                fail[j - k] = -1;
            } else {
                fail[j - k] = i + 1;
            }
        }
        k
    }
}

fn bitstrings_less_than(x: &[i32], y: &[i32]) -> bool {
    for (xi, yi) in x.iter().zip(y.iter()) {
        if xi < yi {
            return true;
        } else if xi > yi {
            return false;
        }
    }
    false
}

fn bitstrings_equal(x: &[i32], y: &[i32]) -> bool {
    x.iter().zip(y.iter()).all(|(a, b)| a == b)
}

struct HamCycle {
    x: Vertex,
    y: Vertex,
    limit: i64,
    visit_f: fn(&[i32], usize),
    length: i64,
}

impl HamCycle {
    fn new(x: Vertex, limit: i64, visit_f: fn(&[i32], usize)) -> Self {
        assert!(x.len() % 2 == 1);
        let n = x.len() / 2;
        
        let mut xs = x.clone();
        let mut skip = 0;
        
        if xs.bits[2 * n] == 1 {
            xs.rev_inv(None, None);
            skip += xs.to_last_vertex() as i64;
            xs.rev_inv(None, None);
            xs.bits[2 * n] = 0;
            skip += 1;
        }
        
        skip += xs.to_first_vertex() as i64;
        assert!(xs.is_first_vertex());
        
        let mut y = xs.clone();
        let mut y_tree = Tree::new(&y);
        
        if skip > 0 && y_tree.flip_tree() {
            if xs.bits[1] == 1 && skip <= 5 {
                skip = 6 - skip;
            }
            let y_string = y_tree.to_bitstring();
            let mut y_vec = y_string;
            y_vec.push(0);
            xs = Vertex::new(y_vec);
            y = xs.clone();
        }
        
        HamCycle {
            x,
            y,
            limit,
            visit_f,
            length: 0,
        }
    }

    fn flip_seq(
        &mut self,
        seq: &[usize],
        dist_to_start: &mut i64,
        final_path: bool,
    ) -> bool {
        if *dist_to_start > 0 || final_path || (self.limit >= 0 && self.length + seq.len() as i64 >= self.limit) {
            for &i in seq {
                if (final_path && *dist_to_start == 0) || (self.limit >= 0 && self.length == self.limit) {
                    return true;
                }
                
                self.y.bits[i] = 1 - self.y.bits[i];
                (self.visit_f)(&self.y.bits, i);
                self.length += 1;
                
                if *dist_to_start > 0 {
                    *dist_to_start -= 1;
                }
            }
        } else {
            for &i in seq {
                self.y.bits[i] = 1 - self.y.bits[i];
                (self.visit_f)(&self.y.bits, i);
            }
            self.length += seq.len() as i64;
        }
        false
    }
}

fn help() {
    println!(
        "./middle [options]  compute middle levels Gray code from [Muetze,Nummenpalo]
-h                  display this help
-n<1,2,...>         list bitstrings of length 2n+1 with weight n or n+1
-l<-1,0,1,2,...>    number of bitstrings to list; -1 for full cycle
-v<0,1>^{2n+1}      initial bitstring (length 2n+1, weight n or n+1)
-s<0,1>             store and print all visited bitstrings (no=0, yes=1)
-p<0,1>             print the flip positions instead of bitstrings (no=0, yes=1)
examples:  ./middle -n2
           ./middle -n2 -v01010
           ./middle -n2 -p1
           ./middle -n10 -l50
           ./middle -n12 -s0"
    );
}

fn opt_n_missing() {
    eprintln!("option -n is mandatory and must come before -v");
}

fn opt_v_error() {
    eprintln!("option -v must be followed by a bitstring of length 2n+1 with weight n or n+1");
}

fn visit_f_empty(_y: &[i32], _i: usize) {}

static mut FLIP_SEQ: Vec<usize> = Vec::new();

fn visit_f_log(y: &[i32], i: usize) {
    unsafe {
        FLIP_SEQ.push(i);
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 0;
    let mut n_set = false;
    let mut limit = -1;
    let mut v = Vec::new();
    let mut v_set = false;
    let mut store_vertices = true;
    let mut print_flip_pos = false;
    
    let mut opts = getopts::Options::new();
    opts.optflag("h", "help", "display this help");
    opts.optopt("n", "", "dimension parameter", "1,2,...");
    opts.optopt("l", "", "limit", "-1,0,1,2,...");
    opts.optopt("v", "", "initial vertex", "bitstring");
    opts.optopt("s", "", "store vertices", "0,1");
    opts.optopt("p", "", "print flip positions", "0,1");
    
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            eprintln!("{}", f);
            help();
            std::process::exit(1);
        }
    };
    
    if matches.opt_present("h") {
        help();
        return;
    }
    
    if let Some(n_str) = matches.opt_str("n") {
        n = match n_str.parse() {
            Ok(num) if num >= 1 => num,
            _ => {
                eprintln!("option -n must be followed by an integer from {{1,2,...}}");
                std::process::exit(1);
            }
        };
        v = vec![0; 2 * n + 1];
        n_set = true;
    }
    
    if let Some(l_str) = matches.opt_str("l") {
        limit = match l_str.parse() {
            Ok(num) if num >= -1 => num,
            _ => {
                eprintln!("option -l must be followed by an integer from {{-1,0,1,2,...}}");
                std::process::exit(1);
            }
        };
    }
    
    if let Some(v_str) = matches.opt_str("v") {
        if !n_set {
            opt_n_missing();
            help();
            std::process::exit(1);
        }
        
        let mut xv = Vec::new();
        let mut num_ones = 0;
        
        for c in v_str.chars() {
            match c {
                '0' => {
                    xv.push(0);
                }
                '1' => {
                    xv.push(1);
                    num_ones += 1;
                }
                _ => {
                    opt_v_error();
                    std::process::exit(1);
                }
            }
            if xv.len() > 2 * n + 1 {
                opt_v_error();
                std::process::exit(1);
            }
        }
        
        if xv.len() != 2 * n + 1 || num_ones < n || num_ones > n + 1 {
            opt_v_error();
            std::process::exit(1);
        }
        
        v = xv;
        v_set = true;
    }
    
    if let Some(s_str) = matches.opt_str("s") {
        store_vertices = match s_str.parse() {
            Ok(0) => false,
            Ok(1) => true,
            _ => {
                eprintln!("option -s must be followed by 0 or 1");
                std::process::exit(1);
            }
        };
    }
    
    if let Some(p_str) = matches.opt_str("p") {
        print_flip_pos = match p_str.parse() {
            Ok(0) => false,
            Ok(1) => true,
            _ => {
                eprintln!("option -p must be followed by 0 or 1");
                std::process::exit(1);
            }
        };
    }
    
    if !n_set {
        opt_n_missing();
        help();
        std::process::exit(1);
    }
    
    if !v_set {
        v = vec![1; n];
        v.extend(vec![0; n + 1]);
    }
    
    let x = Vertex::new(v);
    let visit_f = if store_vertices { visit_f_log } else { visit_f_empty };
    
    unsafe {
        FLIP_SEQ.clear();
    }
    
    let mut hc = HamCycle::new(x.clone(), limit, visit_f);
    
    if store_vertices {
        if limit != 0 {
            println!("{}", x);
        }
        
        unsafe {
            for &i in &FLIP_SEQ[..FLIP_SEQ.len().saturating_sub(1)] {
                hc.y.bits[i] = 1 - hc.y.bits[i];
                if print_flip_pos {
                    println!("{}", i);
                } else {
                    println!("{}", hc.y);
                }
            }
            
            if limit == FLIP_SEQ.len() as i64 {
                println!("output limit reached");
            }
        }
    }
}
