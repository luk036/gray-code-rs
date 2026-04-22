/// A struct to represent and manipulate a vertex of the cube of odd dimension 2n+1.
///
/// This vertex is represented as a bitstring of length 2n+1 where n is
/// the dimension parameter. The bitstring encodes the position within
/// the Gray code structure.
#[derive(Debug, Clone, Default)]
pub struct Vertex {
    bits: Vec<i32>,
}

impl PartialEq for Vertex {
    fn eq(&self, other: &Self) -> bool {
        self.bits == other.bits
    }
}

impl Eq for Vertex {}

impl std::fmt::Display for Vertex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.bits)
    }
}

impl Vertex {
    /// Constructs a new `Vertex` from a bitstring vector.
    ///
    /// # Panics
    ///
    /// Panics if the bitstring length is not odd or less than 3.
    ///
    /// # Arguments
    ///
    /// * `x` - A vector of integers (0 or 1) representing the bitstring.
    pub fn new(x: Vec<i32>) -> Self {
        assert!(x.len() % 2 == 1);
        assert!(x.len() >= 3);
        Vertex { bits: x }
    }

    /// Returns an immutable reference to the underlying bitstring.
    pub fn get_bits(&self) -> &Vec<i32> {
        &self.bits
    }

    /// Returns a mutable reference to the underlying bitstring.
    pub fn get_bits_mut(&mut self) -> &mut Vec<i32> {
        &mut self.bits
    }

    /// Reverses and inverts the bitstring.
    ///
    /// This operation flips each bit (0 becomes 1, 1 becomes 0) except
    /// for the last element, then reverses the order of the first len-2 elements.
    pub fn rev_inv(&mut self) {
        let len = self.bits.len();
        for i in 0..len - 1 {
            self.bits[i] = 1 - self.bits[i];
        }
        self.bits[..len - 2].reverse();
    }

    /// Computes whether this vertex is the first vertex on a path.
    ///
    /// A vertex is considered a "first vertex" if the cumulative height
    /// (computed as sum of 2*bit - 1 for each position) reaches zero
    /// before the end of the bitstring.
    pub fn is_first_vertex(&self) -> bool {
        let mut height = 0;
        for i in &self.bits[..self.bits.len() - 1] {
            height += 2 * i - 1;
            if height == 0 {
                return true;
            }
        }
        false
    }

    /// Converts the vertex to the last vertex position.
    ///
    /// Returns the number of steps taken.
    pub fn to_last_vertex(&mut self) -> usize {
        todo!("to_last_vertex implementation pending")
    }

    /// Converts the vertex to the first vertex position.
    ///
    /// Returns the number of steps taken.
    pub fn to_first_vertex(&mut self) -> usize {
        todo!("to_first_vertex implementation pending")
    }

    /// Flips the bit at the specified index.
    ///
    /// # Arguments
    ///
    /// * `i` - The index of the bit to flip.
    pub fn flip_bit(&mut self, i: usize) {
        if i < self.bits.len() {
            self.bits[i] = 1 - self.bits[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let v = Vertex::new(vec![1, 0, 1]);
        assert_eq!(v.get_bits(), &vec![1, 0, 1]);
    }

    #[test]
    fn test_rev_inv() {
        let mut v = Vertex::new(vec![1, 0, 1]);
        v.rev_inv();
        assert_eq!(v.get_bits(), &vec![0, 1, 1]);
    }

    #[test]
    fn test_is_first_vertex() {
        let v = Vertex::new(vec![1, 0, 1]);
        assert!(v.is_first_vertex());
        let v = Vertex::new(vec![1, 1, 0, 1, 0]);
        assert!(!v.is_first_vertex());
    }
}
