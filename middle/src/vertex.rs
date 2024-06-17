

// A struct to represent and manipulate a vertex of the cube of odd dimension 2n+1
#[derive(Debug, Clone)]
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
    // constructor
    pub fn new(x: Vec<i32>) -> Self {
        assert!(x.len() % 2 == 1);
        assert!(x.len() >= 3);
        Vertex { bits: x }
    }

    pub fn get_bits(&self) -> &Vec<i32> {
        &self.bits
    }

    // reverse and invert bitstring
    pub fn rev_inv(&mut self) {
        let len = self.bits.len();
        for i in 0..len - 1 {
            self.bits[i] = 1 - self.bits[i];
        }
        self.bits[..len - 2].reverse();
    }

    // compute whether vertex is first or last vertex on a path
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
        assert_eq!(v.is_first_vertex(), true);
        let v = Vertex::new(vec![1, 1, 0, 1, 0]);
        assert_eq!(v.is_first_vertex(), false);
    }
}
