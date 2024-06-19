/// Represents the direction of an edge.
#[derive(Debug, PartialEq)]
pub enum EdgeDir {
    Hor,  // Horizontal
    Ver,  // Vertical
    None, // No direction
}

/// Represents an edge in a graph with additional properties like direction, endpoints, and adjacency info.
#[derive(Debug)]
pub struct Edge {
    pub dir: EdgeDir,
    pub tail: i32,
    pub head: i32,
    pub prev: i32,
    pub next: i32,
    pub left: i32,
    pub right: i32,
    pub wall: i32,
}

impl Default for Edge {
    fn default() -> Self {
        Self::new()
    }
}

impl Edge {
    /// Constructs a new `Edge` with default values.
    pub fn new() -> Self {
        Edge {
            dir: EdgeDir::None,
            tail: 0,
            head: 0,
            prev: 0,
            next: 0,
            left: 0,
            right: 0,
            wall: 0,
        }
    }

    /// Initializes the edge with provided parameters.
    pub fn init(
        &mut self,
        dir: EdgeDir,
        tail: i32,
        head: i32,
        prev: i32,
        next: i32,
        left: i32,
        right: i32,
        wall: i32,
    ) {
        self.dir = dir;
        self.tail = tail;
        self.head = head;
        self.prev = prev;
        self.next = next;
        self.left = left;
        self.right = right;
        self.wall = wall;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn testedgeinitialization() {
        let mut edge = Edge::new();
        edge.init(EdgeDir::Hor, 1, 2, 3, 4, 5, 6, 7);

        assert_eq!(edge.dir, EdgeDir::Hor);
        assert_eq!(edge.tail, 1);
        assert_eq!(edge.head, 2);
        assert_eq!(edge.prev, 3);
        assert_eq!(edge.next, 4);
        assert_eq!(edge.left, 5);
        assert_eq!(edge.right, 6);
        assert_eq!(edge.wall, 7);
    }
}
