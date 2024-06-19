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
    dir_: EdgeDir,
    tail_: i32,
    head_: i32,
    prev_: i32,
    next_: i32,
    left_: i32,
    right_: i32,
    wall_: i32,
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
            dir_: EdgeDir::None,
            tail_: 0,
            head_: 0,
            prev_: 0,
            next_: 0,
            left_: 0,
            right_: 0,
            wall_: 0,
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
        self.dir_ = dir;
        self.tail_ = tail;
        self.head_ = head;
        self.prev_ = prev;
        self.next_ = next;
        self.left_ = left;
        self.right_ = right;
        self.wall_ = wall;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_initialization() {
        let mut edge = Edge::new();
        edge.init(EdgeDir::Hor, 1, 2, 3, 4, 5, 6, 7);

        assert_eq!(edge.dir_, EdgeDir::Hor);
        assert_eq!(edge.tail_, 1);
        assert_eq!(edge.head_, 2);
        assert_eq!(edge.prev_, 3);
        assert_eq!(edge.next_, 4);
        assert_eq!(edge.left_, 5);
        assert_eq!(edge.right_, 6);
        assert_eq!(edge.wall_, 7);
    }
}
