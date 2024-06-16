#[derive(Debug, Clone, PartialEq)]
pub enum EdgeDir {
   Hor, Ver, None
}

#[derive(Debug, Clone)]
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
}

impl Edge {
    pub fn init(&mut self, dir: EdgeDir, tail: i32, head: i32, prev: i32, next: i32, left: i32, right: i32, wall: i32) {
        self.dir = dir;
        self.tail = tail;
        self.head = head;
        self.left = left;
        self.right = right;
        self.wall = wall;
        self.prev = prev;
        self.next = next;
    }
}

#[cfg(test)]
mod tests {
    use super::*; // bring everything from parent scope into this one

    #[test]
    fn test_edge_init() {
        let mut edge = Edge::default();
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
