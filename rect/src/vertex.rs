// The type indicates where the T is pointing to.
// The 4 corner vertices of the rectangulation have type 'corner'.
#[derive(Debug, Clone, PartialEq)]
pub enum VertexType {
    Corner,
    Top,
    Bottom,
    Left,
    Right,
    None,
}

#[derive(Debug, Clone)]
pub struct Vertex {
    pub north: i32,
    pub east: i32,
    pub south: i32,
    pub west: i32,
    pub type_: VertexType,
}

impl Default for Vertex {
    fn default() -> Self {
        Self {
            north: 0,
            east: 0,
            south: 0,
            west: 0,
            type_: VertexType::None,
        }
    }
}

impl Vertex {
    pub fn init(&mut self, north: i32, east: i32, south: i32, west: i32) {
        self.north = north;
        self.east = east;
        self.south = south;
        self.west = west;

        let zeros = (self.north == 0) as i32
            + (self.south == 0) as i32
            + (self.west == 0) as i32
            + (self.east == 0) as i32;

        if zeros >= 3 || zeros <= 0 {
            self.type_ = VertexType::None;
        } else if zeros == 2 {
            self.type_ = VertexType::Corner;
        } else if self.south == 0 {
            self.type_ = VertexType::Top;
        } else if self.north == 0 {
            self.type_ = VertexType::Bottom;
        } else if self.east == 0 {
            self.type_ = VertexType::Left;
        } else if self.west == 0 {
            self.type_ = VertexType::Right;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*; // bring everything from parent scope into this one

    #[test]
    fn test_vertex_default() {
        let vertex = Vertex::default();

        assert_eq!(vertex.north, 0);
        assert_eq!(vertex.east, 0);
        assert_eq!(vertex.south, 0);
        assert_eq!(vertex.west, 0);
        assert_eq!(vertex.type_, VertexType::None);
    }

    #[test]
    fn test_vertex_init() {
        let mut vertex = Vertex::default();
        vertex.init(1, 2, 3, 4);

        assert_eq!(vertex.north, 1);
        assert_eq!(vertex.east, 2);
        assert_eq!(vertex.south, 3);
        assert_eq!(vertex.west, 4);
        assert_eq!(vertex.type_, VertexType::None);
    }

    #[test]
    fn test_vertex_init_corner() {
        let mut vertex = Vertex::default();
        vertex.init(0, 2, 3, 4);

        assert_eq!(vertex.type_, VertexType::Bottom);

        vertex.init(1, 0, 3, 4);

        assert_eq!(vertex.type_, VertexType::Left);
    }

    #[test]
    fn test_vertex_init_top() {
        let mut vertex = Vertex::default();
        vertex.init(1, 2, 0, 4);

        assert_eq!(vertex.type_, VertexType::Top);
    }

    // Similar tests can be written for other VertexTypes (Bottom, Left, Right) by changing the init method arguments
}
