// Define an enum to represent the type of a vertex.
#[derive(Debug, PartialEq)]
pub enum VertexType {
    Corner,
    Top,
    Bottom,
    Left,
    Right,
    None,
}

// Define a struct for the Vertex.
#[derive(Debug)]
pub struct Vertex {
    pub north: i32,
    pub east: i32,
    pub south: i32,
    pub west: i32,
    pub type_: VertexType,
}

impl Default for Vertex {
    fn default() -> Self {
        Self::new()
    }
}

impl Vertex {
    // Constructor for Vertex with default type None.
    pub fn new() -> Self {
        Vertex {
            north: 0,
            east: 0,
            south: 0,
            west: 0,
            type_: VertexType::None,
        }
    }

    // Method to initialize the vertex with given coordinates and determine its type.
    pub fn init(&mut self, north: i32, east: i32, south: i32, west: i32) {
        self.north = north;
        self.east = east;
        self.south = south;
        self.west = west;
        let zeros = (self.north == 0) as i32
            + (self.south == 0) as i32
            + (self.west == 0) as i32
            + (self.east == 0) as i32;
        match zeros {
            3..=std::i32::MAX | 0 => self.type_ = VertexType::None,
            2 => self.type_ = VertexType::Corner,
            _ if self.south == 0 => self.type_ = VertexType::Top,
            _ if self.north == 0 => self.type_ = VertexType::Bottom,
            _ if self.east == 0 => self.type_ = VertexType::Left,
            _ if self.west == 0 => self.type_ = VertexType::Right,
            _ => unreachable!(), // This case should never occur due to previous checks.
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vertex_init_corner() {
        let mut vertex = Vertex::new();
        vertex.init(0, 1, 9, 8);
        assert_eq!(vertex.type_, VertexType::Bottom);
    }

    #[test]
    fn test_vertex_init_top() {
        let mut vertex = Vertex::new();
        vertex.init(0, 1, 0, 8);
        assert_eq!(vertex.type_, VertexType::Corner);
    }

    #[test]
    fn test_vertex_init_default() {
        let vertex = Vertex::new();
        assert_eq!(vertex.type_, VertexType::None);
    }
}
