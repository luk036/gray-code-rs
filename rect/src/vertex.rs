/// Represents the type of a vertex in a rectangulation.
///
/// Different vertex types indicate the position and connectivity
/// of a vertex within the rectangular partition.
#[derive(Debug, PartialEq)]
pub enum VertexType {
    /// A corner vertex where two edges meet at a rectangle corner.
    Corner,
    /// A top vertex lying on the top boundary.
    Top,
    /// A bottom vertex lying on the bottom boundary.
    Bottom,
    /// A left vertex lying on the left boundary.
    Left,
    /// A right vertex lying on the right boundary.
    Right,
    /// A generic vertex with no special boundary position.
    None,
}

/// Represents a vertex in a rectangulation with coordinates and type.
///
/// Each vertex has four coordinate values representing its position
/// relative to the four boundaries (north, east, south, west).
#[derive(Debug)]
pub struct Vertex {
    /// The north coordinate of the vertex.
    pub north: i32,
    /// The east coordinate of the vertex.
    pub east: i32,
    /// The south coordinate of the vertex.
    pub south: i32,
    /// The west coordinate of the vertex.
    pub west: i32,
    /// The type of the vertex (Corner, Top, Bottom, Left, Right, None).
    pub type_: VertexType,
}

impl Default for Vertex {
    fn default() -> Self {
        Self::new()
    }
}

impl Vertex {
    /// Constructs a new `Vertex` with default coordinates and type `None`.
    pub fn new() -> Self {
        Vertex {
            north: 0,
            east: 0,
            south: 0,
            west: 0,
            type_: VertexType::None,
        }
    }

    /// Initializes the vertex with given coordinates and determines its type.
    ///
    /// The vertex type is determined by counting how many coordinates are zero:
    /// - 3 or more zeros: `None`
    /// - Exactly 2 zeros: `Corner`
    /// - 1 zero on south boundary: `Top`
    /// - 1 zero on north boundary: `Bottom`
    /// - 1 zero on east boundary: `Left`
    /// - 1 zero on west boundary: `Right`
    ///
    /// # Arguments
    ///
    /// * `north` - The north coordinate.
    /// * `east` - The east coordinate.
    /// * `south` - The south coordinate.
    /// * `west` - The west coordinate.
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
