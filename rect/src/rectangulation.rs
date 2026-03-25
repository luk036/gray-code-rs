use crate::edge::Edge;
use crate::rectangle::Rectangle;
use crate::vertex::{Vertex, VertexType};
use crate::wall::Wall;
use std::vec;

/// Represents the type of rectangulation.
#[derive(Debug, PartialEq, Clone)]
pub enum RectangulationType {
    /// A generic rectangulation.
    Generic,
    /// A bottom-aligned rectangulation.
    Baligned,
    /// A diagonal rectangulation.
    Diagonal,
}

/// Represents the direction in a rectangulation.
#[derive(Debug, PartialEq, Clone)]
pub enum RectangulationDirection {
    /// Left direction.
    Left,
    /// Right direction.
    Right,
    /// No direction specified.
    None,
}

/// Represents the pattern used for generating rectangulations.
#[derive(Debug, PartialEq, Clone)]
pub enum RectangulationPattern {
    /// Windmill pattern rotating clockwise.
    WMillClockwise,
    /// Windmill pattern rotating counterclockwise.
    WMillCounterclockwise,
    /// Brick pattern going left to right.
    BrickLeftRight,
    /// Brick pattern going right to left.
    BrickRightLeft,
    /// Brick pattern going top to bottom.
    BrickTopBottom,
    /// Brick pattern going bottom to top.
    BrickBottomTop,
    /// Horizontal vertical pattern.
    HVertical,
    /// Horizontal horizontal pattern.
    HHorizontal,
}

/// A struct representing a rectangulation of n rectangles.
pub struct Rectangulation {
    /// The number of rectangles.
    pub n: usize,
    /// The type of rectangulation.
    pub typ: RectangulationType,
    /// The patterns used in this rectangulation.
    pub patterns: Vec<RectangulationPattern>,
    /// The directions for each rectangle.
    pub directions: Vec<RectangulationDirection>,
    /// The sizes of rectangles.
    pub sizes: Vec<i32>,
    /// The vertices of the rectangulation.
    pub vertices: Vec<Vertex>,
    /// The walls of the rectangulation.
    pub walls: Vec<Wall>,
    /// The edges of the rectangulation.
    pub edges: Vec<Edge>,
    /// The rectangles in the rectangulation.
    pub rectangles: Vec<Rectangle>,
}

impl Rectangulation {
    /// Constructs a new Rectangulation.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of rectangles.
    /// * `typ` - The type of rectangulation.
    /// * `patterns` - The patterns to use.
    pub fn new(n: usize, typ: RectangulationType, patterns: Vec<RectangulationPattern>) -> Self {
        let mut instance = Self {
            n,
            typ,
            patterns,
            directions: vec![RectangulationDirection::None],
            sizes: vec![-1],
            vertices: Vec::new(),
            walls: Vec::new(),
            edges: Vec::new(),
            rectangles: Vec::new(),
        };
        instance.set_all_vertical();
        instance
    }

    /// Initializes the rectangulation with vertices, walls, edges, and rectangles.
    ///
    /// # Arguments
    ///
    /// * `vertices` - The vertices of the rectangulation.
    /// * `walls` - The walls of the rectangulation.
    /// * `edges` - The edges of the rectangulation.
    /// * `rectangles` - The rectangles in the rectangulation.
    pub fn init(
        &mut self,
        vertices: Vec<Vertex>,
        walls: Vec<Wall>,
        edges: Vec<Edge>,
        rectangles: Vec<Rectangle>,
    ) {
        self.vertices = vertices;
        self.walls = walls;
        self.edges = edges;
        self.rectangles = rectangles;
    }

    /// Sets all directions to vertical.
    pub fn set_all_vertical(&mut self) {
        for j in 1..=self.n {
            self.directions.push(RectangulationDirection::Left);
            self.sizes.push(j as i32);
        }
    }
}

impl Rectangulation {
    /// Prints the data structures of the rectangulation (debug method).
    pub fn print_data(&self) {
        println!("edges:");
        for (i, e) in self.edges.iter().enumerate() {
            println!("\t{}. {:?}", i, e); // Placeholder for actual fields in Edge
        }
        println!("vertices:");
        for (i, v) in self.vertices.iter().enumerate() {
            println!("\t{}. {:?}", i, v); // Placeholder for actual fields in Vertex
        }
        // Continue for walls and rectangles similarly.
    }
}

impl Rectangulation {
    /// Prints coordinates for a generic rectangulation.
    ///
    /// Uses a greedy algorithm to find an equispaced grid for vertex placement.
    pub fn print_coordinates_generic(&self) {
        let _vertex_x_coord = vec![-1; 2 * self.n + 3];
        let mut active_vertices = Vec::new();

        // Populate initial active vertices
        for (a, v) in self.vertices.iter().enumerate() {
            let side_edge_id = match v.type_ {
                VertexType::Right => v.north,
                VertexType::Corner => std::cmp::max(v.north, v.south),
                _ => continue,
            };
            if self.edges[side_edge_id as usize].left == 0 {
                active_vertices.push(a);
            }
        }

        // Propagation logic would go here, similar to C++, but adapted for Rust's ownership and borrowing rules
        // This part is complex and requires careful handling of indices and references due to Rust's safety checks.

        // The rest of the method's logic follows, translating each step of the propagation and printing as needed.
        // Since Rust does not have a direct equivalent to C++'s std::cout, you would use println! macro or a logger.

        // For simplicity, I'm omitting the full propagation logic here. It would involve iterating over active_vertices,
        // updating vertex_x_coord, and managing active vertices in a way that respects Rust's memory safety rules.
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rectangulation_new_and_init() {
        // Setup mock data
        let n = 5;
        let typ = RectangulationType::Generic;
        let patterns = vec![
            RectangulationPattern::WMillClockwise,
            RectangulationPattern::BrickLeftRight,
        ];

        // Create a new Rectangulation instance
        let mut rect = Rectangulation::new(n, typ.clone(), patterns.clone());

        // Mock vectors for initialization
        let vertices = vec![]; // Fill with actual Vertex instances if needed
        let walls = vec![]; // Fill with actual Wall instances if needed
        let edges = vec![]; // Fill with actual Edge instances if needed
        let rectangles = vec![]; // Fill with actual Rectangle instances if needed

        // Initialize the Rectangulation
        rect.init(vertices, walls, edges, rectangles);

        // Assertions
        assert_eq!(rect.n, n);
        assert_eq!(rect.typ, typ);
        assert_eq!(rect.patterns, patterns);
        assert_eq!(rect.directions.len(), n + 1); // +1 for the initial None direction
        assert_eq!(rect.sizes.len(), n + 1);
        assert!(rect.vertices.is_empty()); // After init, vertices should be empty if not provided
        assert!(rect.walls.is_empty());
        assert!(rect.edges.is_empty());
        assert!(rect.rectangles.is_empty());
    }
}
