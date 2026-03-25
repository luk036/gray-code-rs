mod edge;
mod rectangle;
mod vertex;
mod wall;

use std::vec::Vec;
use RectangulationPattern::*;
use RectangulationType::*;
use RectangulationDirection::*;

/// Represents the type of rectangulation.
#[derive(Debug, PartialEq)]
enum RectangulationType {
    /// A generic rectangulation.
    Generic,
    /// A bottom-aligned rectangulation.
    Baligned,
    /// A diagonal rectangulation.
    Diagonal,
}

/// Represents the direction in a rectangulation.
#[derive(Debug, PartialEq)]
enum RectangulationDirection {
    /// Left direction.
    Left,
    /// Right direction.
    Right,
    /// No direction specified.
    None,
}

/// Represents the pattern used for generating rectangulations.
#[derive(Debug, PartialEq)]
enum RectangulationPattern {
    /// Windmill pattern rotating clockwise.
    WMillClockwise,
    /// Windmill pattern rotating counterclockwise.
    WMillCounterclockwise,
    /// Brick pattern going left to right.
    BrickLeftright,
    /// Brick pattern going right to left.
    BrickRightleft,
    /// Brick pattern going top to bottom.
    BrickTopbottom,
    /// Brick pattern going bottom to top.
    BrickBottomtop,
    /// Horizontal vertical pattern.
    HVertical,
    /// Horizontal horizontal pattern.
    HHorizontal,
}

/// A struct representing a rectangulation of n rectangles.
pub struct Rectangulation {
    n: i32,
    type_: RectangulationType,
    patterns: Vec<RectangulationPattern>,
    o: Vec<RectangulationDirection>,
    s: Vec<i32>,
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
    pub fn new(n: i32, type_: RectangulationType, patterns: Vec<RectangulationPattern>) -> Self {
        let mut rectangulation = Rectangulation {
            n,
            type_,
            patterns,
            o: vec![],
            s: vec![],
            vertices: vec![],
            walls: vec![],
            edges: vec![],
            rectangles: vec![],
        };

        // M1 - Initialize
        rectangulation.set_all_vertical();
        rectangulation.o.push(None);
        rectangulation.s.push(-1);
        for j in 1..=n {
            rectangulation.o.push(Left);
            rectangulation.s.push(j);
        }

        rectangulation
    }

    fn set_all_vertical(&mut self) {
        // Implement this method according to your needs
    }

    // Continue with other methods...
}

impl Rectangulation {
    pub fn new(n: i32, type_: RectangulationType, patterns: Vec<RectangulationPattern>) -> Self {
        let mut rectangulation = Rectangulation {
            n,
            type_,
            patterns,
            o: vec![],
            s: vec![],
            vertices: vec![],
            walls: vec![],
            edges: vec![],
            rectangles: vec![],
        };

        // M1 - Initialize
        rectangulation.set_all_vertical();
        rectangulation.o.push(None);
        rectangulation.s.push(-1);
        for j in 1..=n {
            rectangulation.o.push(Left);
            rectangulation.s.push(j);
        }

        rectangulation
    }

    fn set_all_vertical(&mut self) {
        // Implement this method according to your needs
    }

    // Continue with other methods...
}

/// Represents the direction of an edge.
pub enum EdgeDir {
    /// Horizontal direction.
    Hor,
    /// Vertical direction.
    Ver,
    /// No direction.
    None,
}

/// Represents the type of a vertex.
pub enum VertexType {
    /// No specific type.
    None,
    /// A corner vertex.
    Corner,
    /// A bottom vertex.
    Bottom,
    /// A top vertex.
    Top,
    /// A left vertex.
    Left,
    /// A right vertex.
    Right,
}

/// Represents the type of rectangulation.
pub enum RectangulationType {
    /// A generic rectangulation.
    Generic,
    /// A diagonal rectangulation.
    Diagonal,
    /// A bottom-aligned rectangulation.
    Baligned,
}

/// Represents the pattern used for generating rectangulations.
pub enum RectangulationPattern {
    /// Brick pattern going left to right.
    BrickLeftright,
    /// Brick pattern going right to left.
    BrickRightleft,
    /// Brick pattern going bottom to top.
    BrickBottomtop,
    /// Brick pattern going top to bottom.
    BrickTopbottom,
    /// Windmill pattern rotating clockwise.
    WMillClockwise,
    /// Windmill pattern rotating counterclockwise.
    WMillCounterclockwise,
    /// Horizontal vertical pattern.
    HVertical,
    /// Horizontal horizontal pattern.
    HHorizontal,
}

/// Represents the direction in a rectangulation.
pub enum RectangulationDirection {
    /// Left direction.
    Left,
    /// Right direction.
    Right,
}

/// Represents an edge in a graph.
pub struct Edge {
    dir: EdgeDir,
    tail: i32,
    head: i32,
    prev: i32,
    next: i32,
    left: i32,
    right: i32,
    wall: i32,
}

/// Represents a vertex in a rectangulation.
pub struct Vertex {
    north: i32,
    east: i32,
    south: i32,
    west: i32,
    type_: VertexType,
}

/// Represents a rectangle with four corners.
pub struct Rectangle {
    nwest: i32,
    swest: i32,
    neast: i32,
    seast: i32,
}

/// Represents a wall with start and end points.
pub struct Wall {
    first: i32,
    last: i32,
}

/// A struct representing a rectangulation of n rectangles.
pub struct Rectangulation {
    n: i32,
    type_: RectangulationType,
    patterns: Vec<RectangulationPattern>,
    o: Vec<RectangulationDirection>,
    s: Vec<i32>,
    vertices: Vec<Vertex>,
    walls: Vec<Wall>,
    edges: Vec<Edge>,
    rectangles: Vec<Rectangle>,
}

impl Rectangulation {
    /// Constructs a new Rectangulation.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of rectangles.
    /// * `type_` - The type of rectangulation.
    /// * `patterns` - The patterns to use.
    pub fn new(n: i32, type_: RectangulationType, patterns: Vec<RectangulationPattern>) -> Self {
        Rectangulation {
            n,
            type_,
            patterns,
            o: vec![],
            s: vec![],
            vertices: vec![],
            walls: vec![],
            edges: vec![],
            rectangles: vec![],
        }
    }

    /// Initializes the rectangulation with vertices, walls, edges, and rectangles.
    ///
    /// # Arguments
    ///
    /// * `vertices` - The vertices of the rectangulation.
    /// * `walls` - The walls of the rectangulation.
    /// * `edges` - The edges of the rectangulation.
    /// * `rectangles` - The rectangles in the rectangulation.
    pub fn init(&mut self, vertices: &Vec<Vertex>, walls: &Vec<Wall>, edges: &Vec<Edge>, rectangles: &Vec<Rectangle>) {
        self.vertices = vertices.clone();
        self.walls = walls.clone();
        self.edges = edges.clone();
        self.rectangles = rectangles.clone();
    }

    /// Sets all directions to vertical and initializes all data structures.
    pub fn set_all_vertical(&mut self) {
        // initialize the 3n+1 edges, 2n+2 vertices, n rectangles and n+3 walls
        self.edges = vec![Edge { dir: EdgeDir::None, tail: 0, head: 0, prev: 0, next: 0, left: 0, right: 0, wall: 0 }; 3 * self.n as usize + 2];
        self.vertices = vec![Vertex { north: 0, east: 0, south: 0, west: 0, type_: VertexType::None }; 2 * self.n as usize + 3];
        self.rectangles = vec![Rectangle { nwest: 0, swest: 0, neast: 0, seast: 0 }; self.n as usize + 1];
        self.walls = vec![Wall { first: 0, last: 0 }; self.n as usize + 4];

        // set edge, vertex, rectangle and wall "0"
        self.edges.init(EdgeDir::None, 0, 0, 0, 0, 0, 0, 0);
        self.vertices.init(0, 0, 0, 0);
        self.rectangles.init(0, 0, 0, 0);
        self.walls.init(0, 0);

        // set edge properties on the top side of the rectangulation
        for i in 1..self.n {
            self.edges[i].init(EdgeDir::Hor, i, i + 1, i - 1, i + 1, 0, i, 1);
        }
        self.edges[self.n].init(EdgeDir::Hor, self.n, self.n + 1, self.n - 1, 0, 0, self.n, 1);

        // set vertex properties on the top side of the rectangulation
        self.vertices.init(0, 1, 2 * self.n + 1, 0);
        self.vertices[self.n + 1].init(0, 0, 3 * self.n + 1, self.n);
        for i in 2..self.n {
            self.vertices[i].init(0, i, 2 * self.n + i, i - 1);
        }

        // set edge properties on the bottom side of the rectangulation
        for i in 2..self.n {
            self.edges[self.n + i].init(EdgeDir::Hor, self.n + i + 1, self.n + i + 2, self.n + i - 1, i, 0, 2);
        }
        self.edges[2 * self.n].init(EdgeDir::Hor, 2 * self.n + 1, 2 * self.n + 2, 2 * self.n - 1, 0, self.n, 2);

        // set vertex properties on the bottom side of the rectangulation
        self.vertices[self.n + 2].init(2 * self.n + 1, self.n + 1, 0, 0);
        self.vertices[2 * self.n + 2].init(3 * self.n + 1, 0, 0, 2 * self.n);
        for i in 2..self.n {
            self.vertices[self.n + i + 1].init(2 * self.n + i, self.n + i, 0, self.n + i - 1);
        }

        // set edge properties of the vertical edges
        self.edges[2 * self.n + 1].init(EdgeDir::Ver, self.n + 2, 1, 0, 0, 0, 1, 3);
        self.edges[3 * self.n + 1].init(EdgeDir::Ver, 2 * self.n + 2, self.n + 1, 0, 0, self.n, 0, self.n + 3);

        for i in 2..self.n {
            self.edges[2 * self.n + i].init(EdgeDir::Ver, self.n + i + 1, i, 0, 0, i - 1, i, i + 2);
        }

        // set rectangle properties
        for i in 1..self.n {
            self.rectangles[i].init(i + 1, self.n + i + 2, self.n + i + 1, i);
        }

        // set wall parameters
        self.walls.init(1, self.n + 1);
        self.walls.init(self.n + 2, 2 * self.n + 2);
        for i in 1..self.n + 1 {
            self.walls[i + 2].init(self.n + i + 1, i);
        }
    }

    /// Prints the data structures of the rectangulation (debug method).
    pub fn print_data(&self) {
        println!("edges:");
        for (i, e) in self.edges.iter().enumerate() {
            println!("  {}. {}", i, match e.dir {
                EdgeDir::Hor => "Hor",
                EdgeDir::Ver => "Ver",
                EdgeDir::None => "None",
            });
        }

        println!("vertices:");
        for (i, v) in self.vertices.iter().enumerate() {
            println!("  {}. {}", i, match v.type_ {
                VertexType::None => "None",
                VertexType::Corner => "Corner",
                VertexType::Bottom => "Bottom",
                VertexType::Top => "Top",
                VertexType::Left => "Left",
                VertexType::Right => "Right",
            });
        }

        println!("walls:");
        for (i, w) in self.walls.iter().enumerate() {
            println!("  {}. {}", i, w.first);
        }

        println!("rectangles:");
        for (i, r) in self.rectangles.iter().enumerate() {
            println!("  {}. {}", i, r.neast);
        }
    }

    /// Prints coordinates based on the rectangulation type.
    pub fn print_coordinates(&self) {
        match self.type_ {
            RectangulationType::Generic => self.print_coordinates_generic(),
            RectangulationType::Diagonal => self.print_coordinates_diagonal(),
            RectangulationType::Baligned => self.print_coordinates_diagonal(),
        }
    }

    /// Prints coordinates for a generic rectangulation.
    ///
    /// Uses a greedy algorithm to find an equispaced grid for vertex placement.
    pub fn print_coordinates_generic(&self) {
        // Run the greedy algorithm for finding an equispaced grid to place the vertices of the rectangles.
        // For finding the x-coordinates of the grid, we do a sweep from west to east.
        let mut active_vertices = vec[];
        let mut vertex_x_coord = vec![0; 2 * self.n as usize + 3];
        let mut x_value = 0;

        // start with every vertex which lies on the western side as an active vertex
        for (i, v) in self.vertices.iter().enumerate() {
            if v.type_ == VertexType::Right {
                active_vertices.push(i as i32);
            }
        }

        while !active_vertices.is_empty() {
            // set x-coordinate for active vertices
            for &idx in &active_vertices {
                vertex_x_coord[idx as usize] = x_value;
            }
            x_value += 1;

            let mut new_active_vertices = vec[];
            // propagate east
            for &idx in &active_vertices {
                if self.vertices[idx as usize].east != 0 {
                    new_active_vertices.push(self.edges[self.vertices[idx as usize].east as usize].head);
                }
            }

            // propagate north and south
            let mut new_active_vertices_copy = new_active_vertices.clone();
            for &idx in &new_active_vertices_copy {
                let mut propagate_from = self.vertices[idx as usize].clone();
                // propagate north
                while propagate_from.north != 0 {
                    let e = self.edges[propagate_from.north as usize].clone();
                    propagate_from = self.vertices[e.head as usize].clone();
                    new_active_vertices.push(e.head);
                }
                propagate_from = self.vertices[idx as usize].clone();
                // propagate south
                while propagate_from.south != 0 {
                    let e = self.edges[propagate_from.south as usize].clone();
                    propagate_from = self.vertices[e.tail as usize].clone();
                    new_active_vertices.push(e.tail);
                }
            }

            active_vertices = new_active_vertices;
        }

        // For finding the y-coordinates of the grid, we do the same from south to north.
        active_vertices.clear();
        let mut vertex_y_coord = vec![0; 2 * self.n as usize + 3];
        let mut y_value = 0;

        // start with every vertex which lies on the southern side as an active vertex
        for (i, v) in self.vertices.iter().enumerate() {
            if v.type_ == VertexType::Top {
                active_vertices.push(i as i32);
            }
        }

        while !active_vertices.is_empty() {
            // set y-coordinate for active vertices
            for &idx in &active_vertices {
                vertex_y_coord[idx as usize] = y_value;
            }
            y_value += 1;

            let mut new_active_vertices = vec[];
            // propagate north
            for &idx in &active_vertices {
                if self.vertices[idx as usize].north != 0 {
                    new_active_vertices.push(self.edges[self.vertices[idx as usize].north as usize].head);
                }
            }

            // propagate west and east
            let mut new_active_vertices_copy = new_active_vertices.clone();
            for &idx in &new_active_vertices_copy {
                let mut propagate_from = self.vertices[idx as usize].clone();
                // propagate east
                while propagate_from.east != 0 {
                    let e = self.edges[propagate_from.east as usize].clone();
                    propagate_from = self.vertices[e.head as usize].clone();
                    new_active_vertices.push(e.head);
                }
                propagate_from = self.vertices[idx as usize].clone();
                // propagate west
                while propagate_from.west != 0 {
                    let e = self.edges[propagate_from.west as usize].clone();
                    propagate_from = self.vertices[e.tail as usize].clone();
                    new_active_vertices.push(e.tail);
                }
            }

            active_vertices = new_active_vertices;
        }

        let mut is_first = true;
        let mut is_second = true;
        for r in &self.rectangles {
            if is_first {
                is_first = false;
                continue; // skip non-used 0th rectangle
            }
            if !is_second {
                print!(" | "); // don't print separator before first rectangle
            } else {
                is_second = false;
            }
            println!("{} {}", vertex_x_coord[r.swest as usize]);
        }
    }

    /// Prints coordinates for a diagonal rectangulation.
    ///
    /// This is obtained by mirroring the generic case along the main diagonal.
    pub fn print_coordinates_diagonal(&self) {
        // this code is obtained by mirroring the case RectangulationType::Generic along the main diagonal
    }

    /// Advances to the next rectangulation.
    ///
    /// Runs one iteration of the memoryless algorithm.
    /// Returns true if the next rectangulation is not the identity, false otherwise.
    pub fn next(&mut self) -> bool {
        // Run one iteration of the memoryless algorithm, return true if the next rectangulation is not the identity, false otherwise.
        // M3 - Select rectangle
        let j = self.s[self.n as usize];
        if (j == 1) || (self.n == 2 && self.type_ == RectangulationType::Baligned) {
            return false;
        }

        // M4 - Jump rectangle
        match self.type_ {
            RectangulationType::Generic => self.next_generic(j, self.o[j as usize]),
            RectangulationType::Diagonal => self.next_diagonal(j, self.o[j as usize]),
            RectangulationType::Baligned => self.next_baligned(j, self.o[j as usize]),
        }

        // M5 - Update o and s
        self.s[self.n as usize] = self.n;
        if (self.type_ == RectangulationType::Baligned && self.o[j - 1] == RectangulationDirection::Left && self.is_bottom_based(j - 1)) {
            self.o[j - 1] = RectangulationDirection::Right;
            self.s[j - 1] = self.s[j - 2];
            self.s[j - 2] = j - 2;
        }
        if (self.o[j] == RectangulationDirection::Left && self.is_bottom_based(j)) {
            self.o[j] = RectangulationDirection::Right;
            self.s[j] = self.s[j - 1];
            self.s[j - 1] = j - 1;
        }
        if (self.type_ == RectangulationType::Baligned && self.o[j - 1] == RectangulationDirection::Right && self.is_right_based(j - 1)) {
            self.o[j - 1] = RectangulationDirection::Left;
            self.s[j - 1] = self.s[j - 2];
            self.s[j - 2] = j - 2;
        }
        if (self.o[j] == RectangulationDirection::Right && self.is_right_based(j)) {
            self.o[j] = RectangulationDirection::Left;
            self.s[j] = self.s[j - 1];
            self.s[j - 1] = j - 1;
        }

        true
    }

    /// Checks if rectangle j is bottom-based.
    ///
    /// # Arguments
    ///
    /// * `j` - The rectangle index.
    pub fn is_bottom_based(&self, j: i32) -> bool {
        // this code is obtained by mirroring the case RectangulationDirection::Right along the main diagonal
        let a = self.rectangles[j as usize].nwest;
        let alpha = self.vertices[a as usize].west;
        let beta = self.edges[alpha as usize].tail;
        let k = self.edges[beta as usize].left;
        let l = self.edges[alpha as usize].right;
        let m = self.edges[beta as usize].right;
        if self.vertices[beta as usize].type_ == VertexType::top {
            true
        } else {
            false
        }
    }

    /// Checks if rectangle j is right-based.
    ///
    /// # Arguments
    ///
    /// * `j` - The rectangle index.
    pub fn is_right_based(&self, j: i32) -> bool {
        // this code is obtained by mirroring the case RectangulationDirection::Left along the main diagonal
        let a = self.rectangles[j as usize].nwest;
        let alpha = self.vertices[a as usize].east;
        let beta = self.edges[alpha as usize].head;
        let k = self.edges[beta as usize].left;
        let l = self.edges[alpha as usize].right;
        let m = self.edges[beta as usize].right;
        if self.vertices[beta as usize].type_ == VertexType::Left {
            true
        } else {
            false
        }
    }

    /// Removes the head edge from a wall.
    ///
    /// # Arguments
    ///
    /// * `beta` - The edge index to remove.
    pub fn rem_head(&mut self, beta: i32) {
        // Step 1 - Prepare
        let alpha = self.edges[beta as usize].prev;
        let a = self.edges[alpha as usize].tail;
        // Step 2 - Update edges/vertices
        if alpha != 0 {
            self.edges[alpha as usize].next = beta;
        }
        self.edges[beta as usize].prev = alpha;
        self.vertices[a as usize].east = 0;
        self.vertices[a as usize].type_ = VertexType::Right;
        // Step 3 - Update wall
        let x = self.edges[beta as usize].wall;
        self.walls[x as usize].first = a;
    }

    /// Removes the tail edge from a wall.
    ///
    /// # Arguments
    ///
    /// * `beta` - The edge index to remove.
    pub fn rem_tail(&mut self, beta: i32) {
        // Step 1 - Prepare
        let alpha = self.edges[beta as usize].next;
        let a = self.edges[alpha as usize].head;
        // Step 2 - Update edges/vertices
        if alpha != 0 {
            self.edges[alpha as usize].prev = beta;
        }
        self.edges[beta as usize].next = alpha;
        self.vertices[a as usize].west = 0;
        self.vertices[a as usize].type_ = VertexType::Left;
        // Step 3 - Update wall
        let x = self.edges[beta as usize].wall;
        self.walls[x as usize].last = a;
    }
}
        
impl Rectangulation {
    fn init(&mut self, vertices: Vec<Vertex>, walls: Vec<Wall>, edges: Vec<Edge>, rectangles: Vec<Rectangle>) {
        self.edges = edges;
        self.vertices = vertices;
        self.walls = walls;
        self.rectangles = rectangles;
    }

    fn set_all_vertical(&mut self) {
        // initialize the 3n+1 edges, 2n+2 vertices, n rectangles and n+3 walls
        let mut edges: Vec<Edge> = vec![Edge::default(); 3 * self.n + 2];
        let mut vertices: Vec<Vertex> = vec![Vertex::default(); 2 * self.n + 3];
        let mut rectangles: Vec<Rectangle> = vec![Rectangle::default(); self.n + 1];
        let mut walls: Vec<Wall> = vec![Wall::default(); self.n + 4];
        // set edge, vertex, rectangle and wall "0"
        edges[0] = Edge { /* initialize with values */ };
        vertices[0] = Vertex { /* initialize with values */ };
        rectangles[0] = Rectangle { /* initialize with values */ };
        walls[0] = Wall { /* initialize with values */ };
        // set edge properties on the top side of the rectangulation
        for i in 1..self.n {
            edges[i] = Edge { /* initialize with values */ };
        }
        edges[self.n] = Edge { /* initialize with values */ };
        // set vertex properties on the top side of the rectangulation
        vertices[1] = Vertex { /* initialize with values */ };
        vertices[self.n + 1] = Vertex { /* initialize with values */ };
        for i in 2..=self.n {
            vertices[i] = Vertex { /* initialize with values */ };
        }
        // set edge properties on the bottom side of the rectangulation
        for i in 2..self.n {
            edges[self.n + i] = Edge { /* initialize with values */ };
        }
        edges[self.n + 1] = Edge { /* initialize with values */ };
        edges[2 * self.n] = Edge { /* initialize with values */ };
        // set vertex properties on the bottom side of the rectangulation
        vertices[self.n + 2] = Vertex { /* initialize with values */ };
        vertices[2 * self.n + 2] = Vertex { /* initialize with values */ };
        for i in 2..=self.n {
            vertices[self.n + i + 1] = Vertex { /* initialize with values */ };
        }
        // set edge properties of the vertical edges
        edges[2 * self.n + 1] = Edge { /* initialize with values */ };
        edges[3 * self.n + 1] = Edge { /* initialize with values */ };
        for i in 2..=self.n {
            edges[2 * self.n + i] = Edge { /* initialize with values */ };
        }
        // set rectangle properties
        for i in 1..=self.n {
            rectangles[i] = Rectangle { /* initialize with values */ };
        }
        // set wall parameters
        walls[1] = Wall { /* initialize with values */ };
        walls[2] = Wall { /* initialize with values */ };
        for i in 1..=self.n + 1 {
            walls[i + 2] = Wall { /* initialize with values */ };
        }
        self.init(vertices, walls, edges, rectangles);
    }
}

