use std::collections::VecDeque;
use std::fmt;
use std::str::FromStr;
use clap::{Parser, ValueEnum};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum EdgeDir {
    Hor,
    Ver,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum VertexType {
    Corner,
    Top,
    Bottom,
    Left,
    Right,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum RectangulationType {
    Generic,
    Baligned,
    Diagonal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum RectangulationDirection {
    Left,
    Right,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum RectangulationPattern {
    WmillClockwise,
    WmillCounterclockwise,
    BrickLeftRight,
    BrickRightLeft,
    BrickBottomTop,
    BrickTopBottom,
    HVertical,
    HHorizontal,
}

#[derive(Debug, Clone)]
struct Edge {
    dir: EdgeDir,
    tail: usize,
    head: usize,
    prev: usize,
    next: usize,
    left: usize,
    right: usize,
    wall: usize,
}

impl Edge {
    fn new() -> Self {
        Self {
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

    fn init(
        &mut self,
        dir: EdgeDir,
        tail: usize,
        head: usize,
        prev: usize,
        next: usize,
        left: usize,
        right: usize,
        wall: usize,
    ) {
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

#[derive(Debug, Clone)]
struct Vertex {
    north: usize,
    east: usize,
    south: usize,
    west: usize,
    type_: VertexType,
}

impl Vertex {
    fn new() -> Self {
        Self {
            north: 0,
            east: 0,
            south: 0,
            west: 0,
            type_: VertexType::None,
        }
    }

    fn init(&mut self, north: usize, east: usize, south: usize, west: usize) {
        self.north = north;
        self.east = east;
        self.south = south;
        self.west = west;

        let zeros = (self.north == 0) as u8
            + (self.south == 0) as u8
            + (self.west == 0) as u8
            + (self.east == 0) as u8;

        self.type_ = if zeros >= 3 || zeros == 0 {
            VertexType::None
        } else if zeros == 2 {
            VertexType::Corner
        } else if self.south == 0 {
            VertexType::Top
        } else if self.north == 0 {
            VertexType::Bottom
        } else if self.east == 0 {
            VertexType::Left
        } else {
            VertexType::Right
        };
    }
}

#[derive(Debug, Clone)]
struct Wall {
    first: usize,
    last: usize,
}

impl Wall {
    fn new() -> Self {
        Self { first: 0, last: 0 }
    }

    fn init(&mut self, first: usize, last: usize) {
        self.first = first;
        self.last = last;
    }
}

#[derive(Debug, Clone)]
struct Rectangle {
    nwest: usize,
    neast: usize,
    swest: usize,
    seast: usize,
}

impl Rectangle {
    fn new() -> Self {
        Self {
            nwest: 0,
            neast: 0,
            swest: 0,
            seast: 0,
        }
    }

    fn init(&mut self, neast: usize, seast: usize, swest: usize, nwest: usize) {
        self.nwest = nwest;
        self.neast = neast;
        self.swest = swest;
        self.seast = seast;
    }
}

#[derive(Debug, Clone)]
struct Rectangulation {
    n: usize,
    type_: RectangulationType,
    patterns: Vec<RectangulationPattern>,
    o: Vec<RectangulationDirection>,
    s: Vec<usize>,
    vertices: Vec<Vertex>,
    walls: Vec<Wall>,
    edges: Vec<Edge>,
    rectangles: Vec<Rectangle>,
}

impl Rectangulation {
    fn new(n: usize, type_: RectangulationType, patterns: Vec<RectangulationPattern>) -> Self {
        let mut r = Self {
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
        r.set_all_vertical();
        r.o.push(RectangulationDirection::None);
        r.s.push(0); // Using 0 instead of -1 since usize can't be negative
        for j in 1..=n {
            r.o.push(RectangulationDirection::Left);
            r.s.push(j);
        }
        r
    }

    fn set_all_vertical(&mut self) {
        // Initialize data structures
        self.edges = vec![Edge::new(); 3 * self.n + 2];
        self.vertices = vec![Vertex::new(); 2 * self.n + 3];
        self.rectangles = vec![Rectangle::new(); self.n + 1];
        self.walls = vec![Wall::new(); self.n + 4];

        // Initialize edge, vertex, rectangle and wall "0"
        self.edges[0].init(EdgeDir::None, 0, 0, 0, 0, 0, 0, 0);
        self.vertices[0].init(0, 0, 0, 0);
        self.rectangles[0].init(0, 0, 0, 0);
        self.walls[0].init(0, 0);

        // Set edge properties on the top side
        for i in 1..self.n {
            self.edges[i].init(
                EdgeDir::Hor,
                i,
                i + 1,
                i - 1,
                i + 1,
                0,
                i,
                1,
            );
        }
        self.edges[self.n].init(
            EdgeDir::Hor,
            self.n,
            self.n + 1,
            self.n - 1,
            0,
            0,
            self.n,
            1,
        );

        // Set vertex properties on the top side
        self.vertices[1].init(0, 1, 2 * self.n + 1, 0);
        self.vertices[self.n + 1].init(0, 0, 3 * self.n + 1, self.n);
        for i in 2..=self.n {
            self.vertices[i].init(0, i, 2 * self.n + i, i - 1);
        }

        // Set edge properties on the bottom side
        for i in 2..self.n {
            self.edges[self.n + i].init(
                EdgeDir::Hor,
                self.n + i + 1,
                self.n + i + 2,
                self.n + i - 1,
                self.n + i + 1,
                i,
                0,
                2,
            );
        }
        self.edges[self.n + 1].init(
            EdgeDir::Hor,
            self.n + 2,
            self.n + 3,
            0,
            self.n + 2,
            1,
            0,
            2,
        );
        self.edges[2 * self.n].init(
            EdgeDir::Hor,
            2 * self.n + 1,
            2 * self.n + 2,
            2 * self.n - 1,
            0,
            self.n,
            0,
            2,
        );

        // Set vertex properties on the bottom side
        self.vertices[self.n + 2].init(2 * self.n + 1, self.n + 1, 0, 0);
        self.vertices[2 * self.n + 2].init(3 * self.n + 1, 0, 0, 2 * self.n);
        for i in 2..=self.n {
            self.vertices[self.n + i + 1].init(2 * self.n + i, self.n + i, 0, self.n + i - 1);
        }

        // Set edge properties of the vertical edges
        self.edges[2 * self.n + 1].init(
            EdgeDir::Ver,
            self.n + 2,
            1,
            0,
            0,
            0,
            1,
            3,
        );
        self.edges[3 * self.n + 1].init(
            EdgeDir::Ver,
            2 * self.n + 2,
            self.n + 1,
            0,
            0,
            self.n,
            0,
            self.n + 3,
        );
        for i in 2..=self.n {
            self.edges[2 * self.n + i].init(
                EdgeDir::Ver,
                self.n + i + 1,
                i,
                0,
                0,
                i - 1,
                i,
                i + 2,
            );
        }

        // Set rectangle properties
        for i in 1..=self.n {
            self.rectangles[i].init(i + 1, self.n + i + 2, self.n + i + 1, i);
        }

        // Set wall parameters
        self.walls[1].init(1, self.n + 1);
        self.walls[2].init(self.n + 2, 2 * self.n + 2);
        for i in 1..=self.n + 1 {
            self.walls[i + 2].init(self.n + i + 1, i);
        }
    }

    fn print_data(&self) {
        println!("edges:");
        for (i, e) in self.edges.iter().enumerate() {
            let dir_str = match e.dir {
                EdgeDir::Hor => "Hor",
                EdgeDir::Ver => "Ver",
                EdgeDir::None => "None",
            };
            println!(
                "\t{}. {} {} {} {} {} {} {}",
                i, dir_str, e.tail, e.head, e.prev, e.next, e.left, e.right, e.wall
            );
        }

        println!("vertices:");
        for (i, v) in self.vertices.iter().enumerate() {
            let type_str = match v.type_ {
                VertexType::None => "None",
                VertexType::Corner => "corner",
                VertexType::Bottom => "bottom",
                VertexType::Top => "top",
                VertexType::Left => "left",
                VertexType::Right => "right",
            };
            println!(
                "\t{}. {} {} {} {} {}",
                i, v.north, v.east, v.south, v.west, type_str
            );
        }

        println!("walls:");
        for (i, w) in self.walls.iter().enumerate() {
            println!("\t{}. {} {}", i, w.first, w.last);
        }

        println!("rectangles:");
        for (i, r) in self.rectangles.iter().enumerate() {
            println!(
                "\t{}. {} {} {} {}",
                i, r.neast, r.seast, r.swest, r.nwest
            );
        }
    }

    fn print_coordinates(&self) {
        match self.type_ {
            RectangulationType::Generic => self.print_coordinates_generic(),
            RectangulationType::Diagonal | RectangulationType::Baligned => {
                self.print_coordinates_diagonal()
            }
        }
    }

    fn print_coordinates_generic(&self) {
        let mut active_vertices = Vec::new();
        let mut vertex_x_coord = vec![-1; 2 * self.n + 3];

        // Find initial active vertices (western side)
        for (a, v) in self.vertices.iter().enumerate() {
            let side_edge_id = match v.type_ {
                VertexType::Right => v.north,
                VertexType::Corner => v.north.max(v.south),
                _ => continue,
            };

            if self.edges[side_edge_id].left == 0 {
                active_vertices.push(a);
            }
        }

        // Propagate x-coordinates
        let mut x_value = 0;
        while !active_vertices.is_empty() {
            for &idx in &active_vertices {
                vertex_x_coord[idx] = x_value;
            }
            x_value += 1;

            let mut new_active_vertices = Vec::new();
            for &idx in &active_vertices {
                if self.vertices[idx].east != 0 {
                    let alpha = self.vertices[idx].east;
                    new_active_vertices.push(self.edges[alpha].head);
                }
            }

            // Propagate north and south
            let new_active_vertices_copy = new_active_vertices.clone();
            for &idx in &new_active_vertices_copy {
                let mut propagate_from = &self.vertices[idx];
                while propagate_from.north != 0 {
                    let e = &self.edges[propagate_from.north];
                    propagate_from = &self.vertices[e.head];
                    new_active_vertices.push(e.head);
                }

                let mut propagate_from = &self.vertices[idx];
                while propagate_from.south != 0 {
                    let e = &self.edges[propagate_from.south];
                    propagate_from = &self.vertices[e.tail];
                    new_active_vertices.push(e.tail);
                }
            }

            active_vertices = new_active_vertices;
        }

        // Find y-coordinates (southern side)
        let mut active_vertices = Vec::new();
        let mut vertex_y_coord = vec![-1; 2 * self.n + 3];

        for (a, v) in self.vertices.iter().enumerate() {
            let side_edge_id = match v.type_ {
                VertexType::Top => v.east,
                VertexType::Corner => v.east.max(v.west),
                _ => continue,
            };

            if self.edges[side_edge_id].right == 0 {
                active_vertices.push(a);
            }
        }

        // Propagate y-coordinates
        let mut y_value = 0;
        while !active_vertices.is_empty() {
            for &idx in &active_vertices {
                vertex_y_coord[idx] = y_value;
            }
            y_value += 1;

            let mut new_active_vertices = Vec::new();
            for &idx in &active_vertices {
                if self.vertices[idx].north != 0 {
                    let alpha = self.vertices[idx].north;
                    new_active_vertices.push(self.edges[alpha].head);
                }
            }

            // Propagate east and west
            let new_active_vertices_copy = new_active_vertices.clone();
            for &idx in &new_active_vertices_copy {
                let mut propagate_from = &self.vertices[idx];
                while propagate_from.east != 0 {
                    let e = &self.edges[propagate_from.east];
                    propagate_from = &self.vertices[e.head];
                    new_active_vertices.push(e.head);
                }

                let mut propagate_from = &self.vertices[idx];
                while propagate_from.west != 0 {
                    let e = &self.edges[propagate_from.west];
                    propagate_from = &self.vertices[e.tail];
                    new_active_vertices.push(e.tail);
                }
            }

            active_vertices = new_active_vertices;
        }

        // Print rectangle coordinates
        let mut output = Vec::new();
        for r in self.rectangles.iter().skip(1) {
            let swest_x = vertex_x_coord[r.swest];
            let swest_y = vertex_y_coord[r.swest];
            let neast_x = vertex_x_coord[r.neast];
            let neast_y = vertex_y_coord[r.neast];
            output.push(format!("{} {} {} {}", swest_x, swest_y, neast_x, neast_y));
        }
        println!("{}", output.join(" | "));
    }

    fn dfs_bl(
        &self,
        vertex_id: usize,
        val: &mut i32,
        vertex_x_coord: &mut [i32],
        vertex_y_coord: &mut [i32],
    ) {
        let v = &self.vertices[vertex_id];
        let top_vertex = self.edges[v.north].head;
        let right_vertex = self.edges[v.east].head;

        if matches!(
            self.vertices[top_vertex].type_,
            VertexType::Corner | VertexType::Bottom | VertexType::Left
        ) {
            vertex_x_coord[vertex_id] = *val;
            *val += 1;
        } else {
            self.dfs_bl(top_vertex, val, vertex_x_coord, vertex_y_coord);
            vertex_x_coord[vertex_id] = vertex_x_coord[top_vertex];
        }

        if matches!(
            self.vertices[right_vertex].type_,
            VertexType::Corner | VertexType::Bottom | VertexType::Left
        ) {
            vertex_y_coord[vertex_id] = self.n as i32 - *val;
            *val += 1;
        } else {
            self.dfs_bl(right_vertex, val, vertex_x_coord, vertex_y_coord);
            vertex_y_coord[vertex_id] = vertex_y_coord[right_vertex];
        }
    }

    fn dfs_tr(
        &self,
        vertex_id: usize,
        val: &mut i32,
        vertex_x_coord: &mut [i32],
        vertex_y_coord: &mut [i32],
    ) {
        let v = &self.vertices[vertex_id];
        let left_vertex = self.edges[v.west].tail;
        let bottom_vertex = self.edges[v.south].tail;

        if matches!(
            self.vertices[bottom_vertex].type_,
            VertexType::Corner | VertexType::Right | VertexType::Top
        ) {
            vertex_x_coord[vertex_id] = self.n as i32 - *val;
            *val += 1;
        } else {
            self.dfs_tr(bottom_vertex, val, vertex_x_coord, vertex_y_coord);
            vertex_x_coord[vertex_id] = vertex_x_coord[bottom_vertex];
        }

        if matches!(
            self.vertices[left_vertex].type_,
            VertexType::Corner | VertexType::Right | VertexType::Top
        ) {
            vertex_y_coord[vertex_id] = *val;
            *val += 1;
        } else {
            self.dfs_tr(left_vertex, val, vertex_x_coord, vertex_y_coord);
            vertex_y_coord[vertex_id] = vertex_y_coord[left_vertex];
        }
    }

    fn print_coordinates_diagonal(&self) {
        let mut vertex_x_coord = vec![-1; 2 * self.n + 3];
        let mut vertex_y_coord = vec![-1; 2 * self.n + 3];

        // Find bottom-left and top-right corners
        let mut bl = None;
        let mut tr = None;
        for i in 1..2 * self.n + 3 {
            let v = &self.vertices[i];
            if v.north == 0 && v.east == 0 && v.type_ == VertexType::Corner {
                tr = Some(i);
            } else if v.south == 0 && v.west == 0 && v.type_ == VertexType::Corner {
                bl = Some(i);
            }
        }

        let bl = bl.expect("Bottom-left corner not found");
        let tr = tr.expect("Top-right corner not found");

        let mut val = 0;
        self.dfs_bl(bl, &mut val, &mut vertex_x_coord, &mut vertex_y_coord);
        let mut val = 0;
        self.dfs_tr(tr, &mut val, &mut vertex_x_coord, &mut vertex_y_coord);

        // Print rectangle coordinates
        let mut output = Vec::new();
        for r in self.rectangles.iter().skip(1) {
            let swest_x = vertex_x_coord[r.swest];
            let swest_y = vertex_y_coord[r.swest];
            let neast_x = vertex_x_coord[r.neast];
            let neast_y = vertex_y_coord[r.neast];
            output.push(format!("{} {} {} {}", swest_x, swest_y, neast_x, neast_y));
        }
        println!("{}", output.join(" | "));
    }

    fn next(&mut self) -> bool {
        let j = self.s[self.n];
        if j == 1 || (self.n == 2 && self.type_ == RectangulationType::Baligned) {
            return false;
        }

        match self.type_ {
            RectangulationType::Generic => {
                self.next_generic(j, self.o[j]);
                while self.contains_pattern(j) {
                    self.next_generic(j, self.o[j]);
                }
            }
            RectangulationType::Diagonal => {
                self.next_diagonal(j, self.o[j]);
                while self.contains_pattern(j) {
                    self.next_diagonal(j, self.o[j]);
                }
            }
            RectangulationType::Baligned => {
                self.next_baligned(j, self.o[j]);
                while self.contains_pattern(j) {
                    self.next_baligned(j, self.o[j]);
                }
            }
        }

        self.s[self.n] = self.n;

        if self.type_ == RectangulationType::Baligned
            && self.o[j - 1] == RectangulationDirection::Left
            && self.is_bottom_based(j - 1)
        {
            self.o[j - 1] = RectangulationDirection::Right;
            self.s[j - 1] = self.s[j - 2];
            self.s[j - 2] = j - 2;
        }

        if self.o[j] == RectangulationDirection::Left && self.is_bottom_based(j) {
            self.o[j] = RectangulationDirection::Right;
            self.s[j] = self.s[j - 1];
            self.s[j - 1] = j - 1;
        }

        if self.type_ == RectangulationType::Baligned
            && self.o[j - 1] == RectangulationDirection::Right
            && self.is_right_based(j - 1)
        {
            self.o[j - 1] = RectangulationDirection::Left;
            self.s[j - 1] = self.s[j - 2];
            self.s[j - 2] = j - 2;
        }

        if self.o[j] == RectangulationDirection::Right && self.is_right_based(j) {
            self.o[j] = RectangulationDirection::Left;
            self.s[j] = self.s[j - 1];
            self.s[j - 1] = j - 1;
        }

        true
    }

    fn is_bottom_based(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        let alpha = self.vertices[a].south;
        let b = self.rectangles[j].swest;

        if self.edges[alpha].left == 0 {
            return true;
        }

        if self.type_ == RectangulationType::Baligned
            && a == self.rectangles[j - 1].neast
            && b == self.rectangles[j - 1].seast
        {
            let c = self.rectangles[j - 1].nwest;
            let gamma = self.vertices[c].south;
            if self.edges[gamma].left == 0 {
                return true;
            }
        }

        false
    }

    fn is_right_based(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        let alpha = self.vertices[a].east;
        let b = self.rectangles[j].neast;

        if self.edges[alpha].left == 0 {
            return true;
        }

        if self.type_ == RectangulationType::Baligned
            && a == self.rectangles[j - 1].swest
            && b == self.rectangles[j - 1].seast
        {
            let c = self.rectangles[j].nwest;
            let gamma = self.vertices[c].east;
            if self.edges[gamma].left == 0 {
                return true;
            }
        }

        false
    }

    fn rem_head(&mut self, beta: usize) {
        let alpha = self.edges[beta].prev;
        let gamma = self.edges[beta].next;
        let a = self.edges[beta].tail;

        if alpha != 0 {
            self.edges[alpha].next = gamma;
        }
        if gamma != 0 {
            self.edges[gamma].prev = alpha;
        }
        self.edges[gamma].tail = a;

        match self.edges[beta].dir {
            EdgeDir::Hor => self.vertices[a].east = gamma,
            EdgeDir::Ver => self.vertices[a].north = gamma,
            EdgeDir::None => {}
        }

        let x = self.edges[beta].wall;
        if self.edges[beta].head == self.walls[x].last {
            self.walls[x].last = a;
        }
    }

    fn rem_tail(&mut self, beta: usize) {
        let alpha = self.edges[beta].prev;
        let gamma = self.edges[beta].next;
        let a = self.edges[beta].head;

        if alpha != 0 {
            self.edges[alpha].next = gamma;
        }
        if gamma != 0 {
            self.edges[gamma].prev = alpha;
        }
        self.edges[alpha].head = a;

        match self.edges[beta].dir {
            EdgeDir::Hor => self.vertices[a].west = alpha,
            EdgeDir::Ver => self.vertices[a].south = alpha,
            EdgeDir::None => {}
        }

        let x = self.edges[beta].wall;
        if self.edges[beta].tail == self.walls[x].first {
            self.walls[x].first = a;
        }
    }

    fn ins_before(&mut self, beta: usize, a: usize, gamma: usize) {
        let alpha = self.edges[gamma].prev;
        let b = self.edges[gamma].tail;

        self.edges[beta].tail = b;
        self.edges[beta].head = a;
        self.edges[beta].prev = alpha;
        self.edges[beta].next = gamma;
        self.edges[gamma].tail = a;
        self.edges[gamma].prev = beta;

        if alpha != 0 {
            self.edges[alpha].next = beta;
        }

        match self.edges[gamma].dir {
            EdgeDir::Hor => {
                self.edges[beta].dir = EdgeDir::Hor;
                self.vertices[a].west = beta;
                self.vertices[a].east = gamma;
                self.vertices[b].east = beta;
            }
            EdgeDir::Ver => {
                self.edges[beta].dir = EdgeDir::Ver;
                self.vertices[a].south = beta;
                self.vertices[a].north = gamma;
                self.vertices[b].north = beta;
            }
            EdgeDir::None => {}
        }

        self.edges[beta].wall = self.edges[gamma].wall;
    }

    fn ins_after(&mut self, alpha: usize, a: usize, beta: usize) {
        let gamma = self.edges[alpha].next;
        let b = self.edges[alpha].head;

        self.edges[beta].tail = a;
        self.edges[beta].head = b;
        self.edges[beta].prev = alpha;
        self.edges[beta].next = gamma;
        self.edges[alpha].head = a;
        self.edges[alpha].next = beta;

        if gamma != 0 {
            self.edges[gamma].prev = beta;
        }

        match self.edges[alpha].dir {
            EdgeDir::Hor => {
                self.edges[beta].dir = EdgeDir::Hor;
                self.vertices[a].west = alpha;
                self.vertices[a].east = beta;
                self.vertices[b].west = beta;
            }
            EdgeDir::Ver => {
                self.edges[beta].dir = EdgeDir::Ver;
                self.vertices[a].south = alpha;
                self.vertices[a].north = beta;
                self.vertices[b].south = beta;
            }
            EdgeDir::None => {}
        }

        self.edges[beta].wall = self.edges[alpha].wall;
    }

    fn wjump_hor(&mut self, j: usize, dir: RectangulationDirection, alpha: usize) {
        match dir {
            RectangulationDirection::Left => {
                let a = self.rectangles[j].nwest;
                let beta = self.vertices[a].west;
                debug_assert!(alpha == self.edges[beta].prev);
                let k = self.edges[alpha].left;

                self.rem_head(beta);
                self.ins_after(alpha, a, beta);
                self.edges[beta].left = k;
                self.edges[beta].right = j;
            }
            RectangulationDirection::Right => {
                let a = self.rectangles[j].nwest;
                let beta = self.vertices[a].east;
                debug_assert!(alpha == self.edges[beta].next);
                let k = self.edges[alpha].left;
                let alpha_prime = self.edges[beta].prev;
                let l = self.edges[alpha_prime].right;

                self.rem_tail(beta);
                self.ins_before(beta, a, alpha);
                self.edges[beta].left = k;
                self.edges[beta].right = l;
            }
            RectangulationDirection::None => unreachable!(),
        }
    }

    fn wjump_ver(&mut self, j: usize, dir: RectangulationDirection, alpha: usize) {
        match dir {
            RectangulationDirection::Right => {
                let a = self.rectangles[j].nwest;
                let beta = self.vertices[a].north;
                debug_assert!(alpha == self.edges[beta].next);
                let k = self.edges[alpha].left;

                self.rem_tail(beta);
                self.ins_before(beta, a, alpha);
                self.edges[beta].left = k;
                self.edges[beta].right = j;
            }
            RectangulationDirection::Left => {
                let a = self.rectangles[j].nwest;
                let beta = self.vertices[a].south;
                debug_assert!(alpha == self.edges[beta].prev);
                let k = self.edges[alpha].left;
                let alpha_prime = self.edges[beta].next;
                let l = self.edges[alpha_prime].right;

                self.rem_head(beta);
                self.ins_after(alpha, a, beta);
                self.edges[beta].left = k;
                self.edges[beta].right = l;
            }
            RectangulationDirection::None => unreachable!(),
        }
    }

    fn sjump(&mut self, j: usize, d: RectangulationDirection, alpha: usize) {
        match d {
            RectangulationDirection::Left => {
                let a = self.rectangles[j].nwest;
                let b = self.rectangles[j].swest;
                let c = self.rectangles[j].neast;
                let alpha_prime = self.vertices[a].west;
                let beta = self.vertices[a].east;
                let beta_prime = self.vertices[b].west;
                let gamma = self.vertices[c].south;
                let delta = self.vertices[a].south;
                let c_prime = self.edges[beta_prime].tail;
                let k = self.edges[alpha].left;
                let l = self.edges[gamma].right;
                let x = self.edges[delta].wall;

                self.rem_tail(beta);
                self.rem_head(beta_prime);
                self.ins_before(beta, a, alpha);
                self.ins_after(gamma, b, beta_prime);

                self.edges[delta].dir = EdgeDir::Hor;
                self.edges[delta].tail = a;
                self.edges[delta].head = b;
                self.vertices[a].east = delta;
                self.vertices[a].west = 0;
                self.vertices[a].type_ = VertexType::Right;
                self.vertices[b].east = 0;
                self.vertices[b].west = delta;
                self.vertices[b].type_ = VertexType::Left;
                self.walls[x].first = a;
                self.walls[x].last = b;

                self.rectangles[j].neast = b;
                self.rectangles[j].swest = c_prime;
                self.rectangles[j - 1].neast = c;
                self.rectangles[j - 1].swest = a;

                let mut nu = self.vertices[c].west;
                while nu != alpha_prime {
                    self.edges[nu].right = j - 1;
                    nu = self.edges[nu].prev;
                }

                let mut nu = self.vertices[c_prime].north;
                while nu != alpha {
                    self.edges[nu].right = j;
                    nu = self.edges[nu].next;
                }

                self.edges[beta].left = k;
                self.edges[beta].right = j;
                self.edges[beta_prime].left = j - 1;
                self.edges[beta_prime].right = l;
            }
            RectangulationDirection::Right => {
                let a = self.rectangles[j].nwest;
                let b = self.rectangles[j].neast;
                let c = self.rectangles[j].swest;
                let alpha_prime = self.vertices[a].north;
                let beta = self.vertices[a].south;
                let beta_prime = self.vertices[b].north;
                let gamma = self.vertices[c].east;
                let delta = self.vertices[a].east;
                let c_prime = self.edges[beta_prime].head;
                let k = self.edges[alpha].left;
                let l = self.edges[gamma].right;
                let x = self.edges[delta].wall;

                self.rem_head(beta);
                self.rem_tail(beta_prime);
                self.ins_after(alpha, a, beta);
                self.ins_before(beta_prime, b, gamma);

                self.edges[delta].dir = EdgeDir::Ver;
                self.edges[delta].head = a;
                self.edges[delta].tail = b;
                self.vertices[a].south = delta;
                self.vertices[a].north = 0;
                self.vertices[a].type_ = VertexType::Bottom;
                self.vertices[b].south = 0;
                self.vertices[b].north = delta;
                self.vertices[b].type_ = VertexType::Top;
                self.walls[x].last = a;
                self.walls[x].first = b;

                self.rectangles[j].swest = b;
                self.rectangles[j].neast = c_prime;
                self.rectangles[j - 1].swest = c;
                self.rectangles[j - 1].neast = a;

                let mut nu = self.vertices[c].north;
                while nu != alpha_prime {
                    self.edges[nu].right = j - 1;
                    nu = self.edges[nu].next;
                }

                let mut nu = self.vertices[c_prime].west;
                while nu != alpha {
                    self.edges[nu].right = j;
                    nu = self.edges[nu].prev;
                }

                self.edges[beta].left = k;
                self.edges[beta].right = j;
                self.edges[beta_prime].left = j - 1;
                self.edges[beta_prime].right = l;
            }
            RectangulationDirection::None => unreachable!(),
        }
    }

    fn tjump_hor(&mut self, j: usize, dir: RectangulationDirection, alpha: usize) {
        match dir {
            RectangulationDirection::Left => {
                let a = self.rectangles[j].nwest;
                let b = self.edges[alpha].head;
                let c = self.rectangles[j].neast;
                let alpha_prime = self.vertices[a].west;
                let beta = self.vertices[a].east;
                let beta_prime = self.vertices[a].south;
                let gamma = self.vertices[c].south;
                let gamma_prime = self.vertices[b].south;
                let k = self.edges[beta_prime].left;
                let l = self.edges[gamma].right;
                let m = self.edges[alpha].right;
                let x = self.edges[alpha].wall;
                let y = self.edges[gamma_prime].wall;

                self.rem_tail(beta);
                self.rem_tail(beta_prime);
                self.ins_after(alpha, a, beta);
                self.ins_after(gamma, b, beta_prime);

                self.edges[beta].head = b;
                self.edges[gamma_prime].head = a;
                self.vertices[a].south = gamma_prime;
                self.vertices[b].west = beta;
                self.walls[x].last = b;
                self.walls[y].last = a;

                self.rectangles[j].neast = b;
                self.rectangles[k].neast = c;
                self.rectangles[m].neast = a;

                let mut nu = self.vertices[c].west;
                while nu != alpha_prime {
                    self.edges[nu].right = k;
                    nu = self.edges[nu].prev;
                }

                self.edges[beta].left = k;
                self.edges[beta_prime].right = l;
            }
            RectangulationDirection::Right => {
                let a = self.rectangles[j].nwest;
                let b = self.rectangles[j].neast;
                let alpha_prime = self.vertices[a].west;
                let gamma_prime = self.vertices[a].south;
                let beta = self.vertices[a].east;
                let beta_prime = self.vertices[b].north;
                let c = self.edges[beta_prime].head;
                let k = self.edges[beta].left;
                let l = self.edges[alpha].left;
                let m = self.edges[alpha_prime].right;
                let x = self.edges[alpha_prime].wall;
                let y = self.edges[gamma_prime].wall;

                self.rem_tail(beta);
                self.rem_tail(beta_prime);
                self.ins_after(alpha, a, beta);
                self.ins_after(gamma_prime, b, beta_prime);

                self.edges[alpha_prime].head = b;
                self.edges[beta_prime].head = a;
                self.vertices[a].south = beta_prime;
                self.vertices[b].west = alpha_prime;
                self.walls[x].last = b;
                self.walls[y].last = a;

                self.rectangles[j].neast = c;
                self.rectangles[k].neast = a;
                self.rectangles[m].neast = b;

                let mut nu = self.vertices[c].west;
                while nu != alpha {
                    self.edges[nu].right = j;
                    nu = self.edges[nu].prev;
                }

                self.edges[beta].left = l;
                self.edges[beta_prime].right = j;
            }
            RectangulationDirection::None => unreachable!(),
        }
    }

    fn tjump_ver(&mut self, j: usize, dir: RectangulationDirection, alpha: usize) {
        match dir {
            RectangulationDirection::Right => {
                let a = self.rectangles[j].nwest;
                let b = self.edges[alpha].tail;
                let c = self.rectangles[j].swest;
                let alpha_prime = self.vertices[a].north;
                let beta = self.vertices[a].south;
                let beta_prime = self.vertices[a].east;
                let gamma = self.vertices[c].east;
                let gamma_prime = self.vertices[b].east;
                let k = self.edges[beta_prime].left;
                let l = self.edges[gamma].right;
                let m = self.edges[alpha].right;
                let x = self.edges[alpha].wall;
                let y = self.edges[gamma_prime].wall;

                self.rem_head(beta);
                self.rem_head(beta_prime);
                self.ins_before(beta, a, alpha);
                self.ins_before(beta_prime, b, gamma);

                self.edges[beta].tail = b;
                self.edges[gamma_prime].tail = a;
                self.vertices[a].east = gamma_prime;
                self.vertices[a].type_ = VertexType::Right;
                self.vertices[b].north = beta;
                self.vertices[b].type_ = VertexType::Top;
                self.walls[x].first = b;
                self.walls[y].first = a;

                self.rectangles[j].swest = b;
                self.rectangles[k].swest = c;
                self.rectangles[m].swest = a;

                let mut nu = self.vertices[c].north;
                while nu != alpha_prime {
                    self.edges[nu].right = k;
                    nu = self.edges[nu].next;
                }

                self.edges[beta].left = k;
                self.edges[beta_prime].right = l;
            }
            RectangulationDirection::Left => {
                let a = self.rectangles[j].nwest;
                let b = self.rectangles[j].swest;
                let alpha_prime = self.vertices[a].north;
                let gamma_prime = self.vertices[a].east;
                let beta = self.vertices[a].south;
                let beta_prime = self.vertices[b].west;
                let c = self.edges[beta_prime].tail;
                let k = self.edges[beta].left;
                let l = self.edges[alpha].left;
                let m = self.edges[alpha_prime].right;
                let x = self.edges[alpha_prime].wall;
                let y = self.edges[gamma_prime].wall;

                self.rem_head(beta);
                self.rem_head(beta_prime);
                self.ins_before(beta, a, alpha);
                self.ins_before(beta_prime, b, gamma_prime);

                self.edges[alpha_prime].tail = b;
                self.edges[beta_prime].tail = a;
                self.vertices[a].east = beta_prime;
                self.vertices[a].type_ = VertexType::Right;
                self.vertices[b].north = alpha_prime;
                self.vertices[b].type_ = VertexType::Top;
                self.walls[x].first = b;
                self.walls[y].first = a;

                self.rectangles[j].swest = c;
                self.rectangles[k].swest = a;
                self.rectangles[m].swest = b;

                let mut nu = self.vertices[c].north;
                while nu != alpha {
                    self.edges[nu].right = j;
                    nu = self.edges[nu].next;
                }

                self.edges[beta].left = l;
                self.edges[beta_prime].right = j;
            }
            RectangulationDirection::None => unreachable!(),
        }
    }

    fn next_generic(&mut self, j: usize, dir: RectangulationDirection) {
        let a = self.rectangles[j].nwest;

        if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].west;
            let beta = self.vertices[a].south;
            let b = self.edges[beta].tail;
            let c = self.edges[alpha].tail;

            if self.vertices[c].type_ == VertexType::Top {
                let gamma = self.vertices[c].west;
                self.wjump_hor(j, RectangulationDirection::Left, gamma);
            } else if self.vertices[b].type_ == VertexType::Left {
                let gamma = self.vertices[b].west;
                self.tjump_hor(j, RectangulationDirection::Left, gamma);
            } else {
                let gamma = self.vertices[c].south;
                self.sjump(j, RectangulationDirection::Left, gamma);
            }
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].east;
            let b = self.edges[alpha].head;

            if self.vertices[b].type_ == VertexType::Top {
                let gamma = self.vertices[b].east;
                self.wjump_hor(j, RectangulationDirection::Right, gamma);
            } else {
                let k = self.edges[alpha].left;
                let c = self.rectangles[k].nwest;
                let gamma = self.vertices[c].east;
                self.tjump_hor(j, RectangulationDirection::Right, gamma);
            }
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].north;
            let beta = self.vertices[a].east;
            let b = self.edges[beta].head;
            let c = self.edges[alpha].head;

            if self.vertices[c].type_ == VertexType::Left {
                let gamma = self.vertices[c].north;
                self.wjump_ver(j, RectangulationDirection::Right, gamma);
            } else if self.vertices[b].type_ == VertexType::Top {
                let gamma = self.vertices[b].north;
                self.tjump_ver(j, RectangulationDirection::Right, gamma);
            } else {
                let gamma = self.vertices[c].east;
                self.sjump(j, RectangulationDirection::Right, gamma);
            }
        } else if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].south;
            let b = self.edges[alpha].tail;

            if self.vertices[b].type_ == VertexType::Left {
                let gamma = self.vertices[b].south;
                self.wjump_ver(j, RectangulationDirection::Left, gamma);
            } else {
                let k = self.edges[alpha].left;
                let c = self.rectangles[k].nwest;
                let gamma = self.vertices[c].south;
                self.tjump_ver(j, RectangulationDirection::Left, gamma);
            }
        }
    }

    fn next_diagonal(&mut self, j: usize, dir: RectangulationDirection) {
        let a = self.rectangles[j].nwest;

        if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].south;
            let b = self.edges[alpha].tail;

            if self.vertices[b].type_ == VertexType::Left {
                let gamma = self.vertices[b].west;
                self.tjump_hor(j, RectangulationDirection::Left, gamma);
            } else {
                let c = self.rectangles[j - 1].swest;
                let gamma = self.vertices[c].north;
                self.sjump(j, RectangulationDirection::Left, gamma);
            }
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].east;
            let k = self.edges[alpha].left;
            let b = self.rectangles[k].neast;
            let gamma = self.vertices[b].west;
            self.tjump_hor(j, RectangulationDirection::Right, gamma);
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].east;
            let b = self.edges[alpha].head;

            if self.vertices[b].type_ == VertexType::Top {
                let gamma = self.vertices[b].north;
                self.tjump_ver(j, RectangulationDirection::Right, gamma);
            } else {
                let c = self.rectangles[j - 1].neast;
                let gamma = self.vertices[c].west;
                self.sjump(j, RectangulationDirection::Right, gamma);
            }
        } else if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].south;
            let k = self.edges[alpha].left;
            let b = self.rectangles[k].swest;
            let gamma = self.vertices[b].north;
            self.tjump_ver(j, RectangulationDirection::Left, gamma);
        }
    }

    fn next_baligned(&mut self, j: usize, dir: RectangulationDirection) {
        let a = self.rectangles[j].nwest;
        self.unlock(j, dir);

        if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].south;
            let b = self.edges[alpha].tail;

            if self.vertices[b].type_ == VertexType::Left {
                let gamma = self.vertices[b].west;
                self.tjump_hor(j, RectangulationDirection::Left, gamma);
                let a = self.rectangles[j].nwest;
                let alpha = self.vertices[a].south;
                let b = self.edges[alpha].tail;
                let c = self.rectangles[j - 1].swest;
                let gamma = self.vertices[c].north;
                let c_prime = self.rectangles[j].seast;

                if self.vertices[b].type_ == VertexType::Top
                    && (self.vertices[c_prime].type_ == VertexType::Left
                        || (j == self.n && self.edges[gamma].left == 0))
                {
                    self.sjump(j, RectangulationDirection::Left, gamma);
                }

                self.lock(j, EdgeDir::Hor);
            } else {
                let c = self.rectangles[j - 1].swest;
                let gamma = self.vertices[c].north;
                self.sjump(j, RectangulationDirection::Left, gamma);
                let gamma = self.vertices[c].north;
                let k = self.edges[gamma].left;
                let c_prime = self.rectangles[k].swest;
                let gamma_prime = self.vertices[c_prime].north;
                self.tjump_ver(j, RectangulationDirection::Left, gamma_prime);
                let c = self.rectangles[j - 1].swest;
                let gamma = self.vertices[c].north;
                let a = self.edges[gamma].head;

                if self.vertices[a].type_ == VertexType::Bottom {
                    let c_prime = self.rectangles[j - 2].swest;
                    let gamma_prime = self.vertices[c_prime].north;
                    self.sjump(j - 1, RectangulationDirection::Left, gamma_prime);
                }

                self.lock(j - 1, EdgeDir::Hor);
            }
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Bottom {
            let alpha = self.vertices[a].east;
            let k = self.edges[alpha].left;
            let b = self.rectangles[k].neast;
            let gamma = self.vertices[b].west;
            self.tjump_hor(j, RectangulationDirection::Right, gamma);
            let a = self.rectangles[j].nwest;
            let alpha = self.vertices[a].south;
            let b = self.edges[alpha].tail;
            let beta = self.vertices[b].south;
            let gamma = self.vertices[b].west;
            let c = self.edges[beta].tail;
            let c_prime = self.edges[gamma].tail;

            if self.vertices[c].type_ == VertexType::Top && self.vertices[c_prime].type_ == VertexType::Right {
                let gamma_prime = self.vertices[a].west;
                self.sjump(j - 1, RectangulationDirection::Right, gamma_prime);
            }

            self.lock(j, EdgeDir::Ver);
        } else if dir == RectangulationDirection::Right && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].east;
            let b = self.edges[alpha].head;

            if self.vertices[b].type_ == VertexType::Top {
                let gamma = self.vertices[b].north;
                self.tjump_ver(j, RectangulationDirection::Right, gamma);
                let a = self.rectangles[j].nwest;
                let alpha = self.vertices[a].east;
                let b = self.edges[alpha].head;
                let c = self.rectangles[j - 1].neast;
                let gamma = self.vertices[c].west;
                let c_prime = self.rectangles[j].seast;
                let e = self.rectangles[j - 1].nwest;

                if self.vertices[b].type_ == VertexType::Left
                    && (self.vertices[c_prime].type_ == VertexType::Top
                        || (j == self.n
                            && !(self.vertices[e].type_ == VertexType::Right
                                && self.edges[gamma].tail == e)))
                {
                    self.sjump(j, RectangulationDirection::Right, gamma);
                }

                self.lock(j, EdgeDir::Ver);
            } else {
                let c = self.rectangles[j - 1].neast;
                let gamma = self.vertices[c].west;
                self.sjump(j, RectangulationDirection::Right, gamma);
                let gamma = self.vertices[c].west;
                let k = self.edges[gamma].left;
                let c_prime = self.rectangles[k].neast;
                let gamma_prime = self.vertices[c_prime].west;
                self.tjump_hor(j, RectangulationDirection::Right, gamma_prime);
                let c = self.rectangles[j - 1].neast;
                let gamma = self.vertices[c].west;
                let a = self.edges[gamma].tail;

                if self.vertices[a].type_ == VertexType::Right {
                    let c_prime = self.rectangles[j - 2].neast;
                    let gamma_prime = self.vertices[c_prime].west;
                    self.sjump(j - 1, RectangulationDirection::Right, gamma_prime);
                }

                self.lock(j - 1, EdgeDir::Ver);
            }
        } else if dir == RectangulationDirection::Left && self.vertices[a].type_ == VertexType::Right {
            let alpha = self.vertices[a].south;
            let k = self.edges[alpha].left;
            let b = self.rectangles[k].swest;
            let gamma = self.vertices[b].north;
            self.tjump_ver(j, RectangulationDirection::Left, gamma);
            let a = self.rectangles[j].nwest;
            let alpha = self.vertices[a].east;
            let b = self.edges[alpha].head;
            let beta = self.vertices[b].east;
            let gamma = self.vertices[b].north;
            let c = self.edges[beta].head;
            let c_prime = self.edges[gamma].head;

            if self.vertices[c].type_ == VertexType::Left && self.vertices[c_prime].type_ == VertexType::Bottom {
                let gamma_prime = self.vertices[a].north;
                self.sjump(j - 1, RectangulationDirection::Left, gamma_prime);
            }

            self.lock(j, EdgeDir::Hor);
        }
    }

    fn lock(&mut self, j: usize, dir: EdgeDir) {
        match dir {
            EdgeDir::Hor => {
                let a = self.rectangles[j].neast;
                let b = self.rectangles[j].swest;
                let c = self.rectangles[j].seast;
                let alpha = self.vertices[a].west;
                let beta = self.vertices[b].east;

                if self.vertices[b].type_ != VertexType::Right
                    || self.vertices[c].type_ != VertexType::Left
                    || self.edges[beta].head != c
                {
                    return;
                }

                let d = self.rectangles[j + 1].seast;
                if self.vertices[d].type_ == VertexType::Top {
                    self.sjump(j + 1, RectangulationDirection::Right, alpha);
                }
            }
            EdgeDir::Ver => {
                let a = self.rectangles[j].swest;
                let b = self.rectangles[j].neast;
                let c = self.rectangles[j].seast;
                let alpha = self.vertices[a].north;
                let beta = self.vertices[b].south;

                if self.vertices[b].type_ != VertexType::Bottom
                    || self.vertices[c].type_ != VertexType::Top
                    || self.edges[beta].tail != c
                {
                    return;
                }

                let d = self.rectangles[j + 1].seast;
                if self.vertices[d].type_ == VertexType::Left {
                    self.sjump(j + 1, RectangulationDirection::Left, alpha);
                }
            }
            EdgeDir::None => {}
        }
    }

    fn unlock(&mut self, j: usize, dir: RectangulationDirection) {
        match dir {
            RectangulationDirection::Right => {
                let a = self.rectangles[j].neast;
                let b = self.rectangles[j].seast;
                let c = self.rectangles[j].swest;
                let gamma = self.vertices[c].north;

                if self.vertices[a].type_ == VertexType::Bottom && self.vertices[b].type_ == VertexType::Top {
                    self.sjump(j + 1, RectangulationDirection::Left, gamma);
                }
            }
            RectangulationDirection::Left => {
                let a = self.rectangles[j].swest;
                let b = self.rectangles[j].seast;
                let c = self.rectangles[j].neast;
                let gamma = self.vertices[c].west;

                if self.vertices[a].type_ == VertexType::Right && self.vertices[b].type_ == VertexType::Left {
                    self.sjump(j + 1, RectangulationDirection::Right, gamma);
                }
            }
            RectangulationDirection::None => {}
        }
    }

    fn contains_pattern(&self, j: usize) -> bool {
        for p in &self.patterns {
            match p {
                RectangulationPattern::BrickLeftRight => {
                    if self.contains_brick_leftright(j) {
                        return true;
                    }
                }
                RectangulationPattern::BrickRightLeft => {
                    if self.contains_brick_rightleft(j) {
                        return true;
                    }
                }
                RectangulationPattern::BrickBottomTop => {
                    if self.contains_brick_bottomtop(j) {
                        return true;
                    }
                }
                RectangulationPattern::BrickTopBottom => {
                    if self.contains_brick_topbottom(j) {
                        return true;
                    }
                }
                RectangulationPattern::WmillClockwise => {
                    if self.contains_wmill_clockwise(j) {
                        return true;
                    }
                }
                RectangulationPattern::WmillCounterclockwise => {
                    if self.contains_wmill_counterclockwise(j) {
                        return true;
                    }
                }
                RectangulationPattern::HVertical => {
                    if self.contains_h_vertical(j) {
                        return true;
                    }
                }
                RectangulationPattern::HHorizontal => {
                    if self.contains_h_horizontal(j) {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn contains_brick_leftright(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Right {
            return false;
        }

        let alpha = self.vertices[a].south;
        let b = self.edges[alpha].tail;
        self.vertices[b].type_ == VertexType::Left
    }

    fn contains_brick_rightleft(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Right {
            return false;
        }

        let alpha = self.vertices[a].north;
        let b = self.edges[alpha].head;
        self.vertices[b].type_ == VertexType::Left
    }

    fn contains_brick_bottomtop(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Bottom {
            return false;
        }

        let alpha = self.vertices[a].east;
        let b = self.edges[alpha].head;
        self.vertices[b].type_ == VertexType::Top
    }

    fn contains_brick_topbottom(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Bottom {
            return false;
        }

        let alpha = self.vertices[a].west;
        let b = self.edges[alpha].tail;
        self.vertices[b].type_ == VertexType::Top
    }

    fn contains_wmill_clockwise(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Right {
            return false;
        }

        let alpha = self.vertices[a].north;
        let x = self.edges[alpha].wall;
        let b = self.walls[x].last;
        let beta = self.vertices[b].east;
        let y = self.edges[beta].wall;
        let c = self.walls[y].last;
        let gamma = self.vertices[c].south;
        let z = self.edges[gamma].wall;
        let d = self.walls[z].first;
        let delta = self.vertices[d].west;
        self.edges[delta].right == j
    }

    fn contains_wmill_counterclockwise(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Bottom {
            return false;
        }

        let alpha = self.vertices[a].west;
        let x = self.edges[alpha].wall;
        let b = self.walls[x].first;
        let beta = self.vertices[b].south;
        let y = self.edges[beta].wall;
        let c = self.walls[y].first;
        let gamma = self.vertices[c].east;
        let z = self.edges[gamma].wall;
        let d = self.walls[z].last;
        let delta = self.vertices[d].north;
        self.edges[delta].right == j
    }

    fn contains_h_vertical(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Right {
            return false;
        }

        let mut b = self.rectangles[j].swest;

        while self.vertices[b].type_ != VertexType::Bottom
            && self.vertices[b].type_ != VertexType::Corner
        {
            let mut c = b;
            while self.vertices[c].type_ != VertexType::Right
                && self.vertices[c].type_ != VertexType::Corner
            {
                if self.vertices[c].type_ == VertexType::Top {
                    let mut d = c;
                    while d != b
                        && self.vertices[d].type_ != VertexType::Bottom
                        && self.vertices[d].type_ != VertexType::Corner
                    {
                        if self.vertices[d].type_ == VertexType::Left {
                            return true;
                        }
                        let delta = self.vertices[d].north;
                        d = self.edges[delta].head;
                    }
                }
                let gamma = self.vertices[c].west;
                c = self.edges[gamma].tail;
            }
            let beta = self.vertices[b].north;
            b = self.edges[beta].head;
        }

        false
    }

    fn contains_h_horizontal(&self, j: usize) -> bool {
        let a = self.rectangles[j].nwest;
        if self.vertices[a].type_ != VertexType::Bottom {
            return false;
        }

        let mut b = self.rectangles[j].neast;

        while self.vertices[b].type_ != VertexType::Right
            && self.vertices[b].type_ != VertexType::Corner
        {
            let mut c = b;
            while self.vertices[c].type_ != VertexType::Bottom
                && self.vertices[c].type_ != VertexType::Corner
            {
                if self.vertices[c].type_ == VertexType::Left {
                    let mut d = c;
                    while d != b
                        && self.vertices[d].type_ != VertexType::Right
                        && self.vertices[d].type_ != VertexType::Corner
                    {
                        if self.vertices[d].type_ == VertexType::Top {
                            return true;
                        }
                        let delta = self.vertices[d].west;
                        d = self.edges[delta].tail;
                    }
                }
                let gamma = self.vertices[c].north;
                c = self.edges[gamma].head;
            }
            let beta = self.vertices[b].west;
            b = self.edges[beta].tail;
        }

        false
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of rectangles
    #[arg(short = 'n', long)]
    n: usize,

    /// Base type of rectangulations
    #[arg(short = 't', long, value_enum, default_value_t = RectangulationType::Generic)]
    type_: RectangulationType,

    /// Forbidden patterns
    #[arg(short = 'p', long, value_enum, num_args = 0..)]
    patterns: Vec<RectangulationPattern>,

    /// Number of rectangulations to list (-1 for all)
    #[arg(short = 'l', long, default_value_t = -1)]
    steps: i32,

    /// Quiet output
    #[arg(short = 'q', long, default_value_t = false)]
    quiet: bool,

    /// Output number of rectangles
    #[arg(short = 'c', long, default_value_t = false)]
    count: bool,
}

fn main() {
    let args = Args::parse();

    if args.n < 1 {
        eprintln!("option -n must be followed by an integer from {{1,2,...}}");
        std::process::exit(1);
    }

    if args.type_ == RectangulationType::Baligned
        && args
            .patterns
            .iter()
            .any(|p| matches!(p, RectangulationPattern::BrickLeftRight
                | RectangulationPattern::BrickRightLeft
                | RectangulationPattern::BrickBottomTop
                | RectangulationPattern::BrickTopBottom
                | RectangulationPattern::HVertical
                | RectangulationPattern::HHorizontal))
    {
        eprintln!("patterns -p3 to -p8 unavailable for -t3");
        std::process::exit(1);
    }

    let mut num_rectangulations = 0;
    let mut rect = Rectangulation::new(args.n, args.type_, args.patterns);

    if args.steps == 0 {
        println!("output limit reached");
        std::process::exit(0);
    }

    loop {
        num_rectangulations += 1;
        if !args.quiet {
            rect.print_coordinates();
        } else if num_rectangulations % 10_000_000 == 0 {
            print!(".");
            std::io::stdout().flush().unwrap();
        }

        if !rect.next() {
            break;
        }

        if args.steps >= 0 && num_rectangulations >= args.steps as usize {
            println!("output limit reached");
            break;
        }
    }

    if args.count {
        if args.quiet && num_rectangulations >= 10_000_000 {
            println!();
        }
        println!("number of rectangulations: {}", num_rectangulations);
    }
}
