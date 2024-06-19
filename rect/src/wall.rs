/// Represents a wall with a start and end point.
#[derive(Debug, Default)]
pub struct Wall {
    /// The starting index of the wall.
    pub first_: i32,
    /// The ending index of the wall.
    pub last_: i32,
}

impl Wall {
    /// Creates a new instance of `Wall` with default values.
    ///
    /// # Examples
    ///
    /// ```
    /// use rect::wall::Wall;
    ///
    /// let wall = Wall::new();
    /// assert_eq!(wall.first_, 0);
    /// assert_eq!(wall.last_, 0);
    /// ```
    pub fn new() -> Self {
        Wall::default()
    }

    /// Initializes the wall with specified start and end points.
    ///
    /// # Arguments
    ///
    /// * `first` - The starting index of the wall.
    /// * `last` - The ending index of the wall.
    ///
    /// # Examples
    ///
    /// ```
    /// use rect::wall::Wall;
    ///
    /// let mut wall = Wall::new();
    /// wall.init(1, 10);
    /// assert_eq!(wall.first_, 1);
    /// assert_eq!(wall.last_, 10);
    /// ```
    pub fn init(&mut self, first: i32, last: i32) {
        self.first_ = first;
        self.last_ = last;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wall_initialization() {
        let mut wall = Wall::new();
        wall.init(1, 10);

        assert_eq!(wall.first_, 1);
        assert_eq!(wall.last_, 10);
    }
}
