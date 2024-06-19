/// Represents a rectangle with four boundaries defined by their coordinates.
#[derive(Debug, Clone, Default)]
pub struct Rectangle {
    /// The northwest coordinate.
    pub nwest: i32,
    /// The northeast coordinate.
    pub neast: i32,
    /// The southwest coordinate.
    pub swest: i32,
    /// The southeast coordinate.
    pub seast: i32,
}

impl Rectangle {
    /// Constructs a new `Rectangle` with default coordinate values (0).
    pub fn new() -> Self {
        Rectangle {
            nwest: 0,
            neast: 0,
            swest: 0,
            seast: 0,
        }
    }

    /// Initializes the rectangle with specified coordinates.
    ///
    /// # Arguments
    ///
    /// * `neast` - The northeast coordinate.
    /// * `seast` - The southeast coordinate.
    /// * `swest` - The southwest coordinate.
    /// * `nwest` - The northwest coordinate.
    pub fn init(&mut self, neast: i32, seast: i32, swest: i32, nwest: i32) {
        self.nwest = nwest;
        self.neast = neast;
        self.swest = swest;
        self.seast = seast;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rectangle_init() {
        let mut rect = Rectangle::new();
        rect.init(10, 20, 30, 40);

        assert_eq!(
            rect.nwest, 40,
            "The northwest coordinate was not set correctly."
        );
        assert_eq!(
            rect.neast, 10,
            "The northeast coordinate was not set correctly."
        );
        assert_eq!(
            rect.swest, 30,
            "The southwest coordinate was not set correctly."
        );
        assert_eq!(
            rect.seast, 20,
            "The southeast coordinate was not set correctly."
        );
    }
}
