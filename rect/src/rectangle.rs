#[derive(Debug, Clone, Default)]
pub struct Rectangle {
    pub nwest: i32,
    pub neast: i32,
    pub swest: i32,
    pub seast: i32,
}

impl Rectangle {
    pub fn init(&mut self, neast: i32, seast: i32, swest: i32, nwest: i32) {
        self.nwest = nwest;
        self.neast = neast;
        self.swest = swest;
        self.seast = seast;
    }
}

#[cfg(test)]
mod tests {
    use super::*; // This lets us access the Rectangle struct in this scope

    #[test]
    fn test_default() {
        let rect = Rectangle::default();
        assert_eq!(rect.nwest, 0);
        assert_eq!(rect.neast, 0);
        assert_eq!(rect.swest, 0);
        assert_eq!(rect.seast, 0);
    }

    #[test]
    fn test_init() {
        let mut rect = Rectangle::default();
        rect.init(1, 2, 3, 4);
        assert_eq!(rect.nwest, 4);
        assert_eq!(rect.neast, 1);
        assert_eq!(rect.swest, 3);
        assert_eq!(rect.seast, 2);
    }
}
