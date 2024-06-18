#[derive(Debug, Clone, Default)]
pub struct Wall {
    pub first: i32,
    pub last: i32,
}


impl Wall {
    pub fn init(&mut self, first: i32, last: i32) {
        self.first = first;
        self.last = last;
    }
}

#[cfg(test)]
mod tests {
    use super::*; // bring everything from parent scope into this one

    #[test]
    fn test_wall_default() {
        let wall = Wall::default();

        assert_eq!(wall.first, 0);
        assert_eq!(wall.last, 0);
    }

    #[test]
    fn test_wall_init() {
        let mut wall = Wall::default();
        wall.init(3, 5);

        assert_eq!(wall.first, 3);
        assert_eq!(wall.last, 5);
    }
}
