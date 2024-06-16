#[derive(Debug, Clone)]
pub struct Rectangle {
    pub nwest: i32,
    pub neast: i32,
    pub swest: i32,
    pub seast: i32,
}

impl Default for Rectangle {
    fn default() -> Self {
        Self {
            nwest: 0,
            neast: 0,
            swest: 0,
            seast: 0,
        }
    }
}

impl Rectangle {
    pub fn init(&mut self, neast: i32, seast: i32, swest: i32, nwest: i32) {
        self.nwest = nwest;
        self.neast = neast;
        self.swest = swest;
        self.seast = seast;
    }
}