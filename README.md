# gray-code-rs

Rust libraries for Gray code and rectangulation problems.

## Crates

### rect

Rectangulation library - handling rectangles, vertices, edges, and walls.

```rust
use rect::{Rectangle, Vertex, Edge, Wall};

let mut rect = Rectangle::new();
rect.init(10, 20, 30, 40);
```

### middle

Gray code and Hamiltonian cycle library - vertices and trees.

```rust
use middle::{Vertex, Tree};

let v = Vertex::new(vec![1, 0, 1]);
let tree = Tree::new(vec![1, 0, 1, 0, 1]);
```

## Development

```bash
cargo build
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check
```

## License

MIT