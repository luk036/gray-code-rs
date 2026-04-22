# AGENTS.md - Agent Guidelines for gray-code-rs

This file provides coding guidelines for AI agents operating in this Rust workspace.

## Project Overview

- **`rect`** - Library for rectangulation problems (rectangles, vertices, edges, walls)
- **`middle`** - Library for Gray code/Hamiltonian cycle problems (vertices, trees, hamcycle)

Both use Rust edition 2021.

---

## Build, Lint, and Test Commands

```bash
# Run all tests
cargo test --workspace

# Run tests for a specific crate
cargo test -p rect
cargo test -p middle

# Run a single test
cargo test test_rectangle_init
cargo test test_vertex_init_corner
cargo test testedgeinitialization
cargo test test_is_first_vertex

# Build
cargo build
cargo build -p rect
cargo build -p middle
cargo build --release

# Lint
cargo clippy --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings

# Format
cargo fmt --check
cargo fmt

# Documentation
cargo doc --no-deps
cargo doc --no-deps --open
```

---

## Code Style Guidelines

### Documentation

- Use doc comments (`///`) for public types/functions
- Include `# Arguments`, `# Panics`, `# Returns` sections as needed
- Keep line width under 100 characters

### Structs and Enums

- Use PascalCase for types
- Use `#[derive(Debug, Clone, Default)]` for most structs
- Use `#[derive(Debug, PartialEq)]` for enums and structs needing equality

```rust
#[derive(Debug, Clone, Default)]
pub struct Rectangle { pub nwest: i32, pub neast: i32, pub swest: i32, pub seast: i32 }

#[derive(Debug, PartialEq)]
pub enum VertexType { Corner, Top, Bottom, Left, Right, None }
```

### Fields and Variables

- Use snake_case: `nwest`, `neast`, `type_`
- Avoid single-letter names except for indices (i, j, n)
- Prefix unused fields with underscore: `#[allow(dead_code)]`

### Constructor Patterns

- `new()` - returns default/empty instance
- `init()` - multi-step initialization that mutates self

```rust
impl Rectangle {
    pub fn new() -> Self { Rectangle { nwest: 0, neast: 0, swest: 0, seast: 0 } }
    pub fn init(&mut self, neast: i32, seast: i32, swest: i32, nwest: i32) { /* ... */ }
}
```

### Error Handling

- Use `assert!` for invariants that should never fail
- Use `Result<T, E>` for recoverable errors
- Use `unreachable!()` for exhaustive match cases

### Pattern Matching

- Use exhaustive match patterns
- Use `_` for catch-all only when appropriate

### Imports

- Group: std lib → external crates → local modules
- Use `use crate::module::Type;`

---

## Module Organization

- One module per file (`mod foo;` → `foo.rs`)
- Tests in same file: `#[cfg(test)] mod tests { use super::*; ... }`
- Test functions: `fn test_function_name()`

```
rect/src/
├── lib.rs, main.rs
├── rectangle.rs, vertex.rs, edge.rs, wall.rs, rectangulation.rs
```

---

## Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Struct/Enum | PascalCase | `Rectangle`, `VertexType` |
| Enum Variant | PascalCase | `Corner`, `Hor` |
| Function/Method | snake_case | `new()`, `flip_bit()` |
| Field/Variable | snake_case | `nwest`, `type_` |
| Test | snake_case + test_ prefix | `test_rectangle_init` |

Note: Use trailing underscore when field name conflicts with type: `type_`

---

## Clippy and Pre-commit

### Common Allow Attributes

```rust
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
```

### Pre-commit Checklist

1. `cargo fmt` passes
2. `cargo clippy --all-targets -- -D warnings` passes
3. `cargo test --workspace` passes
4. `cargo doc --no-deps` succeeds

---

## Additional Notes

- No custom rustfmt.toml or clippy.toml - using Rust defaults
- No CI/CD beyond GitHub Pages workflow
- No README or LICENSE - add if needed