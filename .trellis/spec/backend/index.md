# Backend Development Guidelines

> **Note**: This project is a mathematical computation library, not a web backend.
> The original Trellis backend templates (database, logging, etc.) are not applicable.
>
> **Actual guidelines are in `core/` and `testing/` directories.**

---

## Redirects

| Original Template | Replaced By | Location |
|-------------------|------------|----------|
| Directory Structure | Architecture | [`core/architecture.md`](../core/architecture.md) |
| Database Guidelines | *(not applicable)* | — |
| Error Handling | Quality Guidelines (error section) | [`core/quality-guidelines.md`](../core/quality-guidelines.md) |
| Quality Guidelines | Quality Guidelines | [`core/quality-guidelines.md`](../core/quality-guidelines.md) |
| Logging Guidelines | *(not applicable)* | — |

---

## Project-Specific Spec Directories

| Directory | Purpose |
|-----------|---------|
| [`core/`](../core/index.md) | Architecture, API conventions, SymPy patterns, quality |
| [`testing/`](../testing/index.md) | Test conventions, Mathematica reference, fixtures |
| [`domain/`](../domain/index.md) | VOA/CFT domain knowledge for AI agents |
| [`guides/`](../guides/index.md) | General thinking guides (code reuse, cross-layer) |
