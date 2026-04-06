## Context

The current implementation computes OPE and normal ordering directly in Python using recursive rules in `src/pyope/api.py`. This is mathematically aligned with `OPEdefs.m`, but the hottest branches are exactly the composite-operator formulas that `OPEdefs.m` already implements efficiently.

We want a second backend without breaking the public API:

- `OPE(A, B)` stays the entry point
- `NO(A, B)` stays the entry point
- existing code using `pyope` should not need to change expression-building patterns

## Goals / Non-Goals

- Goals:
  - Add an explicit compute backend abstraction
  - Support `sympy` and `wolfram` backends
  - Preserve mathematical equivalence for supported operations
  - Keep the Wolfram integration narrow and testable
- Non-Goals:
  - Full bidirectional feature parity for every internal helper in MVP
  - Long-lived kernel management
  - Automatic cost-model based routing

## Decisions

- Decision: introduce a backend dispatcher at the public compute boundary, not inside every helper
  - Why: this minimizes invasive changes and keeps existing Python logic intact as the reference fallback path

- Decision: MVP delegation target is binary `OPE(A, B)` and binary `NO(A, B)` only
  - Why: multi-argument `NO(...)` already reduces to binary primitives; this keeps the new surface area small

- Decision: use a repository-local `OPEdefs.wls` wrapper around `OPEdefs.m`
  - Why: `wolframscript` executes `.wls` scripts naturally, and the wrapper can standardize CLI arguments, initialization, and output formatting

- Decision: use text-based expression serialization for MVP
  - Why: it is easier to inspect and debug than a binary protocol; correctness matters more than transport efficiency at this stage

- Decision: default backend remains `sympy`
  - Why: `pyope` must still work without Wolfram tooling installed

## Proposed Architecture

1. Add backend configuration module
   - current backend getter/setter
   - default `sympy`
   - optional environment-variable override

2. Add backend dispatcher module
   - `compute_ope(left, right)`
   - `compute_no(left, right)`
   - route to `sympy` or `wolfram`

3. Preserve current Python implementation as backend-local functions
   - move or wrap current `_compute_ope(...)`
   - move or wrap current `_NO_binary(...)`

4. Add Wolfram bridge module
   - encode `pyope` expressions to Wolfram syntax
   - materialize operator/OPE definitions needed for the request
   - run `wolframscript OPEdefs.wls ...`
   - parse the returned expression back into `pyope`/SymPy objects

5. Add `OPEdefs/OPEdefs.wls`
   - load `OPEdefs.m`
   - accept operation + payload
   - emit a stable machine-readable textual form

## Alternatives Considered

- Call `wolframscript -code` directly from Python without a wrapper
  - Rejected because quoting, initialization, and debugging become fragile

- Replace internal helpers one by one with subprocess calls
  - Rejected because it would scatter backend logic across hot recursive paths

- Always compute both backends and compare
  - Rejected for MVP because it doubles runtime and is better suited to tests/debug mode

## Risks / Trade-offs

- Text codecs are easier to debug but can be brittle if the expression grammar is underspecified
- Spawning `wolframscript` per request is simpler but may add overhead
- Registry synchronization introduces state-management complexity

## Migration Plan

1. Add spec + design + tasks
2. Introduce backend config and dispatcher with `sympy` default only
3. Add Wolfram wrapper script and backend bridge
4. Enable delegated `OPE`/`NO` for supported cases
5. Add equivalence and failure-path tests

## Open Questions

- Whether backend selection should be process-global only, or also available as a scoped context manager
- Whether unsupported expressions in MVP should hard-fail or transparently fall back to `sympy`
