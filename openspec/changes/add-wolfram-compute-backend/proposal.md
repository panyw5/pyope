# Proposal: Add wolframscript compute backend for intensive OPE operations

## Why

`pyope` currently evaluates all `OPE(...)` and `NO(...)` computations in Python/SymPy. This keeps the implementation self-contained, but composite-operator calculations become expensive quickly, especially for nested normal products and Jacobi-driven recursive expansions. The repository already treats `OPEdefs.m` as the mathematical reference implementation, and `wolframscript` is available as a practical execution path on developer machines.

For performance-oriented work, we need a second execution backend that can delegate computation-intensive operations to Wolfram Language while preserving the existing Python API and mathematical results.

## What Changes

- Add a configurable compute backend layer for core symbolic operations
- Keep the existing Python/SymPy implementation as the default backend
- Add a `wolframscript` backend for MVP support of binary `OPE(A, B)` and `NO(A, B)`
- Add a Wolfram wrapper script that loads `OPEdefs.m` and evaluates requested operations in a CLI-friendly way
- Add Python-side expression encoding/decoding between `pyope` expressions and Wolfram Language input/output
- Add backend availability checks, failure reporting, and a clear fallback strategy
- Add tests that compare the Wolfram backend path against the existing SymPy backend on representative examples

## Scope

### In scope for MVP

- Backend selection/configuration in Python
- Binary `OPE(A, B)` dispatch through the Wolfram backend
- Binary `NO(A, B)` dispatch through the Wolfram backend
- Synchronizing basic operator declarations and registered OPE definitions needed by the delegated computation
- A repository-local `OPEdefs.wls` wrapper for `wolframscript`
- Tests that validate backend equivalence on representative operator expressions

### Out of scope for MVP

- Replacing the full simplification pipeline with Wolfram Language
- Delegating every helper (`bracket`, commute helpers, simplify helpers) independently
- Persistent Wolfram kernel sessions or daemon-style IPC
- Automatic backend selection based on benchmarking heuristics
- Shipping Wolfram Engine as a project dependency

## Impact

- Affected specs: new compute-backend capability; modified normal-product behavior documentation where backend semantics matter
- Affected code: `src/pyope/api.py`, new backend/config modules, expression codec utilities, Wolfram wrapper assets, tests
- Affected workflow: developers can opt into the Wolfram backend when `wolframscript` is installed

## Risks

- Expression serialization mismatches between SymPy/pyope and Wolfram Language
- Runtime failures when `wolframscript` is missing or `OPEdefs.m` cannot be loaded
- Divergence between Python registry state and Wolfram-side declared state
- Increased latency for trivial computations if delegation is overused

## Mitigations

- Limit MVP delegation to binary `OPE` and `NO`
- Keep `sympy` as the default backend and preserve explicit fallback paths
- Add backend equivalence tests against representative composite examples
- Make backend activation explicit and observable via configuration and error messages
