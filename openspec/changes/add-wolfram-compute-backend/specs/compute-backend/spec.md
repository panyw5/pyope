# compute-backend Specification

## ADDED Requirements

### Requirement: Support configurable compute backends for intensive symbolic operations

The system SHALL provide a configurable compute backend for intensive symbolic operations while preserving the public `pyope` API.

#### Scenario: Default backend remains sympy

**Given** the user imports `pyope`
**When** no backend configuration is provided
**Then** `OPE(A, B)` and `NO(A, B)` SHALL use the built-in SymPy-based implementation

#### Scenario: User selects wolfram backend

**Given** `wolframscript` is available on the current machine
**And** the user selects the `wolfram` backend
**When** the user computes `OPE(A, B)` or `NO(A, B)` for an expression supported by the MVP bridge
**Then** the system SHALL delegate the computation through the Wolfram backend
**And** return the result as a regular `pyope`/SymPy expression

#### Scenario: Wolfram backend is unavailable

**Given** the user selects the `wolfram` backend
**When** `wolframscript` is not installed or the wrapper script cannot initialize `OPEdefs.m`
**Then** the system SHALL raise a clear backend-related error
**And** the error message SHALL explain the missing dependency or initialization failure

### Requirement: Preserve mathematical equivalence across backends for supported operations

For operations supported by the MVP bridge, the Wolfram backend SHALL produce results mathematically equivalent to the SymPy backend.

#### Scenario: Equivalent binary OPE result

**Given** a set of registered basic operators and OPE definitions
**When** the user computes the same binary `OPE(A, B)` through the `sympy` backend and the `wolfram` backend
**Then** both results SHALL represent the same OPE poles and coefficients

#### Scenario: Equivalent binary normal-ordered result

**Given** a pair of operators or operator expressions supported by the MVP bridge
**When** the user computes the same binary `NO(A, B)` through the `sympy` backend and the `wolfram` backend
**Then** both results SHALL represent the same mathematical normal-ordered expression

### Requirement: Preserve existing multi-argument normal product semantics

The addition of a second backend SHALL NOT change the public semantics of `NO(A, B, C, ...)` or `normal_product(...)`.

#### Scenario: Multi-argument NO remains right-associated

**Given** the user calls `NO(A, B, C)`
**When** the system evaluates the expression under any supported backend configuration
**Then** the result SHALL be computed as the right-associated composition of binary normal ordering
**And** preserve the same semantics as before the backend abstraction was added
