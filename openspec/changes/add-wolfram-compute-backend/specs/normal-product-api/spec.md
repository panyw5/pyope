## MODIFIED Requirements

### Requirement: Support variable-argument normal product calculation

The system MUST provide a public API function `normal_product()` that computes the nested normal ordered product of multiple operators. The function MUST accept zero or more operator arguments and return their right-associative normal ordered product. Backend selection MUST NOT change the public semantics of the returned expression.

#### Scenario: Compute normal product through backend-backed binary NO

**Given** the pyope library is imported
**And** a backend is configured for binary `NO(A, B)` evaluation
**When** the user calls `normal_product(A, B, C, D)`
**Then** the function SHALL still build the right-associative structure `NO(A, NO(B, NO(C, D)))`
**And** any backend delegation used for the binary steps SHALL preserve the same mathematical result as the default implementation
