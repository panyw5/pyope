## 1. Specification
- [ ] 1.1 Add compute-backend spec describing backend selection and supported operations
- [ ] 1.2 Add/modify normal-product spec notes where backend-backed `NO(A, B)` behavior must remain equivalent
- [ ] 1.3 Validate the OpenSpec change with `openspec validate add-wolfram-compute-backend --strict --no-interactive`

## 2. Backend infrastructure
- [ ] 2.1 Add backend configuration module with default `sympy`
- [ ] 2.2 Add backend dispatcher for binary `OPE` and `NO`
- [ ] 2.3 Refactor current Python compute path behind backend-local functions without changing behavior

## 3. Wolfram bridge
- [ ] 3.1 Add `OPEdefs/OPEdefs.wls` wrapper script for `wolframscript`
- [ ] 3.2 Add Python expression encoder from `pyope` objects to Wolfram input
- [ ] 3.3 Add Python expression decoder from Wolfram output back to `pyope`/SymPy
- [ ] 3.4 Synchronize required operator declarations and registered OPE definitions for delegated computations
- [ ] 3.5 Add backend availability checks and clear error handling

## 4. Integration
- [ ] 4.1 Route `OPE(A, B)` through the selected backend
- [ ] 4.2 Route binary `NO(A, B)` through the selected backend
- [ ] 4.3 Preserve existing multi-argument `NO(...)` behavior through binary composition

## 5. Validation
- [ ] 5.1 Add unit tests for backend selection/configuration
- [ ] 5.2 Add equivalence tests between `sympy` and `wolfram` backends for representative `OPE` examples
- [ ] 5.3 Add equivalence tests between `sympy` and `wolfram` backends for representative `NO` examples
- [ ] 5.4 Add tests for missing `wolframscript` / backend initialization failure paths
- [ ] 5.5 Run focused test suite for backend-related files
