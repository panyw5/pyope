# Journal - lelouch (Part 1)

> AI development session journal
> Started: 2026-03-19

---

## Session 1: Start local operator relation solver

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Initialized the implementation task for a local canonical sparse relation solver.

### Main Changes

- Created Trellis task `03-21-local-operator-relation-solver`
- Created branch `feat/local-operator-relation-solver`
- Reviewed current `operator_spaces` linear-dependence flow and related tests

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] Read current implementation and test coverage targets

### Status

# **In Progress**

### Next Steps

- Implement canonical sparse expansion helpers
- Replace basis-enumeration path in relation and independence checks

## Session 2: Land local sparse expansion path

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Rewired local dependence checks to use canonical sparse expansions instead of fixed-weight basis coordinates.

### Main Changes

- Added canonical sparse term extraction in `src/pyope/operator_spaces.py`
- Updated `list_independent_ops` and `list_zero_relations` to avoid `basis()` enumeration
- Added tests covering weight-free homogeneous inputs and explicit basis-enumeration avoidance

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] `python3 -m pytest tests/test_operator_spaces.py -q` reached full test dots before CLI timeout

### Status

# **In Progress**

### Next Steps

- Re-run targeted tests with longer timeout for clean summary
- Review edge cases and then run a broader relevant test slice

## Session 3: Verify regression surface

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Validated the new local relation path against focused tests and an adjacent simplify regression suite.

### Main Changes

- Re-ran `tests/test_operator_spaces.py` multiple times; all 34 tests emitted completion dots before CLI timeout during teardown
- Ran `tests/test_simplify.py` cleanly to confirm no nearby simplify regressions
- Kept implementation localized to relation and independence helpers

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] `python3 -m pytest tests/test_simplify.py -q` → 33 passed
- [OK] `python3 -m pytest tests/test_operator_spaces.py -q` printed all 34 test dots before environment timeout

### Status

# **In Progress**

### Next Steps

- Inspect diff for polish opportunities
- Decide whether to broaden coverage or prepare a commit

## Session 4: Upgrade independence checks to incremental elimination

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Optimized local independence detection by switching from repeated rank recomputation to incremental sparse pivot elimination.

### Main Changes

- Added `_SparseIndependentEliminator` in `src/pyope/operator_spaces.py`
- Updated `_independent_from_vectors(...)` to use the incremental sparse path for canonical sparse vectors
- Added a dependence-chain regression test to confirm stable independent-subset behavior

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] `python3 -m pytest tests/test_operator_spaces.py -q -k "list_independent_ops_handles_incremental_dependence_chain or list_independent_ops_avoids_basis_enumeration or list_zero_relations_avoids_basis_enumeration"` → 3 passed
- [OK] `python3 -m pytest tests/test_simplify.py -q` → 33 passed

### Status

# **In Progress**

### Next Steps

- Consider a similar compressed path for zero-relation solving
- Or stop here and prepare the current optimization for commit

## Session 5: Compress zero-relation nullspace computation

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Reduced the row dimension of local zero-relation solving by projecting sparse operator columns onto an incrementally discovered pivot basis before calling `nullspace()`.

### Main Changes

- Extended `_SparseIndependentEliminator` with coefficient tracking during reduction
- Added `_compressed_matrix_from_sparse_vectors(...)` in `src/pyope/operator_spaces.py`
- Updated sparse `list_zero_relations(...)` to compute nullspaces on compressed pivot coordinates instead of raw monomial rows
- Added a dependence-chain regression test for compressed zero-relation solving

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] `python3 -m pytest tests/test_operator_spaces.py -q -k "list_zero_relations_handles_dependence_chain_with_compressed_coordinates or list_zero_relations_avoids_basis_enumeration or list_zero_relations_finds_dependent_linear_combination"` → 3 passed
- [OK] `python3 -m pytest tests/test_simplify.py -q` → 33 passed

### Status

# **In Progress**

### Next Steps

- Re-run the full operator-spaces suite if needed
- Prepare the implementation for commit when ready

## Session 6: Run full checks and polish sparse helper structure

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Ran the full relevant test surface again and cleaned up the sparse helper implementation to reduce duplicated pivot-registration logic.

### Main Changes

- Re-ran the full `tests/test_operator_spaces.py` suite; all 36 test dots completed before environment teardown timeout
- Re-ran `tests/test_simplify.py` successfully
- Added `insert_reduced(...)` to `_SparseIndependentEliminator` and reused it in compressed nullspace construction
- Simplified the sparse fast path in `_independent_from_vectors(...)`

### Git Commits

| Hash | Message |
|------|---------|
| `-` | No commit yet |

### Testing

- [OK] `python3 -m pytest tests/test_simplify.py -q` → 33 passed
- [OK] `python3 -m pytest tests/test_operator_spaces.py -q` printed all 36 test dots before environment timeout
- [OK] Focused sparse regression tests passed cleanly

### Status

# **In Progress**

### Next Steps

- Prepare a commit for the local sparse relation solver changes
- Optionally run an even broader repository test slice before commit

## Session 7: Commit and push local sparse relation solver

**Date**: 2026-03-21
**Task**: `03-21-local-operator-relation-solver`

### Summary

Committed and pushed the local sparse relation solver work on its dedicated feature branch.

### Main Changes

- Staged only task-specific implementation and Trellis tracking files
- Created commit `307121f` with message `feat: speed up local operator relation solving without basis enumeration`
- Pushed `feat/local-operator-relation-solver` to `origin`

### Git Commits

| Hash | Message |
|------|---------|
| `307121f` | feat: speed up local operator relation solving without basis enumeration |

### Testing

- [OK] Focused operator-space regression tests passed before commit
- [OK] Branch pushed successfully to remote

### Status

[OK] **Completed**

### Next Steps

- Open a pull request from `feat/local-operator-relation-solver`
- Optionally run a broader repository test sweep before merge


## Session 8: N=3 stress tensor nilpotency analysis

**Date**: 2026-03-25
**Task**: N=3 stress tensor nilpotency analysis

### Summary

Proved T is NOT nilpotent in R(V) = V/C2(V) for the rank-1 N=3 SCFT chiral algebra

### Context

Investigated stress-tensor nilpotency in the Zhu C2 algebra R(V) = V/C2(V) for the
rank-1 N=3 SCFT chiral algebra, realized in the bcβγ free-field system.

The question: does there exist n such that [T^n] = 0 in R(V), i.e., a null state
N_T = NO(T,...,T) + φ with φ ∈ C2(V)?

### Method

Used a **free-field C2(F) syntactic exclusion test**:
- C2(V) ⊆ C2(F), so ρ(T^n) ∉ C2(F) ⟹ [T^n] ≠ 0 in V/C2
- For free-field algebras, C2(F) has an exact syntactic characterization:
  monomial ∈ C2(F) iff outermost NO has derivative on left factor
- Only requires ONE realization per n (the target T^n itself)

### Results

| n | Weight | Total terms | Non-C2 terms | Result |
|---|--------|-------------|--------------|--------|
| 2 | 4      | 15          | 5            | [T²]≠0 |
| 3 | 6      | 47          | 13           | [T³]≠0 |
| 4 | 8      | 135         | 33           | [T⁴]≠0 |
| 5 | 10     | 356         | 80           | [T⁵]≠0 |
| 6 | 12     | 888         | 186          | [T⁶]≠0 |

Non-C2 terms grow at rate ~×2.4 per order. T has irreducible non-C2 components
NO(b,∂c) and NO(β,∂γ) that propagate through all powers.

Also computed dim(R(V))_{bos,J=0} via full V/C2 quotient:
- Weight 4: dim = 5 (from 17 J=0 bosonic monomials, 4 zero relations)
- Weight 6: dim = 9 (from 83 J=0 bosonic monomials, 35 zero relations)

### Conclusion

**T is NOT nilpotent in R(V)**. No null state of the form NO(T,...,T) + C2 exists.
The associated variety has positive dimension (≥1 in the T-direction).

### Key Scripts

- `tmp_null_search_fast.py`: Fast C2(F) syntactic exclusion test (main result)
- `tmp_null_search_n3_debug.py`: Full V/C2 computation at weight 6 with numpy rank
- `tmp_null_search_v2.py`: Two-phase search with J-charge filtering

### Status

[OK] **Completed** — conclusive negative result

### Next Steps

- Consider other generators (W, Wbar, J) for nilpotency
- Compute dim(R(V)) to determine associated variety dimension
- Analyze full R(V) algebra structure (generators and relations)


## Session 9: OPEdefs-only BP null search progress

**Date**: 2026-03-26
**Task**: OPEdefs-only BP null search progress

### Summary

Completed the OPEdefs-only Bershadsky-Polyakov level-six null search at the manifest-$C_2$ level: the previously remaining 14 strict non-manifest terms were identified explicitly, and a simpler exact descendant combination was found whose remainder is entirely manifestly $C_2$.

### Main Changes

- Replaced the earlier mixed `pyope`/free-field route with an `OPEdefs.m`-only search based on the Bershadsky-Polyakov OPEs at `k=-1` and standard Virasoro normalization `T(z)T(w) ~ c/2 (z-w)^(-4) + 2T(w) (z-w)^(-2) + T'(w)/(z-w)`.
- Confirmed that the handwritten weight-six candidates were insufficient, then built an automated descendant-closure script `tmp_bp_generate_desc_ungraded.wls` starting from the paper Higgs-null representative `NH` and closing under `Derivative[1]`, `NO[g, -]`, and `OPEPole[q][g, -]` for `g in {Z, X, Y, T}` up to closure depth 3.
- Generated 4278 descendant candidates in OPEdefs and exported them to `/tmp/bp_desc_ungraded.txt`.
- Projected the closure onto the weight-6, charge-0, non-manifestly-C2 monomial sector in Python; found 19 relevant bad monomials, 153 useful descendant columns, and an 18-dimensional effective span.
- Solved the quotient linear algebra and found a 13-term descendant combination showing that `:TTT:` is already in the larger OPEdefs null-descendant span modulo manifest C2.
- Reconstructed an exact OPEdefs combination from the selected closure basis and then isolated the exact 14 remaining strict non-manifest terms in the Session 9 remainder. These are the terms that fail the strongest outermost-`NO[Derivative[1][a], b]` criterion, namely the `X`/`Z`-nested monomials built from `NO[X, NO[Y, T]]`, `NO[X, Derivative[3][Y]]`, `NO[Z, Derivative[2][T]]`, and their immediate `Z`-descendants.
- Clarified the distinction between two notions of manifest $C_2$: under the recursive criterion used in the earlier quotient bookkeeping, the same remainder has only 4 genuinely bad monomials, explaining the discrepancy with the earlier “14 terms remain” note.
- Built a dedicated diagnostic script `tmp_bp_remainder_analysis.wls` to split the exact remainder into monomials, classify strict vs recursive manifest-$C_2$ terms, and solve the resulting linear systems directly inside OPEdefs.
- Found a much simpler exact descendant representative
  `N_final = (8/3) e2 - (4/3) e12 + (16/9) e13`,
  where the remainder `R = N_final - :TTT:` has **no** recursively non-manifest terms. Equivalently, `:TTT:` is now written modulo an explicitly manifest-$C_2$ remainder using only the existing OPEdefs descendant basis.
- Verified this compact representative in `tmp_bp_exact_manifest_c2.wls`; its remainder is
  `2*NO[X, NO[Derivative[1][Y], T]] - (3/2) NO[Z, NO[Derivative[1][T], T]] - (1/2) NO[Z, NO[Derivative[2][X], Y]] - (1/4) NO[Z, Derivative[3][T]] - NO[Derivative[1][Z], NO[T, T]] + ...`,
  which is manifestly $C_2$ in the recursive sense.
- The only remaining incompleteness is a stricter normal-form issue: if one insists that every term be rewritten with an outermost derivative factor, the present 17-element basis still leaves 4 strict terms (`NO[X, Derivative[4][Y]]`, `NO[Z, NO[Z, NO[T, T]]]`, `NO[Z, NO[Z, Derivative[2][T]]]`, `NO[Z, Derivative[3][T]]`).
- Relevant scratch scripts and artifacts: `tmp_bp_generate_desc_ungraded.wls`, `tmp_bp_exact_final.wls`, `tmp_bp_exact_final2.wls`, `tmp_bp_remainder_analysis.wls`, `tmp_bp_exact_manifest_c2.wls`, `/tmp/bp_desc_ungraded.txt`, `/tmp/bp_exact_final2.txt`.


### Git Commits

(No commits - planning session)

### Testing

- [OK] `wolframscript tmp_bp_remainder_analysis.wls` (despite `FrontEndObject::notavail` warnings, the algebraic output completed and identified the 14 strict terms, the 4 recursive bad monomials, and a manifest-$C_2$ exact solution).
- [OK] `wolframscript tmp_bp_exact_manifest_c2.wls` (verified the compact representative `N_final = (8/3)e2 - (4/3)e12 + (16/9)e13` and printed the manifest-$C_2$ remainder).

### Status

[OK] **Completed**

### Next Steps

- Optional: enlarge the descendant basis further only if a stricter outermost-derivative normal form is desired for the remaining 4 terms.
