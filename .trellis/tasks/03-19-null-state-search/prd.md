# Advance null state search API

## Goal

Design a reusable local-operator API in `pyope` for searching singular descendants and $C_2$ null relations, staying in the language of OPEs, derivatives, and normal-ordered products.

## Requirements

- Define a canonical fixed-weight local-operator basis API with coordinate extraction.
- Define a generic descendant-space API that works for both Virasoro and extended strong-generator VOAs.
- Define a generic singular-vector analyzer based on OPE pole constraints.
- Define a fixed-weight $C_2$ space API with basis and membership checks.
- Define a null-search API that solves $
  \mathcal N \in \mathrm{Desc}_h,
  \quad
  \mathcal N - T^{(n)} \in C_2(V)_h
  $ using linear algebra in basis coordinates.
- Keep the design API-first and avoid mode-algebra-first machinery.

## Acceptance Criteria

- [x] `prd.md` captures the task goal and scope in Trellis format.
- [x] The task has initialized context files for later implement/check agents.
- [x] Relevant specs and existing code patterns are identified.
- [x] The task is activated as the current Trellis task.
- [x] A concrete next implementation plan is ready for the MVP.

## Current Status

- Planning artifacts are complete and the task is active in Trellis.
- Sparse infrastructure for the search stack has landed, including public sparse-term access, `SparseLinearContext`, and tests that avoid ambient basis enumeration.
- The new quotient-first core is implemented in `src/pyope/c2.py` and `src/pyope/null_search.py`, with `AbstractC2Reducer`, `GenericC2Reducer`, `NullSearchResult`, `quotient_precheck(...)`, and structured stress-tensor search metadata.
- Legacy `C2Space` / `C2NullSearcher` compatibility bridges now expose reducer-backed quotient checks and delegate the descendant lift path through the new core where available.
- `src/pyope/descendants.py` and `src/pyope/singularity.py` have been extracted, and the first realization-aware layer now exists in `src/pyope/realizations.py` and `src/pyope/free_field_c2.py`.
- Remaining work focuses on deeper `C_2` subspace progress: removing remaining basis-dependent lift logic, extending realization-aware quotient rules beyond derivative killing, and recovering explicit known null relations instead of only toy lift witnesses.

## Technical Notes

- Use `null-c2-search-plan.md` as the source design document.
- Treat Virasoro as the one-generator special case of the generic API.
- Favor reuse of existing `pyope` primitives such as `NO`, `normal_product`, and `simplify`.
- The natural target module is likely `src/pyope/operator_spaces.py` or a small related module cluster.
