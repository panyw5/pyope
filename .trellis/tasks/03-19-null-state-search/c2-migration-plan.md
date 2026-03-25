# C2 Migration Plan

## Goal

Migrate the older basis-driven `C2Space` / `C2NullSearcher` APIs toward the new sparse reducer-driven `AbstractC2Reducer` / `GenericC2Reducer` / precheck `C2NullSearcher` stack without breaking existing callers.

## Migration Strategy

### Phase A: compatibility first

- Keep legacy `src/pyope/operator_spaces.py` classes public and working.
- Add reducer-backed compatibility methods onto legacy `C2Space` so callers can use quotient-style APIs without immediately switching modules.
- Add a quotient-precheck bridge onto legacy `C2NullSearcher` so old callers can gradually adopt the new precheck workflow.

### Phase B: old API delegates to new core where safe

- Make legacy `C2Space.contains(...)` capable of delegating to `GenericC2Reducer.is_zero_mod_c2(...)` for membership checks.
- Expose `quotient_normal_form(...)`, `is_zero_mod_c2(...)`, and `solve_c2_witness(...)` on legacy `C2Space` as thin compatibility helpers.
- Store an optional reducer on legacy `C2NullSearcher` and default it from the basis/canonicalizer when absent.

### Phase C: result-shape convergence

- Add `quotient_precheck(...)` to legacy `C2NullSearcher` returning the new structured `NullSearchResult`.
- Keep old `search_from_sources(...)` result dict unchanged for now.
- Keep old `search_stress_tensor_nilpotency(...)` behavior unchanged for now, while allowing the new precheck-only searcher to coexist.

### Phase D: later full migration

- Move descendant solving and singular constraints onto sparse/reducer abstractions.
- Retire dense `D x - C y = t` assembly as the primary implementation path.
- Eventually reduce legacy classes to wrappers or aliases once the new searcher supports explicit lift/witness reconstruction.

### Phase E: descendant lift on new core

- Add a new-core `search_from_sources(...)` that first computes the quotient target via the reducer.
- Solve for descendant combinations in the quotient using sparse canonical supports, rather than explicit fixed-weight `C2` bases.
- Reconstruct the final `C2` witness from `null_operator - target` after the descendant lift succeeds.
- Let legacy `operator_spaces.C2NullSearcher.search_from_sources(...)` delegate to this new-core lift path and then adapt the result back to the legacy dict shape.

## Immediate Coding Scope

- Add a reducer factory/helper to legacy `C2Space`.
- Add quotient compatibility methods to legacy `C2Space`.
- Add `quotient_precheck(...)` and reducer plumbing to legacy `C2NullSearcher`.
- Add tests covering old-to-new bridge behavior.
- Add minimal descendant-lift solving on the new searcher and wire legacy search to it.

## Non-Goals For This Session

- Do not remove legacy APIs.
- Do not rewrite descendant solving yet.
- Do not attempt realization-aware acceleration yet.

## Bigger Picture Remaining Work

- Align the migration with `null-c2-search-plan.md` Workstreams 3-5, not just legacy compatibility.
- Extract descendant generation into a dedicated sparse-first module once the new-core API shape stabilizes.
- Extract singular/null OPE constraints into a reusable constraint-system module and let the high-level searcher consume it.
- Converge legacy dict results and new structured result types into one public story.
- Add README/public examples once the stress-tensor and explicit-null APIs stabilize.
- Add realization-aware acceleration only after the generic sparse path and singular-constraint pipeline are stable.

## Progress Notes

- Legacy `C2Space` now exposes reducer-backed `quotient_normal_form(...)`, `is_zero_mod_c2(...)`, and `solve_c2_witness(...)` compatibility helpers.
- Legacy `C2NullSearcher` now exposes `quotient_precheck(...)` and routes descendant lifting through the new quotient-first core before adapting results back to the legacy payload shape.
- `descendants.py` extracted; `operator_spaces.DescendantSpace` now serves as compatibility import surface.
- `singularity.py` extracted; `operator_spaces.SingularVectorAnalyzer` now serves as compatibility import surface.
- `NullSearchResult` now owns legacy-payload adaptation so old searcher dict results share one bridge.
- Next cleanup target is realization/backend layering rather than more `operator_spaces.py` duplication removal.
- `realizations.py` now provides a minimal `RealizationBackend` abstraction plus an `IdentityRealizationBackend` placeholder.
- `free_field_c2.py` now provides a minimal `FreeFieldC2Reducer` shell that safely falls back to `GenericC2Reducer` until true free-field quotient logic is implemented.
- First realized quotient rule is now implemented: derivative factors are killed in the quotient via `DerivativeKillingRealizationBackend` / `DerivativeKillingFreeFieldC2Reducer`.
- Next Workstream 4 step should compare generic and free-field reducers on overlapping tractable examples and then extend beyond derivative-killing rules.
