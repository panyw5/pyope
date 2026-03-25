# Null / C2 Search API Plan

## Goal

Build a reusable API layer on top of `pyope` for programmatic searches of null vectors of the form

$$
\mathcal N = \mathcal N_{\mathrm{desc}} = T^{(n)} + \varphi,
\qquad
\varphi \in C_2,
$$

where `T^{(n)}` denotes the canonical right-nested local operator

$$
T^{(n)} = \mathrm{NO}(T, \mathrm{NO}(T, \dots, \mathrm{NO}(T,T)\dots)).
$$

The implementation should stay in the language of local operators, OPEs, derivatives, and normal-ordered products, rather than introducing mode-algebra-first machinery.

The design must be layered:

- provide generic algorithms for local-operator linear algebra, $C_2(V)$, and $V/C_2(V)$ that apply to general strong-generator VOAs;
- optionally attach a realization-aware backend for high-performance reduction when the VOA is embedded into a free-field system such as $bc\beta\gamma$.

## Core Design Constraints

### 1. No dependence on complete fixed-weight bases

The API must not assume that one can efficiently enumerate a complete basis of weight-$h$ local operators.

This is essential for free-field realizations with zero-weight or negative-weight generators, where complete fixed-weight basis enumeration may be combinatorially useless or require arbitrary cutoffs.

Accordingly:

- the primary representation of an expression is a sparse canonical coefficient map over the monomials that actually occur in that expression;
- linear algebra is performed on demand over the finite union of supports appearing in the current search problem;
- fixed-weight basis enumeration is at most an optional helper for small positive-weight toy examples, not a core dependency of the null-search pipeline.

### 2. Generic VOA first, optimized realization second

The architecture must separate:

- a generic local-operator/$C_2$/quotient layer that works without any realization data;
- a realization backend layer that can override reduction and quotient-normal-form computations when explicit free-field formulas are available.

Correctness must come from the generic layer. Realization-aware code is an acceleration path, not the semantic definition of the API.

### 3. Sparse local coordinates, not ambient-basis coordinates

The word "coordinates" in this task means:

- sparse canonical expansions such as
  `expr -> {monomial: coefficient}`;
- temporary coordinate matrices assembled from a finite family of queried expressions;
- incremental elimination and relation solving on those temporary supports.

It does **not** mean coordinates with respect to a precomputed complete basis of `V_h`.

## Guiding Principles

- Reuse existing `pyope` functionality for `NO`, `normal_product`, `simplify`, canonical normal-order reordering, and local sparse relation solving.
- Treat Virasoro as the one-generator special case of the generic strong-generator framework.
- Represent descendants, singularity constraints, $C_2$ membership, and quotient classes entirely in local-operator language.
- Expose both correctness-oriented generic reducers and realization-accelerated reducers behind the same public abstractions.

## High-Level Layering

The planned API should be organized into the following layers.

### Layer A. Generic local-operator infrastructure

This layer knows about local operators, normal ordering, derivatives, and OPEs, but not about any particular realization.

Key responsibilities:

- canonicalize expressions;
- extract sparse canonical term maps;
- perform incremental sparse elimination;
- generate descendants and singularity constraints.

### Layer B. Generic $C_2$ / quotient infrastructure

This layer is still realization-free.

Key responsibilities:

- generate $C_2$ candidates on demand;
- reduce expressions modulo $C_2$ without requiring a complete basis of `V_h`;
- compute quotient normal forms or obstruction witnesses;
- solve for explicit $C_2$ corrections when possible.

### Layer C. Realization backends

This layer is optional.

If the user provides an explicit realization into free fields, the backend may:

- realize abstract generators as free-field expressions;
- compute sparse expansions in the ambient free-field algebra;
- produce fast quotient normal forms in the ambient free-field $C_2$ quotient.

### Layer D. High-level null search

This layer orchestrates:

1. target construction,
2. quotient precheck,
3. candidate descendant/$C_2$ witness generation,
4. singular/null constraints,
5. explicit null-state reconstruction.

## Core Public Abstractions

### 1. `LocalOperatorCanonicalizer`

Purpose: provide generic canonicalization and sparse local coordinates for arbitrary operator expressions.

Responsibilities:

- canonicalize expressions using `simplify(...)` and canonical normal-order rules;
- extract sparse coefficient maps from canonical expressions;
- optionally expose grading/sector data used for block decomposition.

Expected interface sketch:

```python
class LocalOperatorCanonicalizer:
    def __init__(self, generators, stress_tensor=None, gradings=None): ...
    def canonicalize(self, expr): ...
    def sparse_terms(self, expr) -> dict[Any, sp.Expr]: ...
    def sector_of(self, expr) -> Any: ...
    def nested_stress_tensor(self, n): ...
```

Remarks:

- `sparse_terms(expr)` is the fundamental "coordinate" API.
- It must not call a complete `basis(weight)` enumerator.

### 2. `SparseLinearContext`

Purpose: provide on-demand linear algebra over sparse canonical supports.

Responsibilities:

- maintain an incremental sparse elimination state;
- test linear independence of candidate expressions;
- extract zero relations among a finite family of expressions;
- assemble temporary matrices only from the supports currently present.

Expected interface sketch:

```python
class SparseLinearContext:
    def __init__(self, canonicalizer): ...
    def reduce_vector(self, sparse_terms): ...
    def insert_expr(self, expr) -> bool: ...
    def independent_subset(self, expressions): ...
    def zero_relations(self, expressions): ...
```

This layer should reuse and generalize the recent sparse relation-solver infrastructure already added to `operator_spaces.py`.

### 3. `DescendantGenerator`

Purpose: generate descendant candidates from source operators without relying on ambient basis enumeration.

Responsibilities:

- accept one or more source operators, typically singular vectors found at lower weight;
- generate descendants by derivatives and normal-ordered products with chosen generators;
- deduplicate via canonicalization;
- return sparse-independent spanning candidates.

Expected interface sketch:

```python
class DescendantGenerator:
    def __init__(self, canonicalizer, linear_context): ...
    def generate(self, source, target_weight): ...
    def span(self, sources, target_weight): ...
```

Important design decision:

- Do **not** split this into separate Virasoro-only and generic-VOA implementations.
- Virasoro is just the one-generator case.

### 4. `SingularConstraintSystem`

Purpose: encode null/singularity constraints from OPE pole conditions in a reusable sparse form.

Responsibilities:

- derive positive-mode annihilation constraints from OPEs;
- support the stress-tensor case as one instance of the general rule;
- build linear systems for ansaetze without requiring a complete basis of the ambient space.

Expected interface sketch:

```python
class SingularConstraintSystem:
    def __init__(self, canonicalizer, generators=None, stress_tensor=None): ...
    def positive_mode_constraints(self, expr, generators=None): ...
    def equations_for_ansatz(self, ansatz_terms, generators=None): ...
    def solve_null_ansatz(self, ansatz_terms, generators=None): ...
```

### 5. `AbstractC2Reducer`

Purpose: define the common interface for all $C_2$ / quotient reducers.

Responsibilities:

- reduce expressions modulo $C_2$;
- test whether an expression is zero in the quotient;
- optionally produce an explicit $C_2$ witness.

Expected interface sketch:

```python
class AbstractC2Reducer(ABC):
    def quotient_normal_form(self, expr): ...
    def is_zero_mod_c2(self, expr) -> bool: ...
    def solve_c2_witness(self, expr): ...
```

### 6. `GenericC2Reducer`

Purpose: realization-free implementation of `AbstractC2Reducer` for general VOAs.

Responsibilities:

- generate $C_2$ candidates of the form `NO(d(a), phi)` on demand;
- reduce modulo the span of those candidates using sparse elimination;
- avoid any requirement to enumerate all operators at fixed weight.

Expected interface sketch:

```python
class GenericC2Reducer(AbstractC2Reducer):
    def __init__(self, canonicalizer, linear_context, generator_provider=None): ...
    def candidate_generators_for_expr(self, expr): ...
    def quotient_normal_form(self, expr): ...
    def is_zero_mod_c2(self, expr) -> bool: ...
    def solve_c2_witness(self, expr): ...
```

Implementation idea:

- start from the sparse support of the target expression;
- generate only those `NO(d(a), phi)` candidates likely to cancel terms currently present;
- enlarge the candidate pool iteratively until the remainder stabilizes or vanishes.

This gives a general VOA algorithm even when no realization is available.

### 7. `RealizationBackend`

Purpose: abstract optional realization-aware acceleration.

Responsibilities:

- realize abstract operators in an ambient algebra;
- expose sparse term extraction in the realized system;
- optionally provide a faster quotient normal form than the generic reducer.

Expected interface sketch:

```python
class RealizationBackend(ABC):
    def realize(self, expr): ...
    def sparse_terms(self, expr): ...
    def quotient_normal_form(self, expr): ...
```

### 8. `FreeFieldC2Reducer`

Purpose: high-performance `AbstractC2Reducer` for VOAs with explicit free-field realizations.

This is the main acceleration path for $bc\beta\gamma$-type examples.

Responsibilities:

- realize target expressions in the free-field algebra;
- use Wick/canonical free-field expansion to obtain sparse monomial data;
- compute quotient normal forms in the free-field $C_2$ quotient;
- expose obstruction terms and candidate lifts back to the abstract VOA layer.

Expected interface sketch:

```python
class FreeFieldC2Reducer(AbstractC2Reducer):
    def __init__(self, canonicalizer, realization_backend): ...
    def quotient_normal_form(self, expr): ...
    def is_zero_mod_c2(self, expr) -> bool: ...
    def solve_c2_witness(self, expr): ...
```

For a $bc\beta\gamma$ realization, the quotient logic should exploit that:

- all total derivatives vanish in the $C_2$ quotient;
- realized normal-ordered products descend to supercommutative multiplication;
- grading data such as weight, charge, ghost number, and parity can be used to block-diagonalize the search.

### 9. `C2NullSearcher`

Purpose: orchestrate the full null search

$$
\mathcal N \in \mathrm{Desc}_h,
\qquad
\mathcal N - T^{(n)} \in C_2(V)_h,
$$

using either the generic reducer or a realization-aware accelerated reducer.

Responsibilities:

- build the target $T^{(n)}`;
- run a quotient precheck via the configured `AbstractC2Reducer`;
- if the quotient class is nonzero, return an obstruction immediately;
- if the quotient class vanishes, build a restricted descendant/$C_2$ ansatz;
- impose singular/null constraints and solve for explicit witnesses.

Expected interface sketch:

```python
class C2NullSearcher:
    def __init__(
        self,
        canonicalizer,
        linear_context,
        descendants,
        singular_constraints,
        c2_reducer,
    ): ...

    def quotient_precheck(self, target_expr): ...
    def search_from_sources(self, target_expr, sources): ...
    def search_stress_tensor_nilpotency(self, n, sources=None): ...
```

## Two-Stage Search Strategy

### Stage 1. Quotient precheck

Given

$$
T^{(n)} = \mathrm{NO}(T, \mathrm{NO}(T, \dots)),
$$

compute its quotient normal form using the configured reducer.

- If the quotient class is nonzero, return an obstruction.
- If the quotient class vanishes, continue.

This stage should be cheap relative to full descendant/null solving.

### Stage 2. Lift and null constraints

Construct a restricted ansatz

$$
\mathcal N = T^{(n)} + \sum_i c_i \psi_i,
\qquad
\psi_i \in C_2,
$$

possibly augmented by descendant generators coming from lower singular sources.

Then impose null/singularity constraints through `SingularConstraintSystem`.

The candidate family should be restricted by:

- conformal weight;
- parity;
- charge / ghost number / other available gradings;
- realization-derived support heuristics when available.

## Generic vs Optimized Reducers

The public `C2NullSearcher` must accept any `AbstractC2Reducer` implementation.

This gives the required layered behavior:

- `GenericC2Reducer`: correctness-first, realization-free, general VOA algorithm;
- `FreeFieldC2Reducer`: realization-aware acceleration for explicit free-field embeddings.

The high-level search flow should be the same in both cases.

## Linear Algebra Philosophy

The old dense formulation

$$
D x - C y = t
$$

is still conceptually correct, but in implementation it must be interpreted sparsely:

- `D`, `C`, and `t` are assembled only from the finite supports appearing in the current problem;
- rows correspond to temporary pivot monomials or temporary quotient-support monomials;
- no step requires precomputing a complete basis of `V_h`.

## Implementation Order

### Phase 1: sparse local infrastructure

1. Promote sparse canonical expansion to a public API.
2. Reuse the existing sparse relation-solver machinery as `SparseLinearContext`.
3. Refactor descendant and singular-vector utilities to consume sparse coordinates rather than ambient-basis coordinates.

### Phase 2: generic $C_2$ / quotient layer

1. Introduce `AbstractC2Reducer`.
2. Implement `GenericC2Reducer` with on-demand `NO(d(a), phi)` generation.
3. Add tests on small generic examples and Virasoro toy cases.

### Phase 3: free-field acceleration

1. Introduce `RealizationBackend` and `FreeFieldC2Reducer`.
2. Implement a $bc\beta\gamma$ quotient-normal-form fast path.
3. Add tests showing agreement between the generic and free-field reducers on overlapping tractable examples.

### Phase 4: null search orchestration

1. Implement `C2NullSearcher` against `AbstractC2Reducer`.
2. Add stress-tensor nilpotency search helpers.
3. Add tests for known Virasoro null relations and free-field examples.

## Short-Term MVP

The first end-to-end milestone should be:

- expose sparse canonical term extraction as a stable public API;
- implement `AbstractC2Reducer` plus a minimal `GenericC2Reducer`;
- implement `C2NullSearcher.quotient_precheck(...)` without any ambient basis enumeration;
- verify on Virasoro that the framework can decide whether

$$
\mathrm{NO}(T,T)
$$

is trivial or nontrivial in the quotient.

Only after this generic path works should the free-field acceleration layer be added.

## Trellis Implementation Checklist

This section translates the design into a concrete implementation sequence for follow-up coding sessions.

### Progress Snapshot (2026-03-23)

- Workstream 0 is largely landed: sparse-term access and `SparseLinearContext` exist, and relation tests now cover the no-ambient-basis path.
- Workstream 1 is landed at MVP level in `src/pyope/c2.py` via `AbstractC2Reducer`, `GenericC2Reducer`, quotient normal forms, and witness-returning diagnostics.
- Workstream 2 is landed and already extends past the original precheck-only milestone: `src/pyope/null_search.py` now supports structured prechecks, stress-tensor search metadata, and a quotient-to-descendant lift path.
- Workstream 3 is partially landed: `src/pyope/descendants.py` and `src/pyope/singularity.py` were extracted, but descendant independence and singular-equation assembly are not yet fully sparse-first.
- Workstream 4 is partially landed: `src/pyope/realizations.py` and `src/pyope/free_field_c2.py` exist, and the first realization rule kills derivative factors in the quotient; richer free-field rules and sector pruning remain pending.
- Workstream 5 is partially landed: the new core can return explicit `null_operator` / `c2_remainder` data on toy examples, but known null-relation recovery and stronger witness reconstruction still remain open.

### Workstream 0: stabilize sparse infrastructure

Goal: make sparse canonical expansion the default internal coordinate language.

Tasks:

- [ ] Audit current sparse utilities in `src/pyope/operator_spaces.py` and identify which helpers should become public or semi-public APIs.
- [ ] Expose a stable sparse-term API, likely via `LocalOperatorCanonicalizer.sparse_terms(expr)` or an equivalent method on the current `LocalOperatorBasis`.
- [ ] Wrap `_SparseIndependentEliminator` and related sparse relation helpers into a reusable `SparseLinearContext` abstraction.
- [ ] Ensure these sparse APIs do not call `basis(weight)` internally.
- [ ] Add tests proving sparse relation solving works without fixed-weight basis enumeration.

Suggested file targets:

- `src/pyope/operator_spaces.py`
- `tests/test_operator_spaces.py`

Acceptance checks:

- [ ] Independence and zero-relation tests pass without ambient basis enumeration.
- [ ] Sparse-term extraction works for both abstract and realized expressions.

### Workstream 1: generic $C_2$ reducer

Goal: provide a realization-free quotient algorithm for general VOAs.

Tasks:

- [ ] Create `src/pyope/c2.py`.
- [ ] Define `AbstractC2Reducer`.
- [ ] Implement a first `GenericC2Reducer` that performs on-demand $C_2$ reduction using sparse elimination.
- [ ] Decide on the first minimal heuristic for `candidate_generators_for_expr(expr)`.
- [ ] Implement `quotient_normal_form(expr)`.
- [ ] Implement `is_zero_mod_c2(expr)` in terms of the quotient normal form.
- [ ] Add a minimal witness-returning path for diagnostics, even if `solve_c2_witness(expr)` is initially partial.

Suggested file targets:

- `src/pyope/c2.py`
- `tests/test_c2.py`

Acceptance checks:

- [ ] Simple Virasoro toy examples reduce correctly.
- [ ] The reducer does not require a complete weight-$h$ basis.
- [ ] Nontrivial quotient obstructions are returned in a readable sparse form.

### Workstream 2: null-search precheck API

Goal: expose the first high-level search entry point before full null-state lifting.

Tasks:

- [ ] Create `src/pyope/null_search.py`.
- [ ] Define a `NullSearchResult` or equivalent structured result type.
- [ ] Implement `C2NullSearcher` against `AbstractC2Reducer`.
- [ ] Implement `quotient_precheck(target_expr)`.
- [ ] Implement `search_stress_tensor_nilpotency(n, sources=None)` in precheck-only form.
- [ ] Return diagnostics including target expression, quotient remainder, and obstruction terms.

Suggested file targets:

- `src/pyope/null_search.py`
- `tests/test_null_search.py`

Acceptance checks:

- [ ] The API can answer whether `T^(n)` is zero or nonzero in `V/C_2(V)`.
- [ ] The result object clearly distinguishes success, obstruction, and "needs lift" states.

### Workstream 3: generic descendant and singularity refactor

Goal: move descendant generation and null constraints onto sparse coordinates.

Tasks:

- [ ] Extract descendant logic into `src/pyope/descendants.py`.
- [ ] Replace ambient-basis coordinate dependence with sparse-linear-context dependence.
- [ ] Extract singular/null OPE constraints into `src/pyope/singularity.py`.
- [ ] Implement `SingularConstraintSystem.equations_for_ansatz(...)` in a basis-free sparse style.
- [ ] Verify the Virasoro-only case remains a one-generator specialization of the generic implementation.

Suggested file targets:

- `src/pyope/descendants.py`
- `src/pyope/singularity.py`
- `tests/test_descendants.py`
- `tests/test_singularity.py`

Acceptance checks:

- [ ] Descendant generation no longer depends on `coordinates(..., weight=...)` over a full basis.
- [ ] Singular-vector searches still work on tractable toy examples.

### Workstream 4: free-field realization backend

Goal: add the high-performance path for VOAs with explicit free-field realizations.

Tasks:

- [ ] Create `src/pyope/realizations.py` with a minimal `RealizationBackend` abstraction.
- [ ] Create `src/pyope/free_field_c2.py`.
- [ ] Implement `FreeFieldC2Reducer` for realization-aware quotient reduction.
- [ ] Add the first $bc\beta\gamma$ quotient-normal-form logic.
- [ ] Encode derivative killing and supercommutative product descent in the free-field quotient path.
- [ ] Add sector pruning hooks: weight, parity, ghost number, charge.
- [ ] Add overlap tests comparing generic and free-field reducers on tractable examples.

Suggested file targets:

- `src/pyope/realizations.py`
- `src/pyope/free_field_c2.py`
- `tests/test_free_field_c2.py`

Acceptance checks:

- [ ] Free-field reduction agrees with the generic reducer where both are feasible.
- [ ] The free-field path gives a noticeably smaller or simpler quotient remainder on $bc\beta\gamma$ examples.

### Workstream 5: explicit null-state lifting

Goal: go beyond quotient precheck and construct explicit null states of the form

$$
\mathcal N = T^{(n)} + \varphi,
\qquad
\varphi \in C_2.
$$

Tasks:

- [ ] Extend `C2NullSearcher` to build restricted descendant/$C_2$ ansaetze.
- [ ] Implement the first usable `solve_c2_witness(expr)` path for the configured reducer.
- [ ] Combine quotient-zero information with singular/null constraints.
- [ ] Return explicit `null_expr`, `c2_correction`, and diagnostic coefficient data.
- [ ] Add tests for at least one known Virasoro null relation.

Suggested file targets:

- `src/pyope/null_search.py`
- `tests/test_null_search.py`

Acceptance checks:

- [ ] The searcher can produce an explicit null witness, not just a quotient obstruction.
- [ ] Known toy-model null relations are recovered in canonicalized local-operator form.

## Recommended Next Coding Session Start Point

The next coding session should start from Workstream 0 and Workstream 1 only.

Recommended immediate scope:

- expose sparse-term extraction publicly;
- introduce `SparseLinearContext`;
- create `src/pyope/c2.py` with `AbstractC2Reducer` and a minimal `GenericC2Reducer`;
- stop before realization-specific acceleration.

This keeps the first implementation session small, testable, and aligned with the generic-first design.

## Notes

- This plan is intentionally API-first.
- It does not yet commit to a final file layout, though `src/pyope/operator_spaces.py` or a small cluster of related modules remains the natural landing zone.
- The existing `src/pyope/null_states.py` module is free-field/Fock oriented and should not dictate the generic local-operator/$C_2$ API design.
- Fixed-weight basis enumeration utilities may remain in the codebase for small positive-weight examples, but they should be treated as optional helpers, not foundational infrastructure for this task.
