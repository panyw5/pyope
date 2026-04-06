# Harden Wolfram backend expression protocol

## Goal

Make the `wolfram` backend reliable for notebook-scale `NO(...)` / `OPE(...)` computations in `pyope`, especially for nested normal-ordered expressions appearing in `demo/N_equals_three_algebra.ipynb`.

## What we already know

* `wolframscript` is installed and callable from the optimization worktree.
* The initial path issue for `OPEdefs/OPEdefs.wls` from `demo/` cwd has been fixed.
* The bridge previously emitted malformed derivative payloads and has been partially corrected.
* The bridge now tolerates `sympy.Float`, but notebook-derived computations still expose deeper protocol issues.
* Current failure indicates the Wolfram payload can decode into illegal operator multiplication on the Python side.
* The notebook null-state mismatch was traced to definition drift versus the reference `OPEdefs.m` workflow, not to a raw Wolfram arithmetic error. In particular, `G = NO(b, γ)` is required for the representative `N=3` null state to vanish.

## Observed failure modes

* Wrapper path depended on notebook cwd instead of repository-relative absolute path.
* Derivative payloads were encoded in malformed forms instead of stable `dn(order, operator)` syntax.
* `sympy.Float` values leaked into `_encode_expr(...)` and caused unsupported-expression failures.
* Decoder currently relies on `eval(...)` over a text protocol that still allows invalid operator `Times` structures.
* Pretty-print style objects such as `∂^2 b` must never appear in the machine protocol.

## Requirements

* The Wolfram wrapper output MUST be a strict machine-readable grammar, not a display-oriented rendering.
* Payloads MUST only use Python-side supported constructors such as `MakeOPE(...)`, `NO(...)`, `dn(...)`, addition, and scalar multiplication.
* The bridge MUST NOT emit illegal operator-times-operator-expression products.
* The protocol MUST preserve exact arithmetic whenever possible.
* Notebook-scale expressions from `demo/N_equals_three_algebra.ipynb` MUST round-trip without protocol-level failures.

## Acceptance Criteria

* [ ] `demo/N_equals_three_algebra.ipynb` representative nested `NO(...)` expressions no longer fail on bridge protocol errors.
* [ ] No bridge payload contains pretty-print derivative fragments such as `∂^n b`.
* [ ] No bridge payload decodes into direct operator multiplication forbidden by `pyope` semantics.
* [ ] Focused regression tests cover the current notebook failure shapes.
* [ ] Existing backend tests continue to pass.

## Progress Update (2026-04-06)

### Correctness progress

* The bridge now uses a stricter machine protocol in `src/pyope/wolfram_backend.py` and `OPEdefs/OPEdefs.wls`.
* Python-side text `eval(...)` decoding has been replaced by a constrained AST decoder.
* Exact rationals are preserved across the bridge instead of leaking into floats.
* Notebook-derived regression coverage has been expanded in `tests/test_backend.py`.
* The representative null state from `demo/N_equals_three_algebra.ipynb` is confirmed to simplify to `Zero` through the Wolfram bridge **when the notebook definitions match the reference `OPEdefs.m` setup**.

### Performance findings

* Profiling shows the dominant cost is repeated `wolframscript` process startup, not the underlying `OPEdefs` computation.
* A representative null-state construction under `set_compute_backend("wolfram")` triggered `63` separate Wolfram launches, costing about `57s` total.
* A whole-expression prototype was added:
  * `OPEdefs/OPEdefs.wls` now supports an `EVAL` mode.
  * `src/pyope/wolfram_backend.py` now exposes `wolfram_evaluate_expr(...)` and `simplify_with_wolfram(...)`.
* When the expression is first built structurally under the `sympy` backend and only then sent wholesale to Wolfram, the same representative null-state simplification drops to about `1.4s` and returns `Zero`.

### Current limitation

* `set_compute_backend("wolfram")` still eagerly dispatches every binary `NO/OPE` during expression construction, so it remains slow for large notebook-scale expressions.
* The likely next architectural step is a persistent Wolfram session with registry synchronization, so `set_compute_backend("wolfram")` can remain automatic without paying a process-launch cost on every primitive call.

## Next implementation plan

1. Audit current `PyEncode[...]` output forms in `OPEdefs/OPEdefs.wls`.
2. Define a minimal closed payload grammar for the bridge.
3. Refactor wrapper encoding so non-scalar multiplication is normalized before serialization.
4. Tighten Python-side decoder assumptions to align with the grammar.
5. Add regression tests using notebook-derived nested `NO` expressions.

## Out of Scope

* Long-lived Wolfram kernel sessions.
* Benchmark-driven backend auto-selection.
* Full simplification parity across every internal helper.

## Technical Notes

* Primary user-facing reproduction is in `demo/N_equals_three_algebra.ipynb`.
* Current bridge implementation lives in `src/pyope/wolfram_backend.py` and `OPEdefs/OPEdefs.wls`.
* Existing backend tests cover basic path selection, binary `OPE`, binary `NO`, and linear `NO` structure, but not notebook-scale nested payloads.
