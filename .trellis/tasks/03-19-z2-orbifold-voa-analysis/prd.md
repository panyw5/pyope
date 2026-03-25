# Z2 Orbifold VOA Analysis

## Problem

Study the analogue of the current $\mathbb{Z}_3$ orbifold / quotient problem for a corresponding $\mathbb{Z}_2$ quotient of a free-field VOA built from
$$
b,\ \partial c,\ \beta,\ \gamma.
$$

The purpose is twofold:

- understand the $\mathbb{Z}_2$ case on its own terms;
- use it as a comparison case to clarify why the $\mathbb{Z}_3$ case becomes structurally complicated.

## Motivation

The current $\mathbb{Z}_3$ analysis suggests that direct closure of naive invariant monomials under OPE produces complicated mixtures of descendants, improvements, and possible indecomposable $\mathcal N=2$ blocks. It is plausible that the $\mathbb{Z}_2$ analogue is simpler because:

- the lowest invariant chiral objects are quadratic rather than cubic;
- the quantum numbers are closer to familiar rank-one small $\mathcal N=4$ / $\mathcal N=2$ examples;
- the associated family/block structure may be cleaner and reveal the correct methodology for the $\mathbb{Z}_3$ problem.

## Key Questions

1. What is the natural $\mathbb{Z}_2$ action on the free fields, and what is the resulting invariant VOA?
2. Does the invariant VOA contain a natural $\mathcal N=2$ SCA subalgebra of the same form as in the $\mathbb{Z}_3$ case?
3. What are the lowest invariant generators?
4. Does the invariant algebra close cleanly in terms of a small finite list of $\mathcal N=2$ families?
5. Is the $\mathbb{Z}_2$ case closer to a standard small $\mathcal N=4$-type situation, or still only naturally expressible in $\mathcal N=2$ family/block language?

## Working Expectations

The $\mathbb{Z}_2$ case is expected to be simpler than $\mathbb{Z}_3$ because the obvious invariant fields start at lower degree. Natural candidates include objects such as
$$
J = :\beta\gamma:,
\qquad
W_+ = :\beta^2:,
\qquad
W_- = :\gamma^2:,
$$
and fermionic mixed invariants such as
$$
:\beta b:,
\qquad
:\gamma\partial c:.
$$

This suggests that the $\mathbb{Z}_2$ problem may admit a more transparent organization into finitely many $\mathcal N=2$ families, perhaps even a familiar enhanced structure.

## Proposed Plan

### 1. Define the precise $\mathbb{Z}_2$ action

- Make explicit how $\mathbb{Z}_2$ acts on
  $$
  b,\ c,\ \beta,\ \gamma.
  $$
- Confirm which fields and composites are invariant.

### 2. Identify the natural $\mathcal N=2$ subalgebra

- Check whether the standard fields
  $$
  J=:\beta\gamma:,
  \qquad
  G=:b\gamma:,
  \qquad
  \bar G=:\beta\partial c:,
  \qquad
  T=\frac12(:\partial\beta\,\gamma:-:\beta\partial\gamma:)-:b\partial c:
  $$
  remain invariant and still generate a $\mathcal N=2$ SCA.

### 3. Low-weight invariant search

- Enumerate the lowest-weight $\mathbb{Z}_2$-invariant fields.
- Organize them by $(h,q)$ under the $\mathcal N=2$ subalgebra.

### 4. Sector-by-sector primary search

- In each low-weight / low-charge sector, solve for $\mathcal N=2$ superconformal primaries using `OPEdefs.m`.
- Determine whether the invariant algebra is best described by:
  - clean superprimary families; or
  - mixed indecomposable blocks.

### 5. Rewrite key OPEs in family language

- Focus first on the lowest nontrivial chiral/anti-chiral pair.
- Attempt to express the analogue of
  $$
  W\times \bar W
  $$
  in the form
  $$
  [\mathbf 1] + [J] + [\mathcal U] + \cdots
  $$
  or the appropriate simpler structure.

### 6. Compare against the $\mathbb{Z}_3$ case

- Record which features simplify:
  - number of candidate generators,
  - existence/non-existence of neutral superprimaries,
  - amount of block mixing,
  - closure of low-weight OPEs.
- Use this comparison to refine the methodology for the $\mathbb{Z}_3$ task.

## Deliverables

- A low-weight description of the $\mathbb{Z}_2$ orbifold VOA.
- A candidate list of $\mathcal N=2$ superprimary families or blocks.
- A comparison note explaining whether and why the $\mathbb{Z}_2$ case is simpler than the $\mathbb{Z}_3$ case.

## Practical Notes

- Use `OPEdefs.m` as the authoritative computational backend.
- Prefer $\mathcal N=2$ family/block organization over naive monomial closure.
- Keep this task separate from the $\mathbb{Z}_3$ task, but record cross-references as patterns emerge.
