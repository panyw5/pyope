# Z3 W-algebra for p=3

## Problem

Study section 5.1 of
`papers/[Bonetti, Meneghelli, Rastelli] VOAs labelled by complex reflection groups and 4d SCFTs.md`
for the rank-one complex reflection group
$$
G = \mathbb{Z}_3
$$
with
$$
p=3.
$$

The goal is to push the $\mathcal N=2$ W-algebra analysis beyond the currently established low-weight evidence and obtain a sharper structural understanding of the $p=3$ algebra.

## Primary Goals

- Construct the complete anti-chiral sector generator $\bar W$ and clarify its full $\mathcal N=2$ multiplet structure.
- Check that the proposed strong generators satisfy Jacobi identities and therefore define a consistent VOA presentation.
- Verify that the resulting strong-generator presentation reproduces the expected flavor Schur index.
- Construct the conformal-weight $8$ null state
  $$
  \mathcal N_T = L_{-2}^4\lvert 0\rangle + \varphi,
  \qquad
  \varphi \in C_2.
  $$

## Scope

This task is specifically about the $G=\mathbb Z_3$ example in section 5.1 and its $p=3$ W-algebra structure. It includes:

- free-field and OPE analysis of the candidate strong generators;
- organization into $\mathcal N=2$ multiplets;
- Jacobi/associativity consistency checks;
- index comparison against the expected Schur-side counting;
- $C_2$-aware null-state construction at weight $8$.

This task does not require a full general-$p$ classification unless it directly helps the $p=3$ case.

## Working Questions

### 1. Strong generators

- What is the minimal strong generating set for the $p=3$ algebra?
- Is
  $$
  \{T,J,G,\bar G,W,\bar W\}
  $$
  sufficient, or are additional neutral/non-chiral generators required?
- What is the cleanest expression for $\bar W$ in the chosen realization/basis?

### 2. Jacobi consistency

- Do the singular OPEs among the proposed strong generators close modulo descendants and null states?
- Which Jacobi identities are nontrivial in the $p=3$ case?
- Are there coefficient choices uniquely fixed by Jacobi constraints?

### 3. Flavor Schur index

- What refined counting should be reproduced from the proposed strong-generator algebra?
- Does the vacuum character / filtered counting from the algebra agree with the expected flavor Schur index at low orders?
- If there is a mismatch, does it signal a missing generator or a null relation?

### 4. Weight-8 null state

- Can one explicitly solve for
  $$
  \varphi \in C_2
  $$
  such that
  $$
  \mathcal N_T = L_{-2}^4\lvert 0\rangle + \varphi
  $$
  is null?
- Is this null state implied by Jacobi consistency, by the free-field realization, or by index constraints?

## Acceptance Criteria

- [ ] A concrete candidate for the complete $\bar W$ is written down and normalized.
- [ ] A proposed minimal strong generating set is identified for the $p=3$ algebra.
- [ ] Nontrivial Jacobi identities for the strong generators are checked and documented.
- [ ] Low-order flavor Schur index agreement is verified, or any mismatch is isolated and explained.
- [ ] A weight-$8$ null state of the form
  $$
  \mathcal N_T = L_{-2}^4\lvert 0\rangle + \varphi,
  \qquad
  \varphi \in C_2
  $$
  is explicitly constructed or reduced to a well-posed linear-algebra problem with a reproducible search setup.

## Suggested Plan

### Phase 1: Reconstruct algebra data

- Re-read section 5.1 and collect the exact $p=3$ generator/OPE ansatz.
- Match this against the existing free-field realization and previous $\mathbb Z_3$ notes.
- Fix conventions for $T,J,G,\bar G,W,\bar W$ and charges/weights.

### Phase 2: Strong-generation and closure

- Compute key OPEs involving $W$ and $\bar W$.
- Identify whether additional top states appear after quotienting descendants and obvious composites.
- Decide the minimal strong generating set.

### Phase 3: Jacobi analysis

- Test the Jacobi identities for the proposed generators.
- Record which OPE coefficients are forced.
- Determine whether null relations are needed for consistency.

### Phase 4: Index check

- Compute low-order counting from the candidate algebra.
- Compare with the expected flavor Schur index for the $G=\mathbb Z_3$ theory.

### Phase 5: Weight-8 null search

- Build the weight-$8$ descendant space of the vacuum.
- Build the weight-$8$ $C_2$ space.
- Solve for
  $$
  L_{-2}^4\lvert 0\rangle + \varphi
  $$
  with
  $$
  \varphi \in C_2
  $$
  that becomes null.

## Notes

- Previous related Trellis task: `.trellis/tasks/03-19-z3-orbifold-voa-analysis/prd.md`.
- The earlier task provides useful low-weight evidence, but this task is more focused on the section 5.1 $\mathbb Z_3$/$p=3$ W-algebra bootstrap and the specific weight-$8$ null state.
