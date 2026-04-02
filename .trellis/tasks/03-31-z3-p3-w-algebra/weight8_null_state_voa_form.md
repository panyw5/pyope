# Weight-8 Null State In VOA Notation

## Generator dictionary

In the internal computation labels,

- `SG0 = T`
- `SG1 = J`
- `SG2 = G`
- `SG3 = \bar G`
- `SG4 = W`
- `SG5 = \bar W`

The relevant conformal weights are

- $h_T = 2$
- $h_J = 1$
- $h_G = h_{\bar G} = 3/2$
- $h_W = h_{\bar W} = 3/2$

## Free-field realizations actually used

The computation uses the rank-one $\beta\gamma bc$ system with free fields

- $b$ fermionic of weight $2$
- $c$ fermionic of weight $-1$
- $\beta$ bosonic of weight $3/2$
- $\gamma$ bosonic of weight $-1/2$

and OPEs
$$
b(z)c(w) \sim \frac{1}{z-w},
\qquad
\beta(z)\gamma(w) \sim -\frac{1}{z-w}.
$$

In these conventions, the strong generators used in the weight-$8$ search are

$$
J = 2:bc: + 3:\beta\gamma:,
$$

$$
G = :\gamma b:,
$$

$$
\bar G = 2:(\partial\beta)c: + 3:\beta(\partial c):,
$$

$$
T = -2:b\,\partial c: - \frac{3}{2}:\beta\,\partial\gamma: - :(\partial b)c: - \frac{1}{2}:(\partial\beta)\gamma:,
$$

$$
W = \beta,
$$

$$
\begin{aligned}
\bar W
= {}& :\beta\beta\gamma\gamma\gamma:
+ 2:\beta\gamma\gamma bc:
- 4:\beta(\partial\gamma)\gamma:
- \frac{4}{3}:\gamma b\partial c: \\
&+ \frac{2}{3}:\gamma(\partial b)c:
+ \frac{2}{3}:(\partial\beta)\gamma\gamma:
- \frac{8}{3}:(\partial\gamma)bc:
+ \frac{10}{9}\partial^2\gamma.
\end{aligned}
$$

For reference, the two additional singular-OPE descendants that appear in the strong-generation closure are

$$
GW := \{G,W\}_1 = b,
$$

and

$$
\begin{aligned}
GWbar = {}& -\frac{10}{9}\partial^3 c
- \frac{8}{3}:b(\partial^2 c)c:
- 3:\beta\beta\gamma\gamma\partial c:
+ 4:\beta\gamma b(\partial c)c: \\
&+ 4:\beta\gamma\partial^2 c:
+ 4:\beta(\partial\gamma)(\partial c):
- \frac{2}{3}:(\partial^2\beta)\gamma c:
+ \frac{2}{3}:(\partial b)(\partial c)c: \\
&- 2:(\partial\beta)\beta\gamma\gamma c:
+ \frac{8}{3}:(\partial\beta)(\partial\gamma)c:.
\end{aligned}
$$

## Definition code locations

The same free-field realization is implemented in several task-local scripts.

### Python definitions

The main reusable Python definition block is in
` .trellis/tasks/03-31-z3-p3-w-algebra/character_scan_p3.py `,
function `build_p3_data()`:

- free fields and OPE setup: lines `54-62`
- $J,G,\bar G,T,W,\bar W$: lines `64-83`
- $GW, GWbar$: lines `84-85`

The same expressions are also used or restated in:

- ` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_p3.py `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/compare_nt_terms_pyope.py `

### OPEdefs / Mathematica definitions

The corresponding `OPEdefs.m` definitions used in the final weight-$8$ verification are in
` .trellis/tasks/03-31-z3-p3-w-algebra/verify_weight8_null_relation_opedefs_fixed.wls `:

- free-field OPE setup: lines `7-11`
- $J,G,\bar G,T,W,GW,\bar W,GWbar$: lines `13-40`

The batch `OPEdefs` realization pipelines use the same data indirectly via exported generator definitions in:

- ` .trellis/tasks/03-31-z3-p3-w-algebra/export_weight8_descendants_for_opedefs.py `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/export_weight8_candidates_for_opedefs.py `

## Field-form relation

The weight-$8$ relation found in the pure `OPEdefs.m` realization is

$$
\begin{aligned}
0 = \mathcal N_T^{\mathrm{field}}
&:= T^4
-\frac{2}{3}T^3(JJ)
+\frac{1}{3}T^3\partial J
-\frac{2}{3}T^2\bigl(J(G\bar G)\bigr)
+\frac{1}{2}T^2(G\partial\bar G) \\
&\phantom{:} -\frac{9}{4}T^2(W\partial\bar W)
-\frac{1}{12}T^2(\partial J\,JJ)
-\frac{1}{3}T^2(\partial J\,\partial J)
-\frac{1}{6}T^2(\partial G\,\bar G) \\
&\phantom{:} +\frac{15}{4}T^2((\partial W)\bar W)
+\frac{3}{4}T^2(\partial^2J\,J).
\end{aligned}
$$

Here products are nested normal-ordered products, with the same bracketing as in the computation.

## Vacuum-state form

Applying state-field correspondence to the above gives a natural vacuum descendant representative

$$
\begin{aligned}
0 = \mathcal N_T
&= L_{-2}^4\lvert 0\rangle
-\frac{2}{3}L_{-2}^3 J_{-1}^2\lvert 0\rangle
+\frac{1}{3}L_{-2}^3 J_{-2}\lvert 0\rangle \\
&\phantom{=} -\frac{2}{3}L_{-2}^2 J_{-1} G_{-3/2}\bar G_{-3/2}\lvert 0\rangle
+\frac{1}{2}L_{-2}^2 G_{-3/2}\bar G_{-5/2}\lvert 0\rangle \\
&\phantom{=} -\frac{9}{4}L_{-2}^2 W_{-3/2}\bar W_{-5/2}\lvert 0\rangle
-\frac{1}{12}L_{-2}^2 J_{-2}J_{-1}^2\lvert 0\rangle
-\frac{1}{3}L_{-2}^2 J_{-2}^2\lvert 0\rangle \\
&\phantom{=} -\frac{1}{6}L_{-2}^2 G_{-5/2}\bar G_{-3/2}\lvert 0\rangle
+\frac{15}{4}L_{-2}^2 W_{-5/2}\bar W_{-3/2}\lvert 0\rangle
+\frac{3}{4}L_{-2}^2 J_{-3}J_{-1}\lvert 0\rangle.
\end{aligned}
$$

This is not yet a quasiprimary basis. It is the direct descendant-basis translation of the normal-ordered field relation.

## Remarks

1. The coefficients above are exactly the coefficients extracted from the pure `OPEdefs.m` linear dependence.
2. The state form should be interpreted as the straightforward vacuum-state representative of the nested normal-order relation, not yet as the unique reduced quasiprimary presentation.
3. A next cleanup step would be to rewrite this in a Virasoro-quasiprimary basis at weight $8$.

## First-pass quasiprimary organization

This section records the first step toward a cleaner quasiprimary-basis rewrite. The goal here is not yet to produce the fully reduced weight-$8$ quasiprimary basis, but to separate the obviously quasiprimary composite blocks from the terms that are manifestly Virasoro-descendant noise.

### Basic quasiprimary blocks

The codebase already implements the quasiprimary completion
$$
(AB)_n = \mathrm{qp}(A,B,n)
$$
in
` src/pyope/quasiprimary.py `
via `quasiprimary_product(...)` / `qp(...)`.

For the present algebra, the most useful low-weight quasiprimary composites are:

$$
\Lambda_{JJ} := (JJ)_0 = \mathrm{qp}(J,J,0),
$$

$$
\Lambda_{G\bar G} := (G\bar G)_0 = \mathrm{qp}(G,\bar G,0),
$$

$$
\Lambda_{JT} := (JT)_0 = \mathrm{qp}(J,T,0),
$$

$$
\Lambda_{W\bar W} := (W\bar W)_0 = \mathrm{qp}(W,\bar W,0).
$$

In the $p=3$ algebra, the Higgs null says precisely that

$$
\Lambda_{W\bar W}
= \frac{1}{27}(J\Lambda_{JJ})_0
- \frac{2}{9}\Lambda_{G\bar G}
- \frac{4}{9}\Lambda_{JT}.
$$

This is the main reason the $W\bar W$ sector can eventually be traded for $J$-, $G\bar G$-, and $JT$-type quasiprimary data.

### Weight-$8$ relation, grouped by structure

The current $11$-term null can be grouped as

$$
\mathcal N_T = \mathcal N_T^{\mathrm{core}} + \mathcal N_T^{\mathrm{desc}},
$$

where the structurally central part is

$$
\mathcal N_T^{\mathrm{core}}
= T^4
- \frac{2}{3}T^3\Lambda_{JJ}
- \frac{2}{3}T^2\bigl(J\,\Lambda_{G\bar G}\bigr)
$$

together with the $W\bar W$-sector contribution

$$
- \frac{9}{4}T^2(W\partial\bar W)
+ \frac{15}{4}T^2((\partial W)\bar W),
$$

and the more obviously descendant-looking remainder

$$
\begin{aligned}
\mathcal N_T^{\mathrm{desc}}
= {}& \frac{1}{3}T^3\partial J
+ \frac{1}{2}T^2(G\partial\bar G)
- \frac{1}{12}T^2(\partial J\,JJ)
- \frac{1}{3}T^2(\partial J\,\partial J) \\
&- \frac{1}{6}T^2(\partial G\,\bar G)
+ \frac{3}{4}T^2(\partial^2J\,J).
\end{aligned}
$$

This split is not unique, but it is useful because:

1. $T^4$, $T^3\Lambda_{JJ}$, and $T^2(J\Lambda_{G\bar G})$ are the pieces most closely tied to genuine weight-$8$ quasiprimary structure;
2. the remaining terms visibly contain explicit derivatives and are therefore the first pieces one expects to reshuffle when moving from a descendant basis to a quasiprimary basis.

### What still needs to be done

To obtain the fully cleaned-up quasiprimary-basis presentation, one still needs to:

1. replace the $W\partial\bar W$ and $(\partial W)\bar W$ sector using the Higgs-null family and its derivatives;
2. identify the precise quasiprimary completions of the weight-$4$, weight-$5$, and weight-$6$ composite blocks appearing inside the nested products;
3. subtract the pure $L_{-1}$ / Virasoro-descendant contributions so that the final expression is written in a basis of genuine weight-$8$ quasiprimary operators plus manifest descendants.

So at the present stage, this document contains a **first-pass quasiprimary organization**, not yet the final reduced quasiprimary basis.

## Refinement of the $W\partial\bar W$ sector

The next useful simplification is to separate the symmetric and antisymmetric combinations of

$$
A := :W\,\partial\bar W:,
\qquad
B := :(\partial W)\bar W:.
$$

### Symmetric part

By the Leibniz rule for normal-ordered products,

$$
A + B = \partial\, :W\bar W:.
$$

So the symmetric part is a pure total derivative. In particular, it is not a genuinely new weight-$4$ quasiprimary structure.

### Antisymmetric part

The combination appearing in the weight-$8$ null is

$$
-\frac{9}{4}A + \frac{15}{4}B.
$$

Writing this in terms of

$$
S := A+B,
\qquad
D := B-A,
$$

gives the exact identity

$$
-\frac{9}{4}A + \frac{15}{4}B
= \frac{3}{4}S + 3D
= \frac{3}{4}\partial:W\bar W: + 3\bigl((\partial W)\bar W - W\partial\bar W\bigr).
$$

So the $W\bar W$ derivative sector splits into:

1. a total-derivative piece
   $$
   \frac{3}{4}\partial:W\bar W:,
   $$
2. a genuinely nontrivial antisymmetric piece
   $$
   3\bigl((\partial W)\bar W - W\partial\bar W\bigr).
   $$

### Quasiprimary interpretation

The first term is immediately tied to the Higgs-family quasiprimary block

$$
\Lambda_{W\bar W} = (W\bar W)_0,
$$

because $:W\bar W:$ differs from $\Lambda_{W\bar W}$ only by lower-descendant corrections. Hence

$$
\partial:W\bar W:
$$

belongs to the descendant sector built from the same Higgs multiplet.

The second term,

$$
\Omega_{W\bar W}
:= (\partial W)\bar W - W\partial\bar W,
$$

is the genuinely new part that still needs to be rewritten in a cleaner quasiprimary basis. So after this refinement, the problematic $W\bar W$ contribution is no longer the whole pair of terms, but only the antisymmetric combination $\Omega_{W\bar W}$.

### Updated structural grouping

With this notation, the weight-$8$ null may be reorganized one step further as

$$
\mathcal N_T = \widetilde{\mathcal N}_T^{\mathrm{core}} + \widetilde{\mathcal N}_T^{\mathrm{desc}} + 3T^2\Omega_{W\bar W},
$$

where

$$
\widetilde{\mathcal N}_T^{\mathrm{core}}
= T^4
- \frac{2}{3}T^3\Lambda_{JJ}
- \frac{2}{3}T^2\bigl(J\Lambda_{G\bar G}\bigr),
$$

and the enlarged descendant-looking part now contains

$$
\frac{3}{4}T^2\partial:W\bar W:
$$

instead of the original two separate $W\partial\bar W$ terms.

This is cleaner because one visibly isolated the only genuinely new $W\bar W$-derivative direction from the obvious total-derivative contamination.

## Attempt to split $\Omega_{W\bar W}$ into quasiprimary plus descendants

I then tested whether the antisymmetric block

$$
\Omega_{W\bar W} = (\partial W)\bar W - W\partial\bar W
$$

can be reduced by subtracting the most natural derivative descendants built from the already identified low-weight quasiprimary composites:

$$
\partial\Lambda_{W\bar W},
\qquad
\partial\Lambda_{G\bar G},
\qquad
\partial\Lambda_{JT},
\qquad
\partial^2\Lambda_{JJ},
\qquad
\partial^3J.
$$

Concretely, I asked whether there exist coefficients such that

$$
\Omega_{W\bar W}
+ a_1\partial\Lambda_{W\bar W}
+ a_2\partial\Lambda_{G\bar G}
+ a_3\partial\Lambda_{JT}
+ a_4\partial^2\Lambda_{JJ}
+ a_5\partial^3J
$$

has no fourth- or third-order poles in its OPE with $T$, which is the necessary condition for it to become quasiprimary of weight $4$.

The resulting linear system has **no solution**.

So the present evidence says:

1. $\Omega_{W\bar W}$ is not itself quasiprimary;
2. it also cannot be made quasiprimary by subtracting only the most obvious low-weight descendant combinations listed above;
3. a full quasiprimary-basis rewrite will therefore require a larger mixing space, likely involving additional weight-$4$ composite directions from the $W\bar W$ sector itself, not just descendants of the previously isolated Higgs-family blocks.

In particular, the $W\bar W$ derivative sector is subtler than a simple
$$
\text{quasiprimary} + \partial(\text{known block})
$$
decomposition.
