# Independent recomputation summary for $\mathcal{W}_{\mathbb{Z}_3}$ at $p=3$

## Scope of this rerun

This note records an independent rerun of the core checks requested in
`.trellis/tasks/03-31-z3-p3-w-algebra/prd.md`, without relying on earlier journal conclusions.

The executable entry point is
` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_p3.py `,
and its latest machine-readable output is
` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_p3_results.json `.

## What was recomputed

1. The free-field candidate for $\bar W$ in the rank-one $\beta\gamma bc$ realization.
2. The paper formulas from section 5.1 for the $W\bar W$ identity family.
3. The $p=3$ Higgs null and the two fermionic nulls.
4. A sample of nontrivial Jacobi identities for the candidate strong generators
   $\{T,J,G,\bar G,W,\bar W\}$.
5. The weight-$8$ stress-tensor candidate $T^4$ modulo $C_2$.

## Results

### 1. Paper formulas reproduced

The rerun reproduces all of the following statements directly from the free-field realization:

- $g_{W\bar W} = -20/9$.
- The $J$ contribution in $W\times \bar W$ matches the section 5.1 formula.
- The quasiprimary level-1 term matches
  $$
  \frac{2}{3}T - \frac{1}{3}(JJ)_0.
  $$
- The Higgs null
  $$
  \mathfrak X^{W\bar W}_{3,0}
  = (W\bar W)_0 - \frac{1}{27}(J(JJ)_0)_0 + \frac{2}{9}(G\bar G)_0 + \frac{4}{9}(JT)_0
  $$
  evaluates to zero.
- The fermionic nulls also evaluate to zero.
- The induced relations
  $$
  \{\mathcal G_W\,\mathfrak X^{W\bar W}_{3,0}\}_2 = \frac{2}{3}\mathfrak X^{GW}_{3,-1},
  \qquad
  \{\widetilde{\mathcal G}_{\bar W}\,\mathfrak X^{W\bar W}_{3,0}\}_2 = \frac{4}{3}\mathfrak X^{\widetilde{\mathcal G}\bar W}_{3,+1}
  $$
  are reproduced.

### 2. Jacobi checks passed on the sampled generator set

The following triples were checked with `check_jacobi_identity(...)` and returned zero in all sampled matrix entries:

- $(T,T,T)$
- $(T,J,J)$
- $(T,G,\bar G)$
- $(J,G,\bar G)$
- $(T,W,\bar W)$
- $(J,W,\bar W)$
- $(G,W,\bar W)$
- $(\bar G,W,\bar W)$
- $(W,W,\bar W)$
- $(W,\bar W,\bar W)$

This is not yet a complete bootstrap proof, but it is consistent with the working strong-generator set
$$
\{T,J,G,\bar G,W,\bar W\}.
$$

### 3. Weight-$8$ stress-tensor null is not supported by this rerun

Two independent checks both point against a null state of the form
$$
\mathcal N_T = L_{-2}^4\lvert 0\rangle + \varphi,
\qquad
\varphi \in C_2.
$$

First, the abstract sparse $C_2$ quotient precheck reports

- status: `obstructed`
- quotient remainder: `NO(T,NO(T,NO(T,T)))`
- $C_2$ part: `Zero`

So in the current implementation the weight-$8$ target does **not** reduce to zero modulo the generated $C_2$ span.

Second, a separate free-field syntactic exclusion test on $\rho(T^4)$ finds

- `target_realizes_to_zero = false`
- `all_terms_in_ff_c2 = false`
- `non_c2_term_count = 27`

This means the realized expression contains explicit non-$C_2(F)$ monomials, so there is no immediate free-field evidence for
$$
[T^4] = 0
$$
in the Zhu $C_2$ quotient.

## Current interpretation

At this point the independent rerun supports the following:

- the free-field realization for $p=3$ does reproduce the section 5.1 null relations and the expected low-order OPE data;
- the candidate strong generators remain
  $$
  \{T,J,G,\bar G,W,\bar W\};
  $$
- the specific weight-$8$ stress-tensor null requested in the PRD is **not** supported by the current independent recomputation.

This does not rule out all possible weight-$8$ null states in the algebra, but it does argue against the particular pure-stress-tensor-leading candidate
$$
L_{-2}^4\lvert 0\rangle + C_2.
$$

## Update: OPEdefs-based weight-$8$ search

The later stages of the investigation moved the expensive free-field realization work from
`pyope` to `OPEdefs.m`, because direct benchmarks showed that `OPEdefs.m` is dramatically faster for the large nested normal-ordered expressions appearing at weight $8$.

### Performance comparison

For a representative slow descendant expression from the weight-$8$ search:

- `pyope realize_and_simplify(...)` took about `609.87` seconds;
- `OPEdefs.m` `Expand[...]` took about `0.0214` seconds.

So for this workload the practical strategy is:

1. use `pyope` only for formal generator / descendant enumeration;
2. use `OPEdefs.m` for actual free-field realization and term expansion;
3. do the final linear algebra on cached sparse term tables.

### OPEdefs descendant cache

Using the three known singular vectors from the paper,
$$
X_h,\qquad X_-,\qquad X_+,
$$
I exported all their weight-$8$ ideal descendants to Mathematica syntax and realized them using `OPEdefs.m`.

The resulting cache contains `128` realized descendants.

### OPEdefs candidate cache

I also exported the full set of weight-$8$, bosonic, $J=0$, non-obvious-$C_2$ candidates to `OPEdefs.m` and realized all of them on the same free-field side.

The resulting cache contains `169` realized candidates.

This removed the earlier ambiguity from mixing `pyope`-side and `OPEdefs`-side realizations in the same linear-algebra step.

## Update: quotient result and explicit relation

Using only `OPEdefs`-realized vectors, the weight-$8$ candidate space and known-singular descendant space satisfy:

- candidate count: `169`
- candidate rank: `140`
- known descendant count: `128`
- known descendant rank: `119`
- combined rank: `157`

Hence the quotient dimension is
$$
157 - 119 = 38.
$$

In the same `OPEdefs` vector space, the stress-tensor-leading candidate satisfies an explicit linear relation:
$$
\begin{aligned}
T^4
&-\frac{2}{3}T^3(JJ)
+\frac{1}{3}T^3\partial J
-\frac{2}{3}T^2\bigl(J(G\bar G)\bigr)
+\frac{1}{2}T^2(G\partial\bar G) \\
&-\frac{9}{4}T^2(W\partial\bar W)
-\frac{1}{12}T^2(\partial J\,JJ)
-\frac{1}{3}T^2(\partial J\,\partial J)
-\frac{1}{6}T^2(\partial G\,\bar G) \\
&+\frac{15}{4}T^2((\partial W)\bar W)
+\frac{3}{4}T^2(\partial^2J\,J)
=0.
\end{aligned}
$$

Equivalently, in the internal candidate labels used during the computation:

- `NO(SG0,NO(SG0,NO(SG0,SG0)))`
- `- 2/3 NO(SG0,NO(SG0,NO(SG0,NO(SG1,SG1))))`
- `+ 1/3 NO(SG0,NO(SG0,NO(SG0,∂SG1)))`
- `- 2/3 NO(SG0,NO(SG0,NO(SG1,NO(SG2,SG3))))`
- `+ 1/2 NO(SG0,NO(SG0,NO(SG2,∂SG3)))`
- `- 9/4 NO(SG0,NO(SG0,NO(SG4,∂SG5)))`
- `- 1/12 NO(SG0,NO(SG0,NO(∂SG1,NO(SG1,SG1))))`
- `- 1/3 NO(SG0,NO(SG0,NO(∂SG1,∂SG1)))`
- `- 1/6 NO(SG0,NO(SG0,NO(∂SG2,SG3)))`
- `+ 15/4 NO(SG0,NO(SG0,NO(∂SG4,SG5)))`
- `+ 3/4 NO(SG0,NO(SG0,NO(∂^2SG1,SG1)))`

## Update: faithful-cache check

There was an intermediate concern that the cached sparse term representation might fail to preserve the true `OPEdefs.m` expansion. This was checked explicitly.

For the above $11$-term relation, I compared:

1. direct `Expand[sum of all 11 terms]`, and
2. the sum of the cached term tables for the same $11$ individual candidates.

The prefix-sum comparison matched at every stage:

- prefix length `1`: exact match;
- prefix length `2`: exact match;
- ...
- prefix length `10`: exact match;
- prefix length `11`: both sides reduced to `0`.

So the current `OPEdefs` sparse cache is faithful for this relation. The linear dependence is therefore a genuine free-field null relation, not an artifact of the serialization step.

## Revised interpretation

The earlier provisional statement that the weight-$8$ stress-tensor-leading null was unsupported is no longer correct.

The current best-supported conclusion is instead:

- the $p=3$ free-field realization reproduces the low-order nulls and OPE data from section 5.1;
- the low-order Schur / vacuum-character evidence is reproduced;
- there exists a genuine weight-$8$ relation with leading term $T^4$;
- this relation is visible in the pure `OPEdefs.m` realization and survives the faithful-cache check.

What remains open is to rewrite this relation in the cleanest VOA / quasiprimary basis and connect it as directly as possible to the fourth-order modular differential equation motivation.

## Files produced

- ` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_p3.py `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_p3_results.json `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/character_scan_p3.py `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/character_scan_p3_results.json `
- ` .trellis/tasks/03-31-z3-p3-w-algebra/recompute_summary.md `

## Low-order Schur / character check

I also performed an independent low-order vacuum-character scan using
`character_scan_p3.py`.

The scan computes the superdimension at each fixed conformal weight by:

1. enumerating canonical monomials built from a chosen generator set;
2. realizing them in the free-field system;
3. quotienting by linear dependence in the realized free-field basis;
4. taking bosonic dimension minus fermionic dimension.

### Eight-generator scan matches the paper up to $h=4$

Using the generator set
$$
\{T,J,G,\bar G,W,\bar W,GW,G\bar W
\}
$$
with `GW := {GW}_1` and `GWbar := {\bar G\bar W}_1`, the computed unflavored series agrees with the paper through
$$
q^8.
$$

The matched coefficients are
$$
1 + q^2 + q^4 + 2q^6 - 2q^7 + 3q^8.
$$

Equivalently, in VOA weights:

- $h=0$: $1$
- $h=1$: $1$
- $h=2$: $1$
- $h=3$: $2$
- $h=7/2$: $-2$
- $h=4$: $3$

### Six-generator scan does not reproduce the same counting

If one instead performs the same naive monomial counting using only
$$
\{T,J,G,\bar G,W,\bar W\},
$$
the resulting low-order superdimensions are too large and fail already at
$$
h=2.
$$

So for the purpose of direct low-order character enumeration in the current codebase, the descendants
$$
GW, \qquad GWbar
$$
must be included explicitly as atomic generators in order to reproduce the expected Schur data.

This does **not** by itself prove that the minimal strong generating set is eight-dimensional. It only shows that the present counting implementation, when restricted to direct normal-ordered monomial enumeration, does not yet recover the correct low-order cancellations from the six-generator presentation alone.

## Diagnosis of the six-generator mismatch

I analyzed the mismatch more closely in
` .trellis/tasks/03-31-z3-p3-w-algebra/character_gap_analysis_p3.py `
with output in
` .trellis/tasks/03-31-z3-p3-w-algebra/character_gap_analysis_p3_results.json `.

The main conclusion is:

### The problem is missing descendants, not missing low-order quotient relations

At the first mismatching weights, the bosonic counts from the six-generator and eight-generator scans are the same, while the fermionic counts differ:

- at $h=2$: bosons $3$ vs $3$, fermions $0$ vs $2$;
- at $h=3$: bosons $8$ vs $8$, fermions $4$ vs $6$;
- at $h=4$: bosons $19$ vs $19$, fermions $12$ vs $16$.

So the excess in the six-generator superdimension is caused by **missing fermionic states**, not by extra bosonic states surviving an absent null quotient.

### Which states are missing?

At low weights, the missing states are precisely the towers generated by
$$
GW := \{G,W\}_1, \qquad GWbar := \{\bar G,\bar W\}_1,
$$
and their descendants / composites.

For example:

- at $h=2$, the eight-generator scan contains `GW` and `GWbar`, absent from the six-generator monomial list;
- at $h=3$, the new states include
  $$
  \partial GW,\quad \partial GWbar,\quad (GW\,J),\quad (GWbar\,J);
  $$
- at $h=7/2$ and $h=4$, the new states continue as composites such as
  $$
  (GW\,W),\ (GW\,\bar W),\ (GWbar\,W),\ (GWbar\,\bar W),\ (\partial GW\,J),\dots
  $$

### Representative realized relations

The eight-generator analysis also reveals explicit realized identities showing that these states are not arbitrary additions, but the expected supersymmetry-descendant sector. For instance at weight $3$ one finds a relation equivalent to
$$
3\,\partial GW - 3\,(GW) + (GW\,J) = 0,
$$
in the realized free-field basis, i.e. the new sector is tied to the old one by descendant relations.

### Interpretation

This means the current six-generator counting script is not wrong because it forgot to quotient by the Higgs / fermionic nulls. Those nulls were already checked independently and do vanish.

Instead, the issue is that the current direct monomial enumerator does not automatically close the state space under the supersymmetry descendants generated by singular OPEs. In other words, a **strong-generator presentation** is not yet being converted into the full low-weight module in a way that includes states like
$$
GW,\ GWbar
$$
unless they are manually added as atomic generators.

### Practical consequence

For low-order Schur counting in the present codebase, there are two workable options:

1. include `GW` and `GWbar` explicitly in the scanned generator list;
2. upgrade the enumerator so that, starting from the six strong generators, it also generates the necessary singular-OPE descendants automatically.

The current evidence favors option 2 as the conceptually correct long-term fix.
