# Examples, Tests, and Reference Map

This file answers a practical question: if you are building a pyope skill, which non-source files should you index first?

## Priority Order For Retrieval Material

### First priority

These files best represent intended user-facing behavior:

- `README.md`
- `src/pyope/__init__.py`
- `tests/algebras/test_virasoro_voa.py`
- `tests/core/test_normal_product.py`
- `tests/research/test_operator_spaces.py`
- `tests/research/test_sparse_c2_api.py`

### Second priority

These files are excellent for broader examples and demonstrations:

- `demo/pyope_basic_demo.ipynb`
- `demo/pyope_ope_demo.ipynb`
- `demo/normal_product_demo.ipynb`
- `demo/jacobi_identity_demo.ipynb`
- `demo/operator_spaces_demo.ipynb`
- `demo/virasoro_c2_algebra_demo.ipynb`
- `demo/w_algebra_null_states_demo.ipynb`

### Third priority

These files are important for historical fidelity and validation against the Mathematica origin:

- `OPEdefs/OPEdefs.m`
- `tests/*mathematica*.wls`
- `demo/mathematica_comparison.ipynb`

## Demo Map By Topic

| Topic | Best demo files |
| --- | --- |
| Core package tour | `demo/pyope_basic_demo.ipynb`, `demo/pyope_ope_demo.ipynb` |
| Normal ordering and simplification | `demo/normal_product_demo.ipynb` |
| Jacobi checks | `demo/jacobi_identity_demo.ipynb` |
| Fixed-weight spaces and realizations | `demo/operator_spaces_demo.ipynb` |
| Virasoro examples | `demo/Virasoro.ipynb`, `demo/virasoro_c2_algebra_demo.ipynb` |
| Null states and W-algebra experiments | `demo/w_algebra_null_states_demo.ipynb`, `demo/N_equals_three_algebra.ipynb` |
| Mathematica comparison | `demo/mathematica_comparison.ipynb` |

## Test Map By Topic

| Topic | Best tests |
| --- | --- |
| Generator declaration and Virasoro OPE | `tests/algebras/test_virasoro_voa.py` |
| OPE registry behavior | `tests/foundation/test_registry.py` |
| Illegal operator products | `tests/foundation/test_illegal_operator_mul.py` |
| Normal products and simplification | `tests/core/test_normal_product.py`, `tests/core/test_simplify.py` |
| Jacobi identity normalization | `tests/core/test_jacobi.py` |
| Quasiprimary helpers | `tests/core/test_quasiprimary.py` |
| Fixed-weight operator spaces | `tests/research/test_operator_spaces.py` |
| Descendant generation | `tests/research/test_descendant_space.py`, `tests/research/test_descendants.py` |
| Sparse `C_2` and quotient search | `tests/research/test_sparse_c2_api.py`, `tests/research/test_c2_null_searcher.py` |
| Singularity and primary constraints | `tests/research/test_singular_vector_analyzer.py`, `tests/research/test_singularity.py` |

## Suggested Retrieval Recipes For The Skill

### Recipe 1: beginner asks for a minimal OPE example

Retrieve:

1. `README.md`
2. `src/pyope/__init__.py`
3. `tests/algebras/test_virasoro_voa.py`

### Recipe 2: user asks about normal ordering or why multiplication fails

Retrieve:

1. `src/pyope/api.py`
2. `src/pyope/simplify.py`
3. `tests/core/test_normal_product.py`
4. `tests/foundation/test_illegal_operator_mul.py`

### Recipe 3: user asks about fixed-weight basis enumeration

Retrieve:

1. `src/pyope/operator_spaces.py`
2. `tests/research/test_operator_spaces.py`
3. `demo/operator_spaces_demo.ipynb`

### Recipe 4: user asks about null states or `C_2`

Retrieve:

1. `src/pyope/c2.py`
2. `src/pyope/null_search.py`
3. `src/pyope/singularity.py`
4. `tests/research/test_sparse_c2_api.py`
5. `tests/research/test_singular_vector_analyzer.py`

## What Not To Use As Primary Skill Context

Avoid centering the skill on these materials unless the user explicitly asks for repository archaeology:

- top-level `tmp_*.py` and `tmp_*.wls`
- experimental one-off verification scripts in the repository root
- notebook outputs as authoritative behavior when an equivalent test exists

## Recommended Skill Framing

If you later turn this into a skill, the skill should describe pyope as:

- a SymPy-backed symbolic OPE library
- conceptually close to `OPEdefs.m`
- strongest on operator definitions, OPEs, normal ordering, and algebra checks
- increasingly capable on fixed-weight spaces and null-state workflows, but with evolving research APIs in the `C_2` and realization layers