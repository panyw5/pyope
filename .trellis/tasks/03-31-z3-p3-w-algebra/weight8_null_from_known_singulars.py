from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp

from pyope import (
    Zero,
    extract_scalar_operator,
    get_operator_parity,
    normal_product,
    qp,
    realize_and_simplify,
)
from pyope.operator_spaces import (
    LocalOperatorBasis,
    RealizedGenerator,
    _get_conformal_weight,
)
from pyope.operators import BasisOperator, DerivativeOperator, NormalOrderedOperator, d


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(charscan)


FREE_FIELD_CHARGES = {
    "b_char_p3": -2,
    "c_char_p3": 2,
    "beta_char_p3": 3,
    "gamma_char_p3": -3,
}


def ff_charge(expr):
    if isinstance(expr, BasisOperator) and not isinstance(expr, RealizedGenerator):
        return FREE_FIELD_CHARGES.get(expr.name, 0)
    if isinstance(expr, DerivativeOperator):
        return ff_charge(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return ff_charge(expr.left) + ff_charge(expr.right)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return ff_charge(op)
    if isinstance(expr, sp.Add):
        charges = {ff_charge(arg) for arg in expr.args}
        return next(iter(charges)) if len(charges) == 1 else None
    return 0


def monomial_charge(expr, gen_charges):
    if isinstance(expr, RealizedGenerator):
        return gen_charges.get(expr.name)
    if isinstance(expr, DerivativeOperator):
        return monomial_charge(expr.base, gen_charges)
    if isinstance(expr, NormalOrderedOperator):
        q_left = monomial_charge(expr.left, gen_charges)
        q_right = monomial_charge(expr.right, gen_charges)
        return None if q_left is None or q_right is None else q_left + q_right
    if isinstance(expr, sp.Mul):
        _, inner = extract_scalar_operator(expr)
        return monomial_charge(inner, gen_charges)
    return None


def is_c2_monomial(expr):
    if isinstance(expr, DerivativeOperator):
        return True
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return is_c2_monomial(op)
    if isinstance(expr, NormalOrderedOperator):
        return isinstance(expr.left, DerivativeOperator)
    return False


def sparse_rank(vectors):
    if not vectors:
        return 0
    monomials = sorted({m for vec in vectors for m in vec}, key=sp.srepr)
    if not monomials:
        return 0
    matrix = sp.Matrix(
        [
            [sp.sympify(vec.get(monomial, 0)) for vec in vectors]
            for monomial in monomials
        ]
    )
    return int(matrix.rank())


def ideal_descendants_at_weight(source, target_weight, basis, gen_charges):
    source_weight = _get_conformal_weight(source)
    source_parity = get_operator_parity(source)
    source_charge = monomial_charge(source, gen_charges)
    if source_weight is None or source_charge is None:
        return []

    descendants = set()
    max_deriv = int(target_weight - source_weight)
    for k in range(max_deriv + 1):
        phi_weight = target_weight - source_weight - k
        phi_basis = basis.list(phi_weight)
        for phi in phi_basis:
            if get_operator_parity(phi) != source_parity:
                continue
            if monomial_charge(phi, gen_charges) != -source_charge:
                continue
            expr = basis.canonicalize(d(normal_product(phi, source), k))
            if expr != Zero and _get_conformal_weight(expr) == target_weight:
                descendants.add(expr)
    return sorted(descendants, key=sp.srepr)


def main():
    weight = Fraction(8)
    raw_map, _, free_fields = charscan.build_p3_data()
    closed_raw, closure = charscan.singular_ope_closure(
        [raw_map[name] for name in ["T", "J", "G", "Gbar", "W", "Wbar"]],
        free_fields,
        weight,
    )
    formal_generators = [
        RealizedGenerator(f"FG{i}", realization=expr)
        for i, expr in enumerate(closed_raw)
    ]
    T, J, G, Gbar, W, Wbar, GW, GWbar = formal_generators

    formal_basis = LocalOperatorBasis(formal_generators, max_weight=weight)
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()), max_weight=weight, max_occurence=12
    )
    gen_charges = {g.name: ff_charge(g.realization) for g in formal_generators}

    candidates = [
        op
        for op in formal_basis.list(weight=weight)
        if get_operator_parity(op) == 0
        and monomial_charge(op, gen_charges) == 0
        and not is_c2_monomial(op)
    ]

    JJ0 = qp(J, J, 0)
    X_higgs = (
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gbar, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )
    X_minus = qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2)
    X_plus = qp(Gbar, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2)

    source_data = {}
    source_vectors = []
    for label, source in [
        ("X_higgs", X_higgs),
        ("X_minus", X_minus),
        ("X_plus", X_plus),
    ]:
        closure_desc = ideal_descendants_at_weight(
            source, weight, formal_basis, gen_charges
        )
        realized_vectors = []
        for expr in closure_desc:
            vec = ff_basis.sparse_terms(realize_and_simplify(expr))
            if vec:
                realized_vectors.append(vec)
        source_vectors.extend(realized_vectors)
        source_data[label] = {
            "source": str(source),
            "descendant_count": len(closure_desc),
            "realized_rank": sparse_rank(realized_vectors),
            "sample_descendants": [str(expr) for expr in closure_desc[:20]],
        }

    source_rank = sparse_rank(source_vectors)
    candidate_vectors = [
        ff_basis.sparse_terms(realize_and_simplify(expr)) for expr in candidates
    ]
    candidate_rank = sparse_rank(candidate_vectors)
    combined_rank = sparse_rank(candidate_vectors + source_vectors)

    target = next(
        expr for expr in candidates if str(expr) == "NO(FG0,NO(FG0,NO(FG0,FG0)))"
    )
    rank_without_target = sparse_rank(
        source_vectors
        + [vec for expr, vec in zip(candidates, candidate_vectors) if expr != target]
    )
    rank_with_target = sparse_rank(source_vectors + candidate_vectors)

    payload = {
        "closed_generator_count": len(formal_generators),
        "closure_discoveries": closure,
        "candidate_count": len(candidates),
        "candidate_rank": candidate_rank,
        "known_singular_descendants": source_data,
        "known_singular_descendant_rank": source_rank,
        "combined_rank": combined_rank,
        "quotient_dimension_after_known_singulars": combined_rank - source_rank,
        "target": str(target),
        "target_survives_mod_known_singulars": rank_with_target > rank_without_target,
        "rank_without_target": rank_without_target,
        "rank_with_target": rank_with_target,
    }

    output_path = THIS_DIR / "weight8_null_from_known_singulars.json"
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
