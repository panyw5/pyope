from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import re
import sys

import sympy as sp

from pyope import extract_scalar_operator, get_operator_parity, normal_product, qp
from pyope.operator_spaces import (
    LocalOperatorBasis,
    RealizedGenerator,
    _get_conformal_weight,
)
from pyope.operators import DerivativeOperator, NormalOrderedOperator, d


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(charscan)

EXPORT_PATH = THIS_DIR / "weight8_descendants_for_opedefs.json"
WEIGHT = Fraction(8)


def to_mma(expr) -> str:
    if isinstance(expr, sp.Integer):
        return str(int(expr))
    if isinstance(expr, sp.Rational):
        return f"{expr.p}/{expr.q}" if expr.q != 1 else str(expr.p)
    base = getattr(expr, "base", None)
    order = getattr(expr, "order", None)
    left = getattr(expr, "left", None)
    right = getattr(expr, "right", None)
    if base is not None and order is not None:
        return f"Derivative[{order}][{to_mma(base)}]"
    if left is not None and right is not None:
        return f"NO[{to_mma(left)}, {to_mma(right)}]"
    if isinstance(expr, sp.Symbol):
        return str(expr)
    if isinstance(expr, RealizedGenerator):
        return expr.name
    if isinstance(expr, sp.Add):
        parts = []
        for arg in expr.args:
            text = to_mma(arg)
            if text.startswith("-"):
                parts.append(text)
            else:
                parts.append(f"+{text}")
        joined = "".join(parts)
        return joined[1:] if joined.startswith("+") else joined
    if isinstance(expr, sp.Mul):
        coeff, op = extract_scalar_operator(expr)
        coeff = sp.sympify(coeff)
        op_text = to_mma(op)
        coeff_text = to_mma(coeff)
        if op == 1:
            return coeff_text
        if coeff == 1:
            return op_text
        if coeff == -1:
            return f"-{op_text}"
        return f"{coeff_text}*{op_text}"
    raise TypeError(
        f"Unsupported expression for Mathematica export: {type(expr)} :: {expr}"
    )


def normalize_mma_text(text: str) -> str:
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    text = text.replace("b_char_p3", "bff")
    text = text.replace("c_char_p3", "cff")
    text = text.replace("beta_char_p3", "betaff")
    text = text.replace("gamma_char_p3", "gammaff")
    return text


def ff_charge(expr):
    if isinstance(expr, RealizedGenerator):
        return ff_charge(expr.realization)
    name = getattr(expr, "name", None)
    if name == "b_char_p3":
        return -2
    if name == "c_char_p3":
        return 2
    if name == "beta_char_p3":
        return 3
    if name == "gamma_char_p3":
        return -3
    if isinstance(expr, DerivativeOperator):
        return ff_charge(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return ff_charge(expr.left) + ff_charge(expr.right)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return ff_charge(op)
    if isinstance(expr, sp.Add):
        vals = {ff_charge(arg) for arg in expr.args}
        return next(iter(vals)) if len(vals) == 1 else None
    return 0


def monomial_charge(expr, gen_charges):
    if isinstance(expr, RealizedGenerator):
        return gen_charges.get(expr.name)
    if isinstance(expr, DerivativeOperator):
        return monomial_charge(expr.base, gen_charges)
    if isinstance(expr, NormalOrderedOperator):
        ql = monomial_charge(expr.left, gen_charges)
        qr = monomial_charge(expr.right, gen_charges)
        return None if ql is None or qr is None else ql + qr
    if isinstance(expr, sp.Mul):
        _, inner = extract_scalar_operator(expr)
        return monomial_charge(inner, gen_charges)
    if isinstance(expr, sp.Add):
        vals = {monomial_charge(arg, gen_charges) for arg in expr.args}
        return next(iter(vals)) if len(vals) == 1 else None
    return None


def main():
    raw_map, _, free_fields = charscan.build_p3_data()
    raw_generators, closure = charscan.singular_ope_closure(
        [raw_map[name] for name in ["T", "J", "G", "Gbar", "W", "Wbar"]],
        free_fields,
        WEIGHT,
    )
    formal_generators = [
        RealizedGenerator(f"FG{i}", realization=expr)
        for i, expr in enumerate(raw_generators)
    ]
    formal_basis = LocalOperatorBasis(formal_generators, max_weight=WEIGHT)
    gen_charges = {g.name: ff_charge(g.realization) for g in formal_generators}

    T, J, G, Gbar, W, Wbar, GW, GWbar = formal_generators
    JJ0 = qp(J, J, 0)
    sources = {
        "X_higgs": qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gbar, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2),
        "X_minus": qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2),
        "X_plus": qp(Gbar, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2),
    }

    descendants = {}
    all_exprs = []
    for label, source in sources.items():
        source_weight = _get_conformal_weight(source)
        source_parity = get_operator_parity(source)
        source_charge = monomial_charge(source, gen_charges)
        desc = set()
        for k in range(int(WEIGHT - source_weight) + 1):
            phi_weight = WEIGHT - source_weight - k
            for phi in formal_basis.list(phi_weight):
                if get_operator_parity(phi) != source_parity:
                    continue
                if monomial_charge(phi, gen_charges) != -source_charge:
                    continue
                expr = formal_basis.canonicalize(d(normal_product(phi, source), k))
                if expr != 0 and _get_conformal_weight(expr) == WEIGHT:
                    desc.add(expr)
        ordered = sorted(desc, key=sp.srepr)
        descendants[label] = [str(expr) for expr in ordered]
        all_exprs.extend(ordered)

    unique_exprs = sorted(set(all_exprs), key=sp.srepr)
    payload = {
        "generator_defs": {
            f"FG{i}": normalize_mma_text(to_mma(expr))
            for i, expr in enumerate(raw_generators)
        },
        "sources": {
            name: normalize_mma_text(to_mma(expr)) for name, expr in sources.items()
        },
        "descendants_by_source": descendants,
        "all_descendants": [
            {"expr": str(expr), "mma": normalize_mma_text(to_mma(expr))}
            for expr in unique_exprs
        ],
    }
    EXPORT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps({"descendant_count": len(unique_exprs)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
