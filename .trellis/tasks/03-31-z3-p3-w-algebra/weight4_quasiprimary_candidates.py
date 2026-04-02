from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp

from pyope import extract_scalar_operator, get_operator_parity
from pyope.operator_spaces import LocalOperatorBasis
from pyope.operators import BasisOperator, DerivativeOperator, NormalOrderedOperator


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(charscan)

spec_export = importlib.util.spec_from_file_location(
    "expdesc", THIS_DIR / "export_weight8_descendants_for_opedefs.py"
)
expdesc = importlib.util.module_from_spec(spec_export)
assert spec_export.loader is not None
spec_export.loader.exec_module(expdesc)

OUT = THIS_DIR / "weight4_quasiprimary_candidates.json"

FREE_FIELD_CHARGES = {
    "b_char_p3": -2,
    "c_char_p3": 2,
    "beta_char_p3": 3,
    "gamma_char_p3": -3,
}


def ff_charge(expr):
    if isinstance(expr, BasisOperator):
        return FREE_FIELD_CHARGES.get(expr.name, 0)
    if isinstance(expr, DerivativeOperator):
        return ff_charge(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return ff_charge(expr.left) + ff_charge(expr.right)
    if hasattr(expr, "realization"):
        return ff_charge(expr.realization)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return ff_charge(op)
    if isinstance(expr, sp.Add):
        vals = {ff_charge(arg) for arg in expr.args}
        return next(iter(vals)) if len(vals) == 1 else None
    return 0


def monomial_charge(expr, gen_charges):
    if hasattr(expr, "realization") and hasattr(expr, "name"):
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
    return None


def main():
    raw_map, _, free_fields = charscan.build_p3_data()
    raw_generators, _ = charscan.singular_ope_closure(
        [raw_map[name] for name in ["T", "J", "G", "Gbar", "W", "Wbar"]],
        free_fields,
        Fraction(4),
    )
    generators = charscan.promote_named_generators(raw_generators)
    basis = LocalOperatorBasis(generators, max_weight=Fraction(4))
    gen_charges = {g.name: ff_charge(g.realization) for g in generators}
    ops = [
        op
        for op in basis.list(weight=Fraction(4))
        if get_operator_parity(op) == 0 and monomial_charge(op, gen_charges) == 0
    ]
    payload = {
        "generator_defs": {
            f"SG{i}": expdesc.normalize_mma_text(expdesc.to_mma(expr))
            for i, expr in enumerate(raw_generators)
        },
        "candidates": [
            {"expr": str(expr), "mma": expdesc.normalize_mma_text(expdesc.to_mma(expr))}
            for expr in ops
        ],
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps({"candidate_count": len(ops)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
