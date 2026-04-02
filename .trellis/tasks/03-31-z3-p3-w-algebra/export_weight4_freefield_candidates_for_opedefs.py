from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    OPE,
    One,
    clear_registry,
    extract_scalar_operator,
    get_operator_parity,
)
from pyope.operator_spaces import LocalOperatorBasis
from pyope.operators import DerivativeOperator, NormalOrderedOperator


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec_export = importlib.util.spec_from_file_location(
    "expdesc", THIS_DIR / "export_weight8_descendants_for_opedefs.py"
)
expdesc = importlib.util.module_from_spec(spec_export)
assert spec_export.loader is not None
spec_export.loader.exec_module(expdesc)

OUT = THIS_DIR / "weight4_freefield_candidates_for_opedefs.json"


FREE_FIELD_CHARGES = {
    "bff_w4qp": -2,
    "cff_w4qp": 2,
    "betaff_w4qp": 3,
    "gammaff_w4qp": -3,
}


def q(expr):
    if isinstance(expr, BasisOperator):
        return FREE_FIELD_CHARGES.get(expr.name, 0)
    if isinstance(expr, DerivativeOperator):
        return q(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return q(expr.left) + q(expr.right)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return q(op)
    if isinstance(expr, sp.Add):
        vals = {q(arg) for arg in expr.args}
        return next(iter(vals)) if len(vals) == 1 else None
    return 0


def main():
    clear_registry()
    b = BasisOperator("bff_w4qp", fermionic=True, conformal_weight=Fraction(2))
    c = BasisOperator("cff_w4qp", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasisOperator("betaff_w4qp", conformal_weight=Fraction(3, 2))
    gamma = BasisOperator("gammaff_w4qp", conformal_weight=Fraction(-1, 2))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    basis = LocalOperatorBasis(
        [b, c, beta, gamma], max_weight=Fraction(4), max_occurence=12
    )
    ops = [
        op
        for op in basis.list(weight=Fraction(4), max_occurence=12)
        if get_operator_parity(op) == 0 and q(op) == 0
    ]

    payload = {
        "free_fields": {
            "bff": "bff_w4qp",
            "cff": "cff_w4qp",
            "betaff": "betaff_w4qp",
            "gammaff": "gammaff_w4qp",
        },
        "candidates": [
            {
                "expr": str(expr),
                "mma": expdesc.normalize_mma_text(expdesc.to_mma(expr))
                .replace("bff_w4qp", "bff")
                .replace("cff_w4qp", "cff")
                .replace("betaff_w4qp", "betaff")
                .replace("gammaff_w4qp", "gammaff"),
            }
            for expr in ops
        ],
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps({"candidate_count": len(ops)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
