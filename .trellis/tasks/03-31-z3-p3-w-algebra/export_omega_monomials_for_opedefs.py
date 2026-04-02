from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp

from pyope import NO, d, simplify, extract_scalar_operator


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec_charscan = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec_charscan)
assert spec_charscan.loader is not None
spec_charscan.loader.exec_module(charscan)

spec_export = importlib.util.spec_from_file_location(
    "expdesc", THIS_DIR / "export_weight8_descendants_for_opedefs.py"
)
expdesc = importlib.util.module_from_spec(spec_export)
assert spec_export.loader is not None
spec_export.loader.exec_module(expdesc)

OUT = THIS_DIR / "omega_weight4_monomials_opedefs.json"


def split_terms(expr):
    expanded = sp.expand(expr)
    if isinstance(expanded, sp.Add):
        return list(expanded.args)
    return [expanded]


def main():
    raw_map, _, _ = charscan.build_p3_data()
    W = raw_map["W"]
    Wbar = raw_map["Wbar"]
    omega = simplify(NO(d(W), Wbar) - NO(W, d(Wbar)))

    terms = []
    for term in split_terms(omega):
        coeff, op = extract_scalar_operator(term)
        terms.append(
            {
                "coeff": str(sp.sympify(coeff)),
                "expr": str(op),
                "mma": expdesc.normalize_mma_text(expdesc.to_mma(op)),
            }
        )

    payload = {
        "generator_defs": json.loads(
            (THIS_DIR / "weight8_descendants_for_opedefs.json").read_text()
        )["generator_defs"],
        "terms": terms,
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps({"term_count": len(terms)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
