from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys

from pyope import Zero, get_operator_parity, list_zero_relations, realize_and_simplify
from pyope.operator_spaces import LocalOperatorBasis

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from character_scan_p3 import build_p3_data


def parity_label(op):
    return "bosonic" if get_operator_parity(op) == 0 else "fermionic"


def analyze_weight(weight: Fraction):
    realized_map, free_fields = build_p3_data()
    names6 = ["T", "J", "G", "Gbar", "W", "Wbar"]
    names8 = ["T", "J", "G", "Gbar", "W", "Wbar", "GW", "GWbar"]

    basis6 = LocalOperatorBasis(
        [realized_map[name] for name in names6], max_weight=weight
    )
    basis8 = LocalOperatorBasis(
        [realized_map[name] for name in names8], max_weight=weight
    )
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()), max_weight=weight, max_occurence=12
    )

    ops6 = basis6.list(weight=weight)
    ops8 = basis8.list(weight=weight)

    realized6 = []
    for op in ops6:
        expr = realize_and_simplify(op)
        if expr not in (0, Zero):
            realized6.append((op, expr))

    realized8 = []
    for op in ops8:
        expr = realize_and_simplify(op)
        if expr not in (0, Zero):
            realized8.append((op, expr))

    rels6 = list_zero_relations(
        [expr for _, expr in realized6], ff_basis, weight=weight, max_occurence=12
    )
    rels8 = list_zero_relations(
        [expr for _, expr in realized8], ff_basis, weight=weight, max_occurence=12
    )

    only_in_8 = sorted(
        set(str(op) for op, _ in realized8) - set(str(op) for op, _ in realized6)
    )

    payload = {
        "weight": str(weight),
        "six_candidates": [
            {"op": str(op), "parity": parity_label(op), "realized": str(expr)}
            for op, expr in realized6
        ],
        "eight_candidates": [
            {"op": str(op), "parity": parity_label(op), "realized": str(expr)}
            for op, expr in realized8
        ],
        "only_in_8": only_in_8,
        "six_relation_count": len(rels6),
        "eight_relation_count": len(rels8),
        "six_relations": [
            {
                "terms": [(str(expr), str(coeff)) for expr, coeff in rel["terms"]],
                "coefficients": [str(c) for c in rel["coefficients"]],
            }
            for rel in rels6[:20]
        ],
        "eight_relations": [
            {
                "terms": [(str(expr), str(coeff)) for expr, coeff in rel["terms"]],
                "coefficients": [str(c) for c in rel["coefficients"]],
            }
            for rel in rels8[:20]
        ],
    }
    return payload


def main():
    results = {
        "weight_2": analyze_weight(Fraction(2)),
        "weight_3": analyze_weight(Fraction(3)),
        "weight_7_2": analyze_weight(Fraction(7, 2)),
        "weight_4": analyze_weight(Fraction(4)),
    }
    output_path = Path(__file__).with_name("character_gap_analysis_p3_results.json")
    output_path.write_text(json.dumps(results, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(results, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
