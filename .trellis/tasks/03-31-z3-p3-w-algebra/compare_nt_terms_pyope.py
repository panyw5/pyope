from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import re

from pyope import (
    BasisOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    NO,
    One,
    OPE,
    clear_registry,
    d,
    realize_and_simplify,
    simplify,
)
from pyope.operator_spaces import LocalOperatorBasis


def normalize_pyope_operator(text: str) -> str:
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    return text


def build_data():
    clear_registry()

    b = BasisOperator("b", fermionic=True, conformal_weight=Fraction(2))
    c = BasisOperator("c", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasisOperator("beta", conformal_weight=Fraction(3, 2))
    gamma = BasisOperator("gamma", conformal_weight=Fraction(-1, 2))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = 2 * NO(b, c) + 3 * NO(beta, gamma)
    G = NO(gamma, b)
    Gbar = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))
    T = (
        -2 * NO(b, d(c))
        - Fraction(3, 2) * NO(beta, d(gamma))
        - NO(d(b), c)
        - Fraction(1, 2) * NO(d(beta), gamma)
    )
    W = beta
    GW = b
    Wbar = (
        NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
        + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c))))
        - 4 * NO(beta, NO(d(gamma), gamma))
        - Fraction(4, 3) * NO(gamma, NO(b, d(c)))
        + Fraction(2, 3) * NO(gamma, NO(d(b), c))
        + Fraction(2, 3) * NO(d(beta), NO(gamma, gamma))
        - Fraction(8, 3) * NO(d(gamma), NO(b, c))
        + Fraction(10, 9) * d(d(gamma))
    )
    GWbar = (
        -Fraction(10, 9) * d(d(d(c)))
        - Fraction(8, 3) * NO(b, NO(d(d(c)), c))
        - 3 * NO(beta, NO(beta, NO(gamma, NO(gamma, d(c)))))
        + 4 * NO(beta, NO(gamma, NO(b, NO(d(c), c))))
        + 4 * NO(beta, NO(gamma, d(d(c))))
        + 4 * NO(beta, NO(d(gamma), d(c)))
        - Fraction(2, 3) * NO(d(d(beta)), NO(gamma, c))
        + Fraction(2, 3) * NO(d(b), NO(d(c), c))
        - 2 * NO(d(beta), NO(beta, NO(gamma, NO(gamma, c))))
        + Fraction(8, 3) * NO(d(beta), NO(d(gamma), c))
    )
    return T, J, G, Gbar, W, GW, Wbar, GWbar, [b, c, beta, gamma]


def main():
    T, J, G, Gbar, W, GW, Wbar, GWbar, free_fields = build_data()
    ff_basis = LocalOperatorBasis(free_fields, max_weight=8, max_occurence=12)

    terms = [
        ("term0_T4", NO(T, NO(T, NO(T, T)))),
        ("term1_T3JJ", -Fraction(2, 3) * NO(T, NO(T, NO(T, NO(J, J))))),
        ("term2_T3dJ", Fraction(1, 3) * NO(T, NO(T, NO(T, d(J))))),
        ("term3_T2JGGB", -Fraction(2, 3) * NO(T, NO(T, NO(J, NO(G, Gbar))))),
        ("term4_T2GdGB", Fraction(1, 2) * NO(T, NO(T, NO(G, d(Gbar))))),
        ("term5_T2WdWB", -Fraction(9, 4) * NO(T, NO(T, NO(W, d(Wbar))))),
        ("term6_T2dJJJ", -Fraction(1, 12) * NO(T, NO(T, NO(d(J), NO(J, J))))),
        ("term7_T2dJdJ", -Fraction(1, 3) * NO(T, NO(T, NO(d(J), d(J))))),
        (
            "term8_T2dGGB",
            -Fraction(1, 6) * NO(T, NO(T, NO(d(G), Gbar))),
        ),
        ("term9_T2dWWB", Fraction(15, 4) * NO(T, NO(T, NO(d(W), Wbar)))),
        ("term10_T2d2JJ", Fraction(3, 4) * NO(T, NO(T, NO(d(d(J)), J)))),
    ]

    payload = {}
    for label, expr in terms:
        realized = simplify(realize_and_simplify(expr))
        sparse = ff_basis.sparse_terms(realized)
        serialized = [
            {
                "monomial": normalize_pyope_operator(str(monomial)),
                "coeff": str(coeff),
            }
            for monomial, coeff in sorted(sparse.items(), key=lambda item: str(item[0]))
        ]
        payload[label] = serialized

    output_path = Path(__file__).with_name("compare_nt_terms_pyope.json")
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
