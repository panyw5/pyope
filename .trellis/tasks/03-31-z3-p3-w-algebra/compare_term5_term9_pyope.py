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


def norm(text: str) -> str:
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    return text


def build():
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
    return T, W, Wbar, [b, c, beta, gamma]


def sparse(expr, ff_basis):
    realized = simplify(realize_and_simplify(expr))
    return [
        {"monomial": norm(str(m)), "coeff": str(c)}
        for m, c in sorted(
            ff_basis.sparse_terms(realized).items(), key=lambda item: str(item[0])
        )
    ]


def main():
    T, W, Wbar, free = build()
    ff_basis = LocalOperatorBasis(free, max_weight=8, max_occurence=12)
    data = {
        "Wbar": sparse(Wbar, ff_basis),
        "dWbar": sparse(d(Wbar), ff_basis),
        "NO_W_dWbar": sparse(NO(W, d(Wbar)), ff_basis),
        "NO_dW_Wbar": sparse(NO(d(W), Wbar), ff_basis),
        "term5": sparse(-Fraction(9, 4) * NO(T, NO(T, NO(W, d(Wbar)))), ff_basis),
        "term9": sparse(Fraction(15, 4) * NO(T, NO(T, NO(d(W), Wbar))), ff_basis),
        "term5_plus_term9": sparse(
            -Fraction(9, 4) * NO(T, NO(T, NO(W, d(Wbar))))
            + Fraction(15, 4) * NO(T, NO(T, NO(d(W), Wbar))),
            ff_basis,
        ),
    }
    out = Path(__file__).with_name("compare_term5_term9_pyope.json")
    out.write_text(json.dumps(data, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps({k: len(v) for k, v in data.items()}, indent=2))


if __name__ == "__main__":
    main()
