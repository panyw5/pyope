from fractions import Fraction

from pyope import (
    BasicOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    NO,
    OPE,
    One,
    Zero,
    bracket,
    d,
    qp,
    simplify,
)


def build_p3_data():
    b = BasicOperator("b_demo_p3", fermionic=True, conformal_weight=Fraction(2))
    c = BasicOperator("c_demo_p3", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasicOperator("beta_demo_p3", conformal_weight=Fraction(3, 2))
    gamma = BasicOperator("gamma_demo_p3", conformal_weight=Fraction(-1, 2))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = 2 * NO(b, c) + 3 * NO(beta, gamma)
    G = NO(gamma, b)
    Gt = 2 * NO(d(beta), c) + 3 * NO(beta, d(c))
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
    GW = bracket(G, W, 1)
    GWbar = bracket(Gt, Wbar, 1)
    return J, G, Gt, T, W, Wbar, GW, GWbar


def main():
    J, G, Gt, T, W, Wbar, GW, GWbar = build_p3_data()

    print("g_{W Wbar} =", OPE(W, Wbar).pole(3))
    print("(W Wbar)_1 =", qp(W, Wbar, 1))

    JJ0 = qp(J, J, 0)
    print("(J J)_0 =", JJ0)
    print(
        "identity-family level-1 check =",
        simplify(qp(W, Wbar, 1) - (Fraction(2, 3) * T - Fraction(1, 3) * JJ0)),
    )

    X_higgs = simplify(
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gt, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )
    X_minus = simplify(qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2))
    X_plus = simplify(
        qp(Gt, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2)
    )

    print("X_higgs =", X_higgs)
    print("X_minus =", X_minus)
    print("X_plus =", X_plus)
    print("all vanish =", X_higgs == Zero and X_minus == Zero and X_plus == Zero)


if __name__ == "__main__":
    main()
