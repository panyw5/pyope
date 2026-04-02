from fractions import Fraction

from pyope import (
    BasisOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    NO,
    OPE,
    One,
    bracket,
    compact_family_poles,
    d,
    qp,
    simplify,
)


def build_p4_data():
    p = 4
    b = BasisOperator("b_p4_cmp", fermionic=True, conformal_weight=Fraction(p + 1, 2))
    c = BasisOperator("c_p4_cmp", fermionic=True, conformal_weight=Fraction(1 - p, 2))
    beta = BasisOperator("beta_p4_cmp", conformal_weight=Fraction(p, 2))
    gamma = BasisOperator("gamma_p4_cmp", conformal_weight=Fraction(2 - p, 2))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = Fraction(p, 2) * NO(beta, gamma) + Fraction(p - 1, 2) * NO(b, c)
    G = NO(b, gamma)
    Gt = p * NO(beta, d(c)) + (p - 1) * NO(d(beta), c)
    T = (
        -Fraction(p, 2) * NO(beta, d(gamma))
        + (1 - Fraction(p, 2)) * NO(d(beta), gamma)
        - Fraction(p + 1, 2) * NO(b, d(c))
        + Fraction(1 - p, 2) * NO(d(b), c)
    )
    W = beta
    coeffs = [
        Fraction(-35, 64),
        Fraction(15, 16),
        Fraction(-15, 8),
        Fraction(3, 16),
        Fraction(75, 16),
        Fraction(-45, 16),
        Fraction(45, 16),
        Fraction(15, 4),
        Fraction(75, 16),
        Fraction(-75, 16),
        Fraction(3, 16),
        Fraction(-15, 4),
        Fraction(9, 4),
        Fraction(-45, 4),
        Fraction(9, 4),
        Fraction(-15, 2),
        Fraction(9, 4),
        Fraction(3, 1),
        Fraction(1, 1),
    ]
    monomials = [
        d(gamma, 3),
        NO(gamma, NO(b, d(c, 2))),
        NO(gamma, NO(d(b), d(c))),
        NO(gamma, NO(d(b, 2), c)),
        NO(d(gamma), NO(b, d(c))),
        NO(d(gamma), NO(d(b), c)),
        NO(d(gamma, 2), NO(b, c)),
        NO(beta, NO(d(gamma, 2), gamma)),
        NO(beta, NO(d(gamma), d(gamma))),
        NO(d(beta), NO(d(gamma), gamma)),
        NO(d(beta, 2), NO(gamma, gamma)),
        NO(beta, NO(gamma, NO(gamma, NO(b, d(c))))),
        NO(beta, NO(gamma, NO(gamma, NO(d(b), c)))),
        NO(beta, NO(d(gamma), NO(gamma, NO(b, c)))),
        NO(d(beta), NO(gamma, NO(gamma, NO(b, c)))),
        NO(beta, NO(beta, NO(d(gamma), NO(gamma, gamma)))),
        NO(d(beta), NO(beta, NO(gamma, NO(gamma, gamma)))),
        NO(beta, NO(beta, NO(gamma, NO(gamma, NO(gamma, NO(b, c)))))),
        NO(beta, NO(beta, NO(beta, NO(gamma, NO(gamma, NO(gamma, gamma)))))),
    ]
    Wbar = sum(coeff * monomial for coeff, monomial in zip(coeffs, monomials))
    return J, G, Gt, T, W, Wbar


def main():
    J, G, Gt, T, W, Wbar = build_p4_data()
    actual = {pole: simplify(bracket(W, Wbar, pole)) for pole in [4, 3, 2, 1]}

    g = Fraction(105, 32)
    JJ0 = qp(J, J, 0)
    JJJ0 = qp(J, JJ0, 0, right_weight=2)
    GG0 = qp(G, Gt, 0)
    JT0 = qp(J, T, 0, right_weight=2)

    compact_terms = [
        (g, One, 0),
        (g * Fraction(-4, 7), J, 1),
        (g * Fraction(-2, 7), T, 2),
        (g * Fraction(1, 7), JJ0, 2),
        (g * Fraction(-2, 105), JJJ0, 3),
        (g * Fraction(2, 35), GG0, 3),
        (g * Fraction(6, 35), JT0, 3),
    ]
    expected = compact_family_poles(Fraction(2), Fraction(2), compact_terms)

    print("p = 4 compact-family comparison")
    for pole in [4, 3, 2, 1]:
        print(f"pole {pole} actual   = {actual[pole]}")
        print(f"pole {pole} expected = {expected[pole]}")
        print(f"pole {pole} diff     = {simplify(actual[pole] - expected[pole])}")
        print()


if __name__ == "__main__":
    main()
