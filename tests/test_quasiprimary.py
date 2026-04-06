from fractions import Fraction

from pyope import (
    BasisOperator,
    Bosonic,
    compact_family_poles,
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


def _z3_rank_one_data():
    b = BasisOperator("b_z3_qp", fermionic=True, conformal_weight=Fraction(2))
    c = BasisOperator("c_z3_qp", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasisOperator("beta_z3_qp", conformal_weight=Fraction(3, 2))
    gamma = BasisOperator("gamma_z3_qp", conformal_weight=Fraction(-1, 2))

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

    return {
        "b": b,
        "c": c,
        "beta": beta,
        "gamma": gamma,
        "J": J,
        "G": G,
        "Gt": Gt,
        "T": T,
        "W": W,
        "Wbar": Wbar,
        "GW": bracket(G, W, 1),
        "GWbar": bracket(Gt, Wbar, 1),
    }


def _z2_rank_one_data():
    b = BasisOperator("b_z2_qp", fermionic=True, conformal_weight=Fraction(3, 2))
    c = BasisOperator("c_z2_qp", fermionic=True, conformal_weight=Fraction(-1, 2))
    beta = BasisOperator("beta_z2_qp", conformal_weight=Fraction(1))
    gamma = BasisOperator("gamma_z2_qp", conformal_weight=Fraction(0))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = NO(b, c) + 2 * NO(beta, gamma)
    T = (
        -Fraction(3, 2) * NO(b, d(c))
        - NO(beta, d(gamma))
        - Fraction(1, 2) * NO(d(b), c)
    )
    W = beta
    Wbar = -NO(beta, NO(gamma, gamma)) - NO(gamma, NO(b, c)) + Fraction(3, 2) * d(gamma)

    return {
        "J": J,
        "T": T,
        "W": W,
        "Wbar": Wbar,
    }


def _z4_rank_one_data():
    p = 4
    b = BasisOperator("b_z4_qp", fermionic=True, conformal_weight=Fraction(p + 1, 2))
    c = BasisOperator("c_z4_qp", fermionic=True, conformal_weight=Fraction(1 - p, 2))
    beta = BasisOperator("beta_z4_qp", conformal_weight=Fraction(p, 2))
    gamma = BasisOperator("gamma_z4_qp", conformal_weight=Fraction(2 - p, 2))

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = Fraction(p, 2) * NO(beta, gamma) + Fraction(p - 1, 2) * NO(b, c)
    G = NO(b, gamma)
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

    return {
        "J": J,
        "G": G,
        "T": T,
        "W": W,
        "Wbar": Wbar,
    }


def test_quasiprimary_completion_matches_paper_w_wbar_identity_family():
    data = _z3_rank_one_data()
    J = data["J"]
    T = data["T"]
    W = data["W"]
    Wbar = data["Wbar"]

    assert OPE(W, Wbar).pole(3) == -Fraction(20, 9) * One
    assert simplify(bracket(W, Wbar, 2) - Fraction(4, 3) * J) == Zero

    JJ0 = qp(J, J, 0)
    expected_first_quasiprimary = simplify(Fraction(2, 3) * T - Fraction(1, 3) * JJ0)
    assert simplify(qp(W, Wbar, 1) - expected_first_quasiprimary) == Zero


def test_compact_family_reconstruction_matches_p3_identity_family():
    data = _z3_rank_one_data()
    J = data["J"]
    T = data["T"]
    W = data["W"]
    Wbar = data["Wbar"]

    g = Fraction(-20, 9)
    JJ0 = qp(J, J, 0)
    compact_terms = [
        (g, One, 0),
        (g * Fraction(-3, 5), J, 1),
        (g * Fraction(-3, 10), T, 2),
        (g * Fraction(3, 20), JJ0, 2),
    ]
    reconstructed = compact_family_poles(Fraction(3, 2), Fraction(3, 2), compact_terms)

    for pole in [3, 2, 1]:
        assert simplify(bracket(W, Wbar, pole) - reconstructed[pole]) == Zero


def test_p2_w_wbar_matches_identity_family_from_paper():
    data = _z2_rank_one_data()
    J = data["J"]
    W = data["W"]
    Wbar = data["Wbar"]

    assert OPE(W, Wbar).pole(2) == -Fraction(3, 2) * One
    assert simplify(bracket(W, Wbar, 1) - J) == Zero


def test_p4_wbar_matches_free_field_constraints_and_leading_ope_data():
    data = _z4_rank_one_data()
    J = data["J"]
    G = data["G"]
    W = data["W"]
    Wbar = data["Wbar"]

    assert simplify(bracket(G, Wbar, 3)) == 0
    assert simplify(bracket(G, Wbar, 2)) == 0
    assert simplify(bracket(G, Wbar, 1)) == 0

    for pole in [5, 4, 3, 2, 1]:
        assert simplify(bracket(Wbar, Wbar, pole)) == 0

    assert OPE(W, Wbar).pole(4) == Fraction(105, 32) * One
    assert simplify(bracket(W, Wbar, 3) + Fraction(15, 4) * J) == Zero


def test_quasiprimary_completion_reproduces_p3_higgs_null():
    data = _z3_rank_one_data()
    J = data["J"]
    G = data["G"]
    Gt = data["Gt"]
    T = data["T"]
    W = data["W"]
    Wbar = data["Wbar"]

    JJ0 = qp(J, J, 0)
    X = simplify(
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gt, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )

    assert X == Zero


def test_quasiprimary_completion_reproduces_p3_fermionic_nulls():
    data = _z3_rank_one_data()
    J = data["J"]
    G = data["G"]
    Gt = data["Gt"]
    W = data["W"]
    Wbar = data["Wbar"]
    GW = data["GW"]
    GWbar = data["GWbar"]

    Xminus = simplify(qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2))
    Xplus = simplify(qp(Gt, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2))

    assert Xminus == Zero
    assert Xplus == Zero


def test_p3_higgs_null_ope_generates_fermionic_null_ideal():
    data = _z3_rank_one_data()
    J = data["J"]
    G = data["G"]
    Gt = data["Gt"]
    T = data["T"]
    W = data["W"]
    Wbar = data["Wbar"]
    GW = data["GW"]
    GWbar = data["GWbar"]

    JJ0 = qp(J, J, 0)
    Xhiggs = simplify(
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gt, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )
    Xminus = simplify(qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2))
    Xplus = simplify(qp(Gt, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2))

    assert simplify(bracket(GW, Xhiggs, 2)) == Zero
    assert simplify(bracket(GWbar, Xhiggs, 2)) == Zero
    assert Xminus == Zero
    assert Xplus == Zero
