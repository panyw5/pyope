from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    DerivativeOperator,
    extract_scalar_operator,
    GenericC2Reducer,
    Fermionic,
    MakeOPE,
    NO,
    NO_product,
    NormalOrderedOperator,
    OPE,
    One,
    Zero,
    bracket,
    check_jacobi_identity,
    clear_registry,
    d,
    make_realized,
    qp,
    QuotientC2NullSearcher,
    realize_and_simplify,
    simplify,
)
from pyope.descendants import DescendantSpace
from pyope.operator_spaces import LocalOperatorBasis, LocalOperatorCanonicalizer


def build_p3_data():
    clear_registry()

    b = BasisOperator("b_recompute_p3", fermionic=True, conformal_weight=Fraction(2))
    c = BasisOperator("c_recompute_p3", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasisOperator("beta_recompute_p3", conformal_weight=Fraction(3, 2))
    gamma = BasisOperator("gamma_recompute_p3", conformal_weight=Fraction(-1, 2))

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
    GW = bracket(G, W, 1)
    GWbar = bracket(Gbar, Wbar, 1)

    return {
        "raw": {
            "J": J,
            "G": G,
            "Gbar": Gbar,
            "T": T,
            "W": W,
            "Wbar": Wbar,
            "GW": GW,
            "GWbar": GWbar,
        },
        "realized": {
            name: value
            for name, value in zip(
                ["T", "W", "J", "G", "Gbar", "GW", "Wbar", "GWbar"],
                make_realized([T, W, J, G, Gbar, GW, Wbar, GWbar]),
            )
        },
    }


def paper_checks(data):
    raw = data["raw"]
    J = raw["J"]
    G = raw["G"]
    Gbar = raw["Gbar"]
    T = raw["T"]
    W = raw["W"]
    Wbar = raw["Wbar"]
    GW = raw["GW"]
    GWbar = raw["GWbar"]

    JJ0 = qp(J, J, 0)
    checks = {}
    checks["g_WWbar_matches_paper"] = (
        simplify(OPE(W, Wbar).pole(3) + Fraction(20, 9) * One) == Zero
    )
    checks["J_term_matches_paper"] = (
        simplify(bracket(W, Wbar, 2) - Fraction(4, 3) * J) == Zero
    )
    checks["T_JJ_term_matches_paper"] = (
        simplify(qp(W, Wbar, 1) - (Fraction(2, 3) * T - Fraction(1, 3) * JJ0)) == Zero
    )

    X_higgs = simplify(
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gbar, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )
    X_minus = simplify(qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2))
    X_plus = simplify(
        qp(Gbar, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2)
    )

    checks["X_higgs_vanishes"] = X_higgs == Zero
    checks["X_minus_vanishes"] = X_minus == Zero
    checks["X_plus_vanishes"] = X_plus == Zero
    checks["GW_on_higgs_matches_paper"] = (
        simplify(bracket(GW, X_higgs, 2) - Fraction(2, 3) * X_minus) == Zero
    )
    checks["GWbar_on_higgs_matches_paper"] = (
        simplify(bracket(GWbar, X_higgs, 2) - Fraction(4, 3) * X_plus) == Zero
    )
    checks["G_annihilates_Wbar"] = all(
        simplify(bracket(G, Wbar, pole)) in (Zero, 0) for pole in (1, 2, 3)
    )

    payload = {
        "checks": checks,
        "X_higgs": str(X_higgs),
        "X_minus": str(X_minus),
        "X_plus": str(X_plus),
    }
    return payload


def jacobi_checks(data):
    gens = {
        "T": data["raw"]["T"],
        "J": data["raw"]["J"],
        "G": data["raw"]["G"],
        "Gbar": data["raw"]["Gbar"],
        "W": data["raw"]["W"],
        "Wbar": data["raw"]["Wbar"],
    }
    triples = [
        ("T", "T", "T"),
        ("T", "J", "J"),
        ("T", "G", "Gbar"),
        ("J", "G", "Gbar"),
        ("T", "W", "Wbar"),
        ("J", "W", "Wbar"),
        ("G", "W", "Wbar"),
        ("Gbar", "W", "Wbar"),
        ("W", "W", "Wbar"),
        ("W", "Wbar", "Wbar"),
    ]
    results = {}
    for names in triples:
        matrix = check_jacobi_identity(
            *(gens[name] for name in names), simplify_func=simplify
        )
        nonzero_entries = []
        for m_index, row in enumerate(matrix, start=1):
            for n_index, value in enumerate(row, start=1):
                if value != Zero:
                    nonzero_entries.append(
                        {
                            "m": m_index,
                            "n": n_index,
                            "value": str(value),
                        }
                    )
        results["-".join(names)] = {
            "ok": len(nonzero_entries) == 0,
            "nonzero_entries": nonzero_entries,
        }
    return results


def stress_tensor_c2_precheck(data):
    realized = data["realized"]
    generators = [
        realized["T"],
        realized["W"],
        realized["J"],
        realized["G"],
        realized["Gbar"],
        realized["GW"],
        realized["Wbar"],
        realized["GWbar"],
    ]
    basis = LocalOperatorBasis(generators, max_weight=8)
    canonicalizer = LocalOperatorCanonicalizer(
        generators, stress_tensor=realized["T"], max_weight=8
    )
    descendants = DescendantSpace(canonicalizer)
    reducer = GenericC2Reducer(canonicalizer)
    searcher = QuotientC2NullSearcher(
        canonicalizer=canonicalizer,
        descendants=descendants,
        c2_reducer=reducer,
    )
    target = canonicalizer.nested_stress_tensor(4)
    precheck = searcher.quotient_precheck(target)
    result = {
        "status": precheck.status,
        "target_weight": str(precheck.target_weight),
        "quotient_remainder": str(precheck.quotient_remainder),
        "c2_part": str(precheck.c2_remainder),
        "num_c2_generators": len(precheck.witness.generators),
        "basis_weight_8_dim": len(basis.list(weight=8)),
    }
    return result


def free_field_c2_exclusion(data):
    raw = data["raw"]
    target = NO_product(raw["T"], raw["T"], raw["T"], raw["T"])
    target_real = realize_and_simplify(target)

    def is_ff_c2_monomial(op):
        if isinstance(op, DerivativeOperator):
            return True
        if isinstance(op, NormalOrderedOperator):
            return isinstance(op.left, DerivativeOperator)
        return False

    if target_real in (0, Zero):
        return {
            "target_realizes_to_zero": True,
            "all_terms_in_ff_c2": True,
            "non_c2_term_count": 0,
            "sample_non_c2_terms": [],
        }

    terms = list(target_real.args) if isinstance(target_real, sp.Add) else [target_real]
    non_c2_terms = []
    for term in terms:
        op = extract_scalar_operator(term)[1] if isinstance(term, sp.Mul) else term
        if not is_ff_c2_monomial(op):
            non_c2_terms.append(str(term))

    return {
        "target_realizes_to_zero": False,
        "all_terms_in_ff_c2": len(non_c2_terms) == 0,
        "non_c2_term_count": len(non_c2_terms),
        "sample_non_c2_terms": non_c2_terms[:8],
    }


def main():
    data = build_p3_data()
    summary = {
        "paper_checks": paper_checks(data),
        "jacobi_checks": jacobi_checks(data),
        "stress_tensor_c2_precheck": stress_tensor_c2_precheck(data),
        "free_field_c2_exclusion": free_field_c2_exclusion(data),
    }

    output_path = Path(__file__).with_name("recompute_p3_results.json")
    output_path.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
