from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    NO,
    One,
    OPE,
    Zero,
    bracket,
    clear_registry,
    d,
    get_operator_parity,
    list_independent_ops,
    make_realized,
    realize_and_simplify,
    simplify,
)
from pyope.operator_spaces import (
    LocalOperatorBasis,
    RealizedGenerator,
    _get_conformal_weight,
)


EXPECTED_SERIES = {
    0: 1,
    2: 1,
    4: 1,
    6: 2,
    7: -2,
    8: 3,
    9: -2,
    10: 4,
    11: -4,
    12: 6,
    13: -6,
    14: 8,
    15: -8,
}


def build_p3_data():
    clear_registry()

    b = BasisOperator("b_char_p3", fermionic=True, conformal_weight=Fraction(2))
    c = BasisOperator("c_char_p3", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasisOperator("beta_char_p3", conformal_weight=Fraction(3, 2))
    gamma = BasisOperator("gamma_char_p3", conformal_weight=Fraction(-1, 2))

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

    raw_map = {
        "T": T,
        "W": W,
        "J": J,
        "G": G,
        "Gbar": Gbar,
        "GW": GW,
        "Wbar": Wbar,
        "GWbar": GWbar,
    }
    realized = make_realized(
        [raw_map[name] for name in ["T", "W", "J", "G", "Gbar", "GW", "Wbar", "GWbar"]]
    )
    names = ["T", "W", "J", "G", "Gbar", "GW", "Wbar", "GWbar"]
    realized_map = {name: value for name, value in zip(names, realized)}
    free_fields = {"b": b, "c": c, "beta": beta, "gamma": gamma}
    return raw_map, realized_map, free_fields


def weight_grid(max_weight: Fraction):
    current = Fraction(0)
    while current <= max_weight:
        yield current
        current += Fraction(1, 2)


def _nonzero_realized(expressions):
    return [
        expr
        for expr in (realize_and_simplify(op) for op in expressions)
        if expr not in (0, Zero)
    ]


def promote_named_generators(raw_generators, prefix: str = "SG"):
    return [
        RealizedGenerator(f"{prefix}{index}", realization=generator)
        for index, generator in enumerate(raw_generators)
    ]


def normalized_sparse_signature(expr, basis: LocalOperatorBasis):
    sparse = basis.sparse_terms(expr)
    if not sparse:
        return "Zero"
    pivot = max(sparse, key=sp.srepr)
    pivot_coeff = sp.sympify(sparse[pivot])
    normalized = {
        monomial: sp.simplify(coeff / pivot_coeff) for monomial, coeff in sparse.items()
    }
    items = sorted(normalized.items(), key=lambda item: sp.srepr(item[0]))
    return repr([(sp.srepr(monomial), str(coeff)) for monomial, coeff in items])


def _is_new_singular_atom(
    candidate, current_raw_generators, ff_basis, max_weight: Fraction
):
    candidate_weight = _get_conformal_weight(candidate)
    if (
        candidate_weight is None
        or candidate_weight <= 0
        or candidate_weight > max_weight
    ):
        return False

    candidate_realized = realize_and_simplify(candidate)
    if candidate_realized in (0, Zero):
        return False

    candidate_parity = get_operator_parity(candidate)
    current_realized_generators = promote_named_generators(
        current_raw_generators, prefix="CHK"
    )
    basis = LocalOperatorBasis(current_realized_generators, max_weight=max_weight)
    same_weight_ops = basis.list(weight=candidate_weight)
    same_parity_realized = _nonzero_realized(
        [op for op in same_weight_ops if get_operator_parity(op) == candidate_parity]
    )
    old_rank = len(
        list_independent_ops(
            same_parity_realized,
            ff_basis,
            weight=candidate_weight,
            max_occurence=12,
        )
    )
    new_rank = len(
        list_independent_ops(
            [*same_parity_realized, candidate_realized],
            ff_basis,
            weight=candidate_weight,
            max_occurence=12,
        )
    )
    return new_rank > old_rank


def singular_ope_closure(seed_raw_generators, free_fields, max_weight: Fraction):
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()),
        max_weight=max_weight,
        max_occurence=12,
    )
    raw_generators = list(seed_raw_generators)
    discovered = []
    seen_realizations = {
        normalized_sparse_signature(realize_and_simplify(generator), ff_basis)
        for generator in raw_generators
        if realize_and_simplify(generator) not in (0, Zero)
    }

    while True:
        additions = []
        current_snapshot = list(raw_generators)
        for left in current_snapshot:
            left_weight = _get_conformal_weight(left)
            if left_weight is None:
                continue
            for right in current_snapshot:
                right_weight = _get_conformal_weight(right)
                if right_weight is None:
                    continue
                max_pole = int(sp.floor(left_weight + right_weight))
                for pole in range(1, max_pole + 1):
                    candidate = simplify(bracket(left, right, pole))
                    if candidate in (0, Zero):
                        continue
                    realized_candidate = realize_and_simplify(candidate)
                    if realized_candidate in (0, Zero):
                        continue
                    signature = normalized_sparse_signature(
                        realized_candidate, ff_basis
                    )
                    if signature in seen_realizations:
                        continue
                    if not _is_new_singular_atom(
                        candidate,
                        raw_generators,
                        ff_basis,
                        max_weight,
                    ):
                        continue
                    additions.append(
                        {
                            "candidate": candidate,
                            "from": (str(left), str(right), pole),
                            "weight": str(_get_conformal_weight(candidate)),
                            "parity": get_operator_parity(candidate),
                        }
                    )
                    seen_realizations.add(signature)

        if not additions:
            break

        for item in additions:
            raw_generators.append(item["candidate"])
            discovered.append(
                {
                    "operator": str(item["candidate"]),
                    "from": item["from"],
                    "weight": item["weight"],
                    "parity": item["parity"],
                }
            )

    return raw_generators, discovered


def scan_character_from_raw_generators(
    raw_generators, free_fields, max_weight: Fraction
):
    generators = promote_named_generators(raw_generators)
    basis = LocalOperatorBasis(generators, max_weight=max_weight)
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()),
        max_weight=max_weight,
        max_occurence=12,
    )

    rows = []
    for weight in weight_grid(max_weight):
        print(f"[scan] generator_count={len(generators)} weight={weight}", flush=True)
        ops = basis.list(weight=weight)
        bosonic_ops = [op for op in ops if get_operator_parity(op) == 0]
        fermionic_ops = [op for op in ops if get_operator_parity(op) == 1]

        bosonic_realized = _nonzero_realized(bosonic_ops)
        fermionic_realized = _nonzero_realized(fermionic_ops)

        indep_bosonic_realized = list_independent_ops(
            bosonic_realized,
            ff_basis,
            weight=weight,
            max_occurence=12,
        )
        indep_fermionic_realized = list_independent_ops(
            fermionic_realized,
            ff_basis,
            weight=weight,
            max_occurence=12,
        )

        bosonic = len(indep_bosonic_realized)
        fermionic = len(indep_fermionic_realized)
        superdim = bosonic - fermionic
        rows.append(
            {
                "weight": str(weight),
                "double_weight": int(2 * weight),
                "candidate_count": len(ops),
                "independent_count": bosonic + fermionic,
                "bosonic_count": bosonic,
                "fermionic_count": fermionic,
                "superdimension": superdim,
            }
        )
    return rows


def scan_character(
    generator_names,
    raw_map,
    free_fields,
    max_weight: Fraction,
    auto_close_singular: bool = False,
):
    raw_generators = [raw_map[name] for name in generator_names]
    closure = []
    if auto_close_singular:
        raw_generators, closure = singular_ope_closure(
            raw_generators, free_fields, max_weight
        )
    rows = scan_character_from_raw_generators(raw_generators, free_fields, max_weight)
    return rows, closure, [str(generator) for generator in raw_generators]


def compare_with_expected(rows):
    comparison = []
    for row in rows:
        power = row["double_weight"]
        expected = EXPECTED_SERIES.get(power)
        comparison.append(
            {
                "double_weight": power,
                "weight": row["weight"],
                "computed": row["superdimension"],
                "expected": expected,
                "matches": expected == row["superdimension"]
                if expected is not None
                else None,
            }
        )
    return comparison


def main():
    max_weight = Fraction(4)
    raw_map, realized_map, free_fields = build_p3_data()

    rows6_naive, closure6_naive, gens6_naive = scan_character(
        ["T", "J", "G", "Gbar", "W", "Wbar"],
        raw_map,
        free_fields,
        max_weight,
        auto_close_singular=False,
    )
    rows6_closed, closure6_closed, gens6_closed = scan_character(
        ["T", "J", "G", "Gbar", "W", "Wbar"],
        raw_map,
        free_fields,
        max_weight,
        auto_close_singular=True,
    )
    rows8, closure8, gens8 = scan_character(
        ["T", "J", "G", "Gbar", "W", "Wbar", "GW", "GWbar"],
        raw_map,
        free_fields,
        max_weight,
        auto_close_singular=False,
    )

    payload = {
        "max_weight": str(max_weight),
        "expected_series": EXPECTED_SERIES,
        "six_generator_scan_naive": rows6_naive,
        "six_generator_comparison_naive": compare_with_expected(rows6_naive),
        "six_generator_scan_closed": rows6_closed,
        "six_generator_comparison_closed": compare_with_expected(rows6_closed),
        "six_generator_closed_atoms": gens6_closed,
        "six_generator_closure_discovered": closure6_closed,
        "eight_generator_scan": rows8,
        "eight_generator_comparison": compare_with_expected(rows8),
        "eight_generator_atoms": gens8,
    }

    output_path = Path(__file__).with_name("character_scan_p3_results.json")
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
