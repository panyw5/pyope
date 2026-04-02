from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import time

import sympy as sp

from pyope import (
    Zero,
    extract_scalar_operator,
    get_operator_parity,
    realize_and_simplify,
)
from pyope.operators import BasisOperator, DerivativeOperator, NormalOrderedOperator
from pyope.operator_spaces import (
    LocalOperatorBasis,
    RealizedGenerator,
)


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(charscan)

CACHE_PATH = THIS_DIR / "search_weight8_t4_null_cache.json"
RESULT_PATH = THIS_DIR / "search_weight8_t4_null_results.json"
WEIGHT = Fraction(8)
TARGET_STR = "NO(SG0,NO(SG0,NO(SG0,SG0)))"

FREE_FIELD_CHARGES = {
    "b_char_p3": -2,
    "c_char_p3": 2,
    "beta_char_p3": 3,
    "gamma_char_p3": -3,
}


def free_field_j_charge(expr):
    if isinstance(expr, BasisOperator) and not isinstance(expr, RealizedGenerator):
        return FREE_FIELD_CHARGES.get(expr.name, 0)
    if isinstance(expr, DerivativeOperator):
        return free_field_j_charge(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return free_field_j_charge(expr.left) + free_field_j_charge(expr.right)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return free_field_j_charge(op)
    if isinstance(expr, sp.Add):
        charges = {free_field_j_charge(arg) for arg in expr.args}
        return next(iter(charges)) if len(charges) == 1 else None
    return 0


def monomial_j_charge(expr, gen_charges):
    if isinstance(expr, RealizedGenerator):
        return gen_charges.get(expr.name)
    if isinstance(expr, DerivativeOperator):
        return monomial_j_charge(expr.base, gen_charges)
    if isinstance(expr, NormalOrderedOperator):
        q_left = monomial_j_charge(expr.left, gen_charges)
        q_right = monomial_j_charge(expr.right, gen_charges)
        return None if q_left is None or q_right is None else q_left + q_right
    if isinstance(expr, sp.Mul):
        _, inner = extract_scalar_operator(expr)
        return monomial_j_charge(inner, gen_charges)
    return None


def is_c2_monomial(expr):
    if isinstance(expr, DerivativeOperator):
        return True
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return is_c2_monomial(op)
    if isinstance(expr, NormalOrderedOperator):
        return isinstance(expr.left, DerivativeOperator)
    return False


def encode_sparse_dict(vector):
    return [[str(monomial), str(coeff)] for monomial, coeff in vector.items()]


def decode_sparse_dict(entries):
    return {
        monomial_repr: sp.sympify(coeff_str) for monomial_repr, coeff_str in entries
    }


class StringSparseEliminator:
    def __init__(self) -> None:
        self._pivot_vectors: dict[str, dict[str, sp.Expr]] = {}
        self._pivot_order: list[str] = []

    def _leading_term(self, vector: dict[str, sp.Expr]):
        if not vector:
            return None
        monomial = max(vector)
        return monomial, vector[monomial]

    def reduce_with_coefficients(self, vector: dict[str, sp.Expr]):
        reduced = {
            key: sp.sympify(value) for key, value in vector.items() if value != 0
        }
        coefficients: dict[str, sp.Expr] = {}

        while reduced:
            leading = self._leading_term(reduced)
            if leading is None:
                break
            pivot_monomial, pivot_coeff = leading
            pivot_vector = self._pivot_vectors.get(pivot_monomial)
            if pivot_vector is None:
                break
            factor = sp.sympify(pivot_coeff / pivot_vector[pivot_monomial])
            coefficients[pivot_monomial] = sp.sympify(
                coefficients.get(pivot_monomial, 0) + factor
            )
            for monomial, coeff in pivot_vector.items():
                updated = sp.sympify(reduced.get(monomial, 0) - factor * coeff)
                if updated == 0:
                    reduced.pop(monomial, None)
                else:
                    reduced[monomial] = updated

        return reduced, {
            key: value for key, value in coefficients.items() if value != 0
        }

    def reduce(self, vector: dict[str, sp.Expr]):
        reduced, _ = self.reduce_with_coefficients(vector)
        return reduced

    def insert_reduced(self, reduced: dict[str, sp.Expr]) -> bool:
        leading = self._leading_term(reduced)
        if leading is None:
            return False
        pivot_monomial, pivot_coeff = leading
        normalized = {
            monomial: sp.sympify(coeff / pivot_coeff)
            for monomial, coeff in reduced.items()
        }
        self._pivot_vectors[pivot_monomial] = normalized
        self._pivot_order.append(pivot_monomial)
        return True


def save_cache(payload):
    CACHE_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")


def load_cache():
    if not CACHE_PATH.exists():
        return None
    return json.loads(CACHE_PATH.read_text())


def build_problem():
    raw_map, _, free_fields = charscan.build_p3_data()
    raw_generators, discovered = charscan.singular_ope_closure(
        [raw_map[name] for name in ["T", "J", "G", "Gbar", "W", "Wbar"]],
        free_fields,
        WEIGHT,
    )
    generators = charscan.promote_named_generators(raw_generators)
    gen_basis = LocalOperatorBasis(generators, max_weight=WEIGHT)
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()),
        max_weight=WEIGHT,
        max_occurence=12,
    )
    gen_charges = {g.name: free_field_j_charge(g.realization) for g in generators}
    all_ops = gen_basis.list(weight=WEIGHT)
    candidates = [
        op
        for op in all_ops
        if get_operator_parity(op) == 0
        and monomial_j_charge(op, gen_charges) == 0
        and not is_c2_monomial(op)
    ]
    target_index = next(i for i, op in enumerate(candidates) if str(op) == TARGET_STR)
    return {
        "generators": generators,
        "discovered": discovered,
        "ff_basis": ff_basis,
        "candidates": candidates,
        "target_index": target_index,
    }


def realize_with_cache(candidates, ff_basis):
    cache = load_cache() or {}
    realized_entries = cache.get("realized_entries", [])
    start_index = len(realized_entries)
    print(f"cached realized entries: {start_index}/{len(candidates)}", flush=True)

    for index in range(start_index, len(candidates)):
        t0 = time.time()
        realized_expr = realize_and_simplify(candidates[index])
        sparse = ff_basis.sparse_terms(realized_expr)
        realized_entries.append(
            {
                "candidate": str(candidates[index]),
                "zero": realized_expr in (0, Zero),
                "sparse": encode_sparse_dict(sparse),
            }
        )
        if (index + 1) % 5 == 0 or time.time() - t0 > 5:
            save_cache({"realized_entries": realized_entries})
            print(
                f"realized {index + 1}/{len(candidates)} in cache update",
                flush=True,
            )

    save_cache({"realized_entries": realized_entries})
    return realized_entries


def solve_target_relation(candidates, target_index, realized_entries):
    eliminator = StringSparseEliminator()
    pivot_sources = {}
    source_data = {}
    target_vector = None
    target_expression = str(candidates[target_index])

    for index, entry in enumerate(realized_entries):
        sparse = decode_sparse_dict(entry["sparse"])
        if index == target_index:
            target_vector = sparse
            continue
        if not sparse:
            continue
        reduced, previous_coeffs = eliminator.reduce_with_coefficients(sparse)
        leading = eliminator._leading_term(reduced)
        if leading is None:
            continue
        pivot_sources[leading[0]] = index
        source_data[index] = {
            "pivot": leading[0],
            "leading_coeff": sp.sympify(leading[1]),
            "previous_coeffs": previous_coeffs,
        }
        eliminator.insert_reduced(reduced)

    if target_vector is None:
        raise ValueError("Target vector was not computed")

    reduced_target, coefficients = eliminator.reduce_with_coefficients(target_vector)

    expansion_cache = {}

    def pivot_to_source_coeffs(pivot):
        source_index = pivot_sources[pivot]
        if source_index in expansion_cache:
            return expansion_cache[source_index]
        data = source_data[source_index]
        result = {source_index: sp.sympify(1 / data["leading_coeff"])}
        for prev_pivot, coeff in data["previous_coeffs"].items():
            prev_expansion = pivot_to_source_coeffs(prev_pivot)
            factor = sp.sympify(-coeff / data["leading_coeff"])
            for prev_index, prev_coeff in prev_expansion.items():
                result[prev_index] = sp.sympify(
                    result.get(prev_index, 0) + factor * prev_coeff
                )
        result = {idx: coeff for idx, coeff in result.items() if coeff != 0}
        expansion_cache[source_index] = result
        return result

    source_coefficients = {}
    for pivot, coeff in coefficients.items():
        expansion = pivot_to_source_coeffs(pivot)
        for source_index, source_coeff in expansion.items():
            source_coefficients[source_index] = sp.sympify(
                source_coefficients.get(source_index, 0) + coeff * source_coeff
            )
    source_coefficients = {
        index: coeff for index, coeff in source_coefficients.items() if coeff != 0
    }

    result = {
        "target": target_expression,
        "target_index": target_index,
        "reduced_target_zero": reduced_target == {},
        "residual_support_size": len(reduced_target),
        "residual_support": encode_sparse_dict(reduced_target),
        "pivot_coefficients": [
            {
                "pivot_monomial": str(pivot),
                "source_index": pivot_sources[pivot],
                "source_operator": str(candidates[pivot_sources[pivot]]),
                "coeff": str(coefficients[pivot]),
            }
            for pivot in sorted(coefficients, key=sp.srepr)
        ],
    }

    if reduced_target == {}:
        normalized_terms = [{"op": target_expression, "coeff": "1"}]
        for source_index in sorted(source_coefficients):
            normalized_terms.append(
                {
                    "op": str(candidates[source_index]),
                    "coeff": str(sp.simplify(-source_coefficients[source_index])),
                }
            )
        result["normalized_relation"] = normalized_terms
        result["source_coefficients"] = [
            {
                "source_index": index,
                "source_operator": str(candidates[index]),
                "coeff": str(source_coefficients[index]),
            }
            for index in sorted(source_coefficients)
        ]

    return result


def main():
    problem = build_problem()
    generators = problem["generators"]
    discovered = problem["discovered"]
    candidates = problem["candidates"]
    ff_basis = problem["ff_basis"]
    target_index = problem["target_index"]

    print(f"closed generator count: {len(generators)}", flush=True)
    print(f"closure discoveries: {len(discovered)}", flush=True)
    print(f"weight-8 non-C2 bosonic J=0 candidates: {len(candidates)}", flush=True)
    print(f"target index: {target_index}", flush=True)

    realized_entries = realize_with_cache(candidates, ff_basis)
    solution = solve_target_relation(candidates, target_index, realized_entries)

    summary = {
        "closed_generator_count": len(generators),
        "closure_discoveries": discovered,
        "candidate_count": len(candidates),
        "target": str(candidates[target_index]),
        "target_index": target_index,
        "cache_entries": len(realized_entries),
        **solution,
    }
    RESULT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
