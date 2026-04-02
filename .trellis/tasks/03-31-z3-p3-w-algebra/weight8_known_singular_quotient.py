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
    normal_product,
    qp,
    realize_and_simplify,
)
from pyope.operator_spaces import (
    LocalOperatorBasis,
    RealizedGenerator,
    _get_conformal_weight,
)
from pyope.operators import BasisOperator, DerivativeOperator, NormalOrderedOperator, d


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

spec = importlib.util.spec_from_file_location(
    "charscan", THIS_DIR / "character_scan_p3.py"
)
charscan = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(charscan)

SEARCH_CACHE = THIS_DIR / "search_weight8_t4_null_cache.json"
DESC_CACHE = THIS_DIR / "weight8_known_singular_desc_cache.json"
RESULT_PATH = THIS_DIR / "weight8_known_singular_quotient.json"
WEIGHT = Fraction(8)
TARGET_STR = "NO(SG0,NO(SG0,NO(SG0,SG0)))"
SKIP_EXPRESSIONS = {
    "20*NO(FG0,NO(∂FG1,NO(FG1,NO(FG1,NO(FG1,FG1)))))/9 + 2*NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG2,∂FG3)))))/9 + NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG4,∂FG5))))) + 2*NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(∂FG2,FG3)))))/9 + NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(∂FG4,FG5))))) + 4*NO(NO(∂FG1,NO(FG1,NO(FG1,FG1))),NO(FG4,FG5)) + 4*NO(∂FG0,NO(FG1,NO(FG1,NO(FG1,NO(FG1,FG1)))))/9 - 7*NO(∂FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,FG1))))))/27 + 8*NO(∂FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG2,FG3)))))/9",
    "4*NO(FG0,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,FG1))))))/9 - NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,FG1)))))))/27 + 2*NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG2,FG3))))))/9 + NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG1,NO(FG4,FG5))))))",
}

FREE_FIELD_CHARGES = {
    "b_char_p3": -2,
    "c_char_p3": 2,
    "beta_char_p3": 3,
    "gamma_char_p3": -3,
}


def ff_charge(expr):
    if isinstance(expr, BasisOperator) and not isinstance(expr, RealizedGenerator):
        return FREE_FIELD_CHARGES.get(expr.name, 0)
    if isinstance(expr, DerivativeOperator):
        return ff_charge(expr.base)
    if isinstance(expr, NormalOrderedOperator):
        return ff_charge(expr.left) + ff_charge(expr.right)
    if isinstance(expr, sp.Mul):
        _, op = extract_scalar_operator(expr)
        return ff_charge(op)
    if isinstance(expr, sp.Add):
        charges = {ff_charge(arg) for arg in expr.args}
        return next(iter(charges)) if len(charges) == 1 else None
    return 0


def monomial_charge(expr, gen_charges):
    if isinstance(expr, RealizedGenerator):
        return gen_charges.get(expr.name)
    if isinstance(expr, DerivativeOperator):
        return monomial_charge(expr.base, gen_charges)
    if isinstance(expr, NormalOrderedOperator):
        q_left = monomial_charge(expr.left, gen_charges)
        q_right = monomial_charge(expr.right, gen_charges)
        return None if q_left is None or q_right is None else q_left + q_right
    if isinstance(expr, sp.Mul):
        _, inner = extract_scalar_operator(expr)
        return monomial_charge(inner, gen_charges)
    if isinstance(expr, sp.Add):
        values = {monomial_charge(arg, gen_charges) for arg in expr.args}
        return next(iter(values)) if len(values) == 1 else None
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
    return {monomial: sp.sympify(coeff) for monomial, coeff in entries}


def sparse_rank(vectors):
    if not vectors:
        return 0
    monomials = sorted({m for vec in vectors for m in vec})
    if not monomials:
        return 0
    matrix = sp.Matrix(
        [
            [sp.sympify(vec.get(monomial, 0)) for vec in vectors]
            for monomial in monomials
        ]
    )
    return int(matrix.rank())


def build_problem():
    raw_map, _, free_fields = charscan.build_p3_data()
    raw_generators, closure = charscan.singular_ope_closure(
        [raw_map[name] for name in ["T", "J", "G", "Gbar", "W", "Wbar"]],
        free_fields,
        WEIGHT,
    )
    generators = charscan.promote_named_generators(raw_generators)
    formal_generators = [
        RealizedGenerator(f"FG{i}", realization=expr)
        for i, expr in enumerate(raw_generators)
    ]
    formal_basis = LocalOperatorBasis(formal_generators, max_weight=WEIGHT)
    ff_basis = LocalOperatorBasis(
        list(free_fields.values()), max_weight=WEIGHT, max_occurence=12
    )
    gen_charges = {g.name: ff_charge(g.realization) for g in formal_generators}

    candidate_cache = json.loads(SEARCH_CACHE.read_text())
    candidate_entries = candidate_cache["realized_entries"]
    candidate_names = [entry["candidate"] for entry in candidate_entries]
    target_index = candidate_names.index(TARGET_STR)

    T, J, G, Gbar, W, Wbar, GW, GWbar = formal_generators
    JJ0 = qp(J, J, 0)
    X_higgs = (
        qp(W, Wbar, 0)
        - Fraction(1, 27) * qp(J, JJ0, 0, right_weight=2)
        + Fraction(2, 9) * qp(G, Gbar, 0)
        + Fraction(4, 9) * qp(J, T, 0, right_weight=2)
    )
    X_minus = qp(G, W, 0) - Fraction(1, 3) * qp(J, GW, 0, right_weight=2)
    X_plus = qp(Gbar, Wbar, 0) + Fraction(1, 3) * qp(J, GWbar, 0, right_weight=2)

    sources = {
        "X_higgs": X_higgs,
        "X_minus": X_minus,
        "X_plus": X_plus,
    }

    descendants = {}
    for label, source in sources.items():
        source_weight = _get_conformal_weight(source)
        source_parity = get_operator_parity(source)
        source_charge = monomial_charge(source, gen_charges)
        desc = set()
        for k in range(int(WEIGHT - source_weight) + 1):
            phi_weight = WEIGHT - source_weight - k
            for phi in formal_basis.list(phi_weight):
                if get_operator_parity(phi) != source_parity:
                    continue
                if monomial_charge(phi, gen_charges) != -source_charge:
                    continue
                expr = formal_basis.canonicalize(d(normal_product(phi, source), k))
                if expr != Zero and _get_conformal_weight(expr) == WEIGHT:
                    desc.add(expr)
        descendants[label] = sorted(desc, key=sp.srepr)

    combined_descendants = sorted(
        set().union(*[set(items) for items in descendants.values()]), key=sp.srepr
    )

    return {
        "closure": closure,
        "ff_basis": ff_basis,
        "candidate_entries": candidate_entries,
        "target_index": target_index,
        "sources": sources,
        "descendants": descendants,
        "combined_descendants": combined_descendants,
    }


def realize_descendants_with_cache(descendants, ff_basis):
    cache = (
        json.loads(DESC_CACHE.read_text()) if DESC_CACHE.exists() else {"entries": []}
    )
    entries = cache.get("entries", [])
    cached_names = {entry["expr"] for entry in entries}
    todo = [
        expr
        for expr in descendants
        if str(expr) not in cached_names and str(expr) not in SKIP_EXPRESSIONS
    ]
    print(f"cached descendant entries: {len(entries)}/{len(descendants)}", flush=True)
    for index, expr in enumerate(todo, start=1):
        DESC_CACHE.write_text(
            json.dumps(
                {
                    "entries": entries,
                    "in_progress": str(expr),
                    "remaining": len(todo) - index + 1,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n"
        )
        print(f"starting descendant {index}/{len(todo)}: {expr}", flush=True)
        t0 = time.time()
        realized = realize_and_simplify(expr)
        sparse = ff_basis.sparse_terms(realized)
        entries.append(
            {
                "expr": str(expr),
                "zero": realized in (0, Zero),
                "sparse": encode_sparse_dict(sparse),
            }
        )
        if index % 1 == 0 or time.time() - t0 > 2:
            DESC_CACHE.write_text(
                json.dumps(
                    {"entries": entries, "last_completed": str(expr)},
                    indent=2,
                    ensure_ascii=True,
                )
                + "\n"
            )
            print(f"desc realized +{index}/{len(todo)}", flush=True)
    DESC_CACHE.write_text(
        json.dumps({"entries": entries, "completed": True}, indent=2, ensure_ascii=True)
        + "\n"
    )
    entry_map = {
        entry["expr"]: decode_sparse_dict(entry["sparse"]) for entry in entries
    }
    return entry_map


def main():
    problem = build_problem()
    desc_map = realize_descendants_with_cache(
        problem["combined_descendants"], problem["ff_basis"]
    )
    candidate_vectors = [
        decode_sparse_dict(entry["sparse"]) for entry in problem["candidate_entries"]
    ]
    desc_vectors = [
        desc_map[str(expr)]
        for expr in problem["combined_descendants"]
        if desc_map.get(str(expr))
    ]

    source_stats = {}
    for label, exprs in problem["descendants"].items():
        vecs = [desc_map[str(expr)] for expr in exprs if desc_map.get(str(expr))]
        source_stats[label] = {
            "descendant_count": len(exprs),
            "realized_rank": sparse_rank(vecs),
            "sample_descendants": [str(expr) for expr in exprs[:20]],
        }

    source_rank = sparse_rank(desc_vectors)
    candidate_rank = sparse_rank(candidate_vectors)
    combined_rank = sparse_rank(candidate_vectors + desc_vectors)

    target_index = problem["target_index"]
    rank_without_target = sparse_rank(
        desc_vectors + [v for i, v in enumerate(candidate_vectors) if i != target_index]
    )
    rank_with_target = sparse_rank(desc_vectors + candidate_vectors)

    payload = {
        "candidate_count": len(candidate_vectors),
        "candidate_rank": candidate_rank,
        "known_singular_descendants": source_stats,
        "known_singular_descendant_count": len(problem["combined_descendants"]),
        "known_singular_descendant_rank": source_rank,
        "combined_rank": combined_rank,
        "quotient_dimension_after_known_singulars": combined_rank - source_rank,
        "target": problem["candidate_entries"][target_index]["candidate"],
        "target_survives_mod_known_singulars": rank_with_target > rank_without_target,
        "rank_without_target": rank_without_target,
        "rank_with_target": rank_with_target,
    }
    RESULT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
