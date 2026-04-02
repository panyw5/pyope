from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
CANDIDATE_CACHE = THIS_DIR / "search_weight8_t4_null_cache.json"
DESC_CACHE = THIS_DIR / "weight8_descendants_opedefs_cache.json"
OUTPUT_PATH = THIS_DIR / "target_vs_descendants.json"
TARGET_STR = "NO(SG0,NO(SG0,NO(SG0,SG0)))"


def normalize_common(text: str) -> str:
    text = text.replace("betaff", "beta")
    text = text.replace("gammaff", "gamma")
    text = text.replace("bff", "b")
    text = text.replace("cff", "c")
    text = text.replace("beta_char_p3", "beta")
    text = text.replace("gamma_char_p3", "gamma")
    text = text.replace("b_char_p3", "b")
    text = text.replace("c_char_p3", "c")
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    text = text.replace("[", "(").replace("]", ")")
    text = text.replace(" ", "")
    return text


def decode_candidate_vector(entries):
    vector = {}
    for monomial, coeff in entries:
        vector[normalize_common(monomial)] = sp.sympify(coeff)
    return {key: value for key, value in vector.items() if value != 0}


def decode_opedefs_vector(entries):
    vector = {}
    for entry in entries:
        vector[normalize_common(entry["monomial"])] = sp.sympify(entry["coeff"])
    return {key: value for key, value in vector.items() if value != 0}


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

    def insert(self, vector: dict[str, sp.Expr]) -> bool:
        reduced, _ = self.reduce_with_coefficients(vector)
        leading = self._leading_term(reduced)
        if leading is None:
            return False
        pivot_monomial, pivot_coeff = leading
        self._pivot_vectors[pivot_monomial] = {
            monomial: sp.sympify(coeff / pivot_coeff)
            for monomial, coeff in reduced.items()
        }
        self._pivot_order.append(pivot_monomial)
        return True


def main():
    candidate_cache = json.loads(CANDIDATE_CACHE.read_text())
    descendant_cache = json.loads(DESC_CACHE.read_text())
    candidate_entries = candidate_cache["realized_entries"]
    target_entry = next(
        entry for entry in candidate_entries if entry["candidate"] == TARGET_STR
    )
    target_vector = decode_candidate_vector(target_entry["sparse"])

    eliminator = StringSparseEliminator()
    descendant_count = 0
    for _, payload in sorted(descendant_cache["entries"].items()):
        vector = decode_opedefs_vector(payload["terms"])
        if vector:
            eliminator.insert(vector)
            descendant_count += 1

    reduced, coeffs = eliminator.reduce_with_coefficients(target_vector)
    payload = {
        "descendant_vector_count": descendant_count,
        "target_in_descendant_span": reduced == {},
        "residual_support_size": len(reduced),
        "residual_support_sample": [
            {"monomial": monomial, "coeff": str(coeff)}
            for monomial, coeff in list(sorted(reduced.items()))[:20]
        ],
        "pivot_coeff_count": len(coeffs),
    }
    OUTPUT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
