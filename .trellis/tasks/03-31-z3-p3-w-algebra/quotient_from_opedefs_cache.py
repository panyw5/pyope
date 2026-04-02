from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
CANDIDATE_CACHE = THIS_DIR / "search_weight8_t4_null_cache.json"
DESC_CACHE = THIS_DIR / "weight8_descendants_opedefs_cache.json"
OUTPUT_PATH = THIS_DIR / "weight8_known_singular_quotient_opedefs.json"
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

    def _leading_term(self, vector: dict[str, sp.Expr]):
        if not vector:
            return None
        monomial = max(vector)
        return monomial, vector[monomial]

    def reduce(self, vector: dict[str, sp.Expr]):
        reduced = {
            key: sp.sympify(value) for key, value in vector.items() if value != 0
        }
        while reduced:
            leading = self._leading_term(reduced)
            if leading is None:
                break
            pivot_monomial, pivot_coeff = leading
            pivot_vector = self._pivot_vectors.get(pivot_monomial)
            if pivot_vector is None:
                break
            factor = sp.sympify(pivot_coeff / pivot_vector[pivot_monomial])
            for monomial, coeff in pivot_vector.items():
                updated = sp.sympify(reduced.get(monomial, 0) - factor * coeff)
                if updated == 0:
                    reduced.pop(monomial, None)
                else:
                    reduced[monomial] = updated
        return reduced

    def insert(self, vector: dict[str, sp.Expr]) -> bool:
        reduced = self.reduce(vector)
        leading = self._leading_term(reduced)
        if leading is None:
            return False
        pivot_monomial, pivot_coeff = leading
        self._pivot_vectors[pivot_monomial] = {
            monomial: sp.sympify(coeff / pivot_coeff)
            for monomial, coeff in reduced.items()
        }
        return True

    @property
    def rank(self) -> int:
        return len(self._pivot_vectors)


def main():
    candidate_cache = json.loads(CANDIDATE_CACHE.read_text())
    descendant_cache = json.loads(DESC_CACHE.read_text())

    candidate_entries = candidate_cache["realized_entries"]
    descendant_entries = descendant_cache["entries"]

    candidate_vectors = [
        decode_candidate_vector(entry["sparse"]) for entry in candidate_entries
    ]
    descendant_vectors = [
        decode_opedefs_vector(payload["terms"])
        for _, payload in sorted(descendant_entries.items())
    ]

    candidate_elim = StringSparseEliminator()
    for vector in candidate_vectors:
        candidate_elim.insert(vector)
    candidate_rank = candidate_elim.rank

    descendant_elim = StringSparseEliminator()
    for vector in descendant_vectors:
        descendant_elim.insert(vector)
    descendant_rank = descendant_elim.rank

    combined_elim = StringSparseEliminator()
    for vector in descendant_vectors + candidate_vectors:
        combined_elim.insert(vector)
    combined_rank = combined_elim.rank

    target_index = next(
        i
        for i, entry in enumerate(candidate_entries)
        if entry["candidate"] == TARGET_STR
    )
    without_target_elim = StringSparseEliminator()
    for vector in descendant_vectors:
        without_target_elim.insert(vector)
    for i, vector in enumerate(candidate_vectors):
        if i == target_index:
            continue
        without_target_elim.insert(vector)
    rank_without_target = without_target_elim.rank
    rank_with_target = combined_rank

    payload = {
        "candidate_count": len(candidate_vectors),
        "candidate_rank": candidate_rank,
        "known_singular_descendant_count": len(descendant_vectors),
        "known_singular_descendant_rank": descendant_rank,
        "combined_rank": combined_rank,
        "quotient_dimension_after_known_singulars": combined_rank - descendant_rank,
        "target": candidate_entries[target_index]["candidate"],
        "target_survives_mod_known_singulars": rank_with_target > rank_without_target,
        "rank_without_target": rank_without_target,
        "rank_with_target": rank_with_target,
    }
    OUTPUT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
