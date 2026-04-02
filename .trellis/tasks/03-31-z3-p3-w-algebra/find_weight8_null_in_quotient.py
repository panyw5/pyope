from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
CANDIDATE_CACHE = THIS_DIR / "search_weight8_t4_null_cache.json"
DESC_CACHE = THIS_DIR / "weight8_descendants_opedefs_cache.json"
OUTPUT_PATH = THIS_DIR / "weight8_null_relation_in_quotient.json"
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


class ProvenanceEliminator:
    def __init__(self) -> None:
        self._pivot_vectors: dict[str, dict[str, sp.Expr]] = {}
        self._pivot_provenance: dict[str, dict[int, sp.Expr]] = {}

    def _leading_term(self, vector: dict[str, sp.Expr]):
        if not vector:
            return None
        monomial = max(vector)
        return monomial, vector[monomial]

    def reduce_with_provenance(
        self, vector: dict[str, sp.Expr], provenance: dict[int, sp.Expr] | None = None
    ):
        reduced = {
            key: sp.sympify(value) for key, value in vector.items() if value != 0
        }
        prov = {
            key: sp.sympify(value)
            for key, value in (provenance or {}).items()
            if value != 0
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
            pivot_prov = self._pivot_provenance[pivot_monomial]
            for source_index, coeff in pivot_prov.items():
                updated = sp.sympify(prov.get(source_index, 0) - factor * coeff)
                if updated == 0:
                    prov.pop(source_index, None)
                else:
                    prov[source_index] = updated
        return reduced, prov

    def insert_source(self, source_index: int, vector: dict[str, sp.Expr]) -> bool:
        reduced, prov = self.reduce_with_provenance(
            vector, {source_index: sp.Integer(1)}
        )
        leading = self._leading_term(reduced)
        if leading is None:
            return False
        pivot_monomial, pivot_coeff = leading
        self._pivot_vectors[pivot_monomial] = {
            monomial: sp.sympify(coeff / pivot_coeff)
            for monomial, coeff in reduced.items()
        }
        self._pivot_provenance[pivot_monomial] = {
            idx: sp.sympify(coeff / pivot_coeff) for idx, coeff in prov.items()
        }
        return True


def main():
    candidate_cache = json.loads(CANDIDATE_CACHE.read_text())
    descendant_cache = json.loads(DESC_CACHE.read_text())

    candidate_entries = candidate_cache["realized_entries"]
    descendant_entries = sorted(descendant_cache["entries"].items())

    target_index = next(
        i
        for i, entry in enumerate(candidate_entries)
        if entry["candidate"] == TARGET_STR
    )
    target_vector = decode_candidate_vector(candidate_entries[target_index]["sparse"])

    sources = []
    for expr_str, payload in descendant_entries:
        sources.append(
            {
                "kind": "descendant",
                "label": expr_str,
                "vector": decode_opedefs_vector(payload["terms"]),
            }
        )
    for i, entry in enumerate(candidate_entries):
        if i == target_index:
            continue
        sources.append(
            {
                "kind": "candidate",
                "label": entry["candidate"],
                "vector": decode_candidate_vector(entry["sparse"]),
            }
        )

    eliminator = ProvenanceEliminator()
    inserted = 0
    for idx, source in enumerate(sources):
        if eliminator.insert_source(idx, source["vector"]):
            inserted += 1

    reduced, provenance = eliminator.reduce_with_provenance(
        target_vector, {-1: sp.Integer(1)}
    )
    relation_terms = []
    for source_index, coeff in sorted(
        provenance.items(), key=lambda item: (item[0] != -1, item[0])
    ):
        if source_index == -1:
            relation_terms.append(
                {
                    "kind": "target",
                    "label": TARGET_STR,
                    "coeff": str(coeff),
                }
            )
        else:
            relation_terms.append(
                {
                    "kind": sources[source_index]["kind"],
                    "label": sources[source_index]["label"],
                    "coeff": str(coeff),
                }
            )

    payload = {
        "source_count": len(sources),
        "pivot_rank": inserted,
        "target_relation_found": reduced == {},
        "residual_support_size": len(reduced),
        "relation_term_count": len(relation_terms),
        "relation_terms": relation_terms[:200],
        "relation_terms_full_path": str(OUTPUT_PATH),
    }
    if reduced == {}:
        full_payload = {
            **payload,
            "relation_terms_full": relation_terms,
        }
    else:
        full_payload = {
            **payload,
            "residual_support_sample": [
                {"monomial": monomial, "coeff": str(coeff)}
                for monomial, coeff in list(sorted(reduced.items()))[:50]
            ],
        }
    OUTPUT_PATH.write_text(json.dumps(full_payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
