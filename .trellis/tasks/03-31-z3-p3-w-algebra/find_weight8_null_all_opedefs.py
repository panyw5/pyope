from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
CAND_CACHE = THIS_DIR / "weight8_candidates_opedefs_cache.json"
DESC_CACHE = THIS_DIR / "weight8_descendants_opedefs_cache.json"
OUT_PATH = THIS_DIR / "weight8_null_all_opedefs.json"
TARGET_STR = "NO(SG0,NO(SG0,NO(SG0,SG0)))"


def normalize_common(text: str) -> str:
    text = text.replace("betaff", "beta")
    text = text.replace("gammaff", "gamma")
    text = text.replace("bff", "b")
    text = text.replace("cff", "c")
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    text = text.replace("[", "(").replace("]", ")")
    text = text.replace(" ", "")
    return text


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

    def reduce_with_provenance(self, vector: dict[str, sp.Expr], provenance=None):
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
    cand_cache = json.loads(CAND_CACHE.read_text())
    desc_cache = json.loads(DESC_CACHE.read_text())

    sources = []
    for expr_str, payload in sorted(desc_cache["entries"].items()):
        sources.append(
            {
                "kind": "descendant",
                "label": expr_str,
                "vector": decode_opedefs_vector(payload["terms"]),
            }
        )

    target_vector = None
    for expr_str, payload in sorted(cand_cache["entries"].items()):
        vector = decode_opedefs_vector(payload["terms"])
        if expr_str == TARGET_STR:
            target_vector = vector
            continue
        sources.append(
            {
                "kind": "candidate",
                "label": expr_str,
                "vector": vector,
            }
        )

    if target_vector is None:
        raise ValueError("Target candidate not found in OPEdefs cache")

    eliminator = ProvenanceEliminator()
    pivot_rank = 0
    for idx, source in enumerate(sources):
        if eliminator.insert_source(idx, source["vector"]):
            pivot_rank += 1

    reduced, prov = eliminator.reduce_with_provenance(
        target_vector, {-1: sp.Integer(1)}
    )
    relation_terms = []
    for idx, coeff in sorted(prov.items(), key=lambda item: (item[0] != -1, item[0])):
        if idx == -1:
            relation_terms.append(
                {"kind": "target", "label": TARGET_STR, "coeff": str(coeff)}
            )
        else:
            relation_terms.append(
                {
                    "kind": sources[idx]["kind"],
                    "label": sources[idx]["label"],
                    "coeff": str(coeff),
                }
            )

    payload = {
        "source_count": len(sources),
        "pivot_rank": pivot_rank,
        "target_relation_found": reduced == {},
        "residual_support_size": len(reduced),
        "relation_term_count": len(relation_terms),
        "relation_terms": relation_terms[:200],
    }
    if reduced == {}:
        payload["relation_terms_full"] = relation_terms
    else:
        payload["residual_support_sample"] = [
            {"monomial": monomial, "coeff": str(coeff)}
            for monomial, coeff in list(sorted(reduced.items()))[:50]
        ]
    OUT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
