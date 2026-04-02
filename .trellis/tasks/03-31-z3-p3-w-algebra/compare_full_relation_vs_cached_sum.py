from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
FULL_REL = THIS_DIR / "full_weight8_relation_opedefs.json"
CAND_CACHE = THIS_DIR / "weight8_candidates_opedefs_cache.json"
RELATION_JSON = THIS_DIR / "weight8_null_relation_in_quotient.json"
OUT_PATH = THIS_DIR / "compare_full_relation_vs_cached_sum.json"


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


def vector_from_opedefs_terms(entries):
    out = {}
    for entry in entries:
        key = normalize_common(entry["monomial"])
        out[key] = sp.sympify(out.get(key, 0) + sp.sympify(entry["coeff"]))
    return {k: v for k, v in out.items() if v != 0}


def main():
    full_rel = json.loads(FULL_REL.read_text())
    cand_cache = json.loads(CAND_CACHE.read_text())
    relation = json.loads(RELATION_JSON.read_text())

    full_vec = vector_from_opedefs_terms(full_rel["terms"])
    cand_entries = cand_cache["entries"]

    sum_vec = {}
    for term in relation["relation_terms_full"]:
        label = term["label"]
        coeff = sp.sympify(term["coeff"])
        if label not in cand_entries:
            raise KeyError(f"Missing candidate in OPEdefs cache: {label}")
        vec = vector_from_opedefs_terms(cand_entries[label]["terms"])
        for monomial, value in vec.items():
            sum_vec[monomial] = sp.sympify(sum_vec.get(monomial, 0) + coeff * value)
    sum_vec = {k: v for k, v in sum_vec.items() if v != 0}

    all_keys = set(full_vec) | set(sum_vec)
    diff = {}
    for key in all_keys:
        value = sp.sympify(full_vec.get(key, 0) - sum_vec.get(key, 0))
        if value != 0:
            diff[key] = value

    payload = {
        "full_relation_term_count": len(full_vec),
        "cached_sum_term_count": len(sum_vec),
        "difference_term_count": len(diff),
        "difference_sample": [
            {"monomial": key, "coeff": str(diff[key])} for key in sorted(diff)[:50]
        ],
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
