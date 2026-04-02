from __future__ import annotations

import json
from pathlib import Path
import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
PREFIX_JSON = THIS_DIR / "prefix_relations_opedefs.json"
CAND_CACHE = THIS_DIR / "weight8_candidates_opedefs_cache.json"
REL_JSON = THIS_DIR / "weight8_null_relation_in_quotient.json"
OUT_PATH = THIS_DIR / "check_prefix_faithfulness.json"


def vec(entries):
    out = {}
    for entry in entries:
        key = entry["monomial"]
        out[key] = sp.sympify(out.get(key, 0) + sp.sympify(entry["coeff"]))
    return {k: v for k, v in out.items() if v != 0}


def main():
    prefix = json.loads(PREFIX_JSON.read_text())
    cand = json.loads(CAND_CACHE.read_text())["entries"]
    rel_terms = json.loads(REL_JSON.read_text())["relation_terms_full"]

    findings = []
    running = {}
    for i, term in enumerate(rel_terms, start=1):
        label = term["label"]
        coeff = sp.sympify(term["coeff"])
        for monomial, value in vec(cand[label]["terms"]).items():
            running[monomial] = sp.sympify(running.get(monomial, 0) + coeff * value)
        running = {k: v for k, v in running.items() if v != 0}

        direct = vec(prefix[str(i)])
        diff = {}
        for key in set(running) | set(direct):
            delta = sp.sympify(direct.get(key, 0) - running.get(key, 0))
            if delta != 0:
                diff[key] = delta
        findings.append(
            {
                "prefix_len": i,
                "direct_terms": len(direct),
                "cached_sum_terms": len(running),
                "diff_terms": len(diff),
                "added_term": label,
            }
        )

    OUT_PATH.write_text(json.dumps(findings, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(findings, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
