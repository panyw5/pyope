from __future__ import annotations

import json
from pathlib import Path
import re


THIS_DIR = Path(__file__).resolve().parent
PYOPE_PATH = THIS_DIR / "compare_nt_terms_pyope.json"
OPEDEFS_PATH = THIS_DIR / "compare_nt_terms_opedefs.json"
OUTPUT_PATH = THIS_DIR / "compare_nt_terms_diff.json"


def normalize_opedefs_text(text: str) -> str:
    text = text.replace("Î²", "beta")
    text = text.replace("Î³", "gamma")
    text = text.replace("[", "(")
    text = text.replace("]", ")")
    text = text.replace(" ", "")
    return text


def load_map(entries, normalizer):
    data = {}
    for entry in entries:
        key = normalizer(entry["monomial"])
        data[key] = entry["coeff"]
    return data


def main():
    pyope = json.loads(PYOPE_PATH.read_text())
    opedefs = json.loads(OPEDEFS_PATH.read_text())

    summary = {}
    for term_name in sorted(pyope):
        p_map = load_map(
            pyope[term_name],
            lambda s: s.replace("[", "(").replace("]", ")").replace(" ", ""),
        )
        o_map = load_map(opedefs[term_name], normalize_opedefs_text)
        only_p = sorted(set(p_map) - set(o_map))
        only_o = sorted(set(o_map) - set(p_map))
        coeff_mismatches = []
        for key in sorted(set(p_map) & set(o_map)):
            if p_map[key] != o_map[key]:
                coeff_mismatches.append(
                    {
                        "monomial": key,
                        "pyope": p_map[key],
                        "opedefs": o_map[key],
                    }
                )
        summary[term_name] = {
            "pyope_term_count": len(p_map),
            "opedefs_term_count": len(o_map),
            "matches_exactly": not only_p and not only_o and not coeff_mismatches,
            "only_in_pyope": only_p[:20],
            "only_in_opedefs": only_o[:20],
            "coeff_mismatches": coeff_mismatches[:20],
        }

    OUTPUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
