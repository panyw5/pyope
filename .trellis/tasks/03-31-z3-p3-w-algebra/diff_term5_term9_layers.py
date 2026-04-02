from __future__ import annotations

import json
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent


def normalize(text: str) -> str:
    return text.replace("[", "(").replace("]", ")").replace(" ", "")


def load_map(entries):
    return {normalize(entry["monomial"]): entry["coeff"] for entry in entries}


def main():
    pyope = json.loads((THIS_DIR / "compare_term5_term9_pyope.json").read_text())
    opedefs = json.loads((THIS_DIR / "compare_term5_term9_opedefs.json").read_text())

    result = {}
    for key in [
        "Wbar",
        "dWbar",
        "NO_W_dWbar",
        "NO_dW_Wbar",
        "term5",
        "term9",
        "term5_plus_term9",
    ]:
        p_map = load_map(pyope[key])
        o_map = load_map(opedefs[key])
        only_p = sorted(set(p_map) - set(o_map))
        only_o = sorted(set(o_map) - set(p_map))
        coeff_mismatches = []
        for monomial in sorted(set(p_map) & set(o_map)):
            if p_map[monomial] != o_map[monomial]:
                coeff_mismatches.append(
                    {
                        "monomial": monomial,
                        "pyope": p_map[monomial],
                        "opedefs": o_map[monomial],
                    }
                )
        result[key] = {
            "pyope_terms": len(p_map),
            "opedefs_terms": len(o_map),
            "matches_exactly": not only_p and not only_o and not coeff_mismatches,
            "only_in_pyope": only_p[:20],
            "only_in_opedefs": only_o[:20],
            "coeff_mismatches": coeff_mismatches[:20],
        }

    out = THIS_DIR / "diff_term5_term9_layers.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
