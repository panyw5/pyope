from __future__ import annotations

import json
from pathlib import Path
import re

import sympy as sp


THIS_DIR = Path(__file__).resolve().parent
POLES_JSON = THIS_DIR / "weight4_t_poles_opedefs.json"
OUT_PATH = THIS_DIR / "weight4_quasiprimaries.json"


def normalize(text: str) -> str:
    text = text.replace("betaff", "beta")
    text = text.replace("gammaff", "gamma")
    text = text.replace("bff", "b")
    text = text.replace("cff", "c")
    text = re.sub(r"∂\^(\d+)([A-Za-z][A-Za-z0-9_]*)", r"Derivative[\1][\2]", text)
    text = re.sub(r"∂([A-Za-z][A-Za-z0-9_]*)", r"Derivative[1][\1]", text)
    text = text.replace("[", "(").replace("]", ")")
    text = text.replace(" ", "")
    return text


def vec(entries):
    out = {}
    for entry in entries:
        key = normalize(entry["monomial"])
        out[key] = sp.sympify(out.get(key, 0) + sp.sympify(entry["coeff"]))
    return {k: v for k, v in out.items() if v != 0}


def main():
    poles = json.loads(POLES_JSON.read_text())
    labels = list(poles.keys())
    pole4 = {label: vec(poles[label]["pole4"]) for label in labels}
    pole3 = {label: vec(poles[label]["pole3"]) for label in labels}

    mons4 = sorted({m for v in pole4.values() for m in v})
    mons3 = sorted({m for v in pole3.values() for m in v})

    rows = []
    for mon in mons4:
        row = [sp.sympify(pole4[label].get(mon, 0)) for label in labels]
        if any(x != 0 for x in row):
            rows.append(row)
    for mon in mons3:
        row = [sp.sympify(pole3[label].get(mon, 0)) for label in labels]
        if any(x != 0 for x in row):
            rows.append(row)

    M = sp.Matrix(rows)
    nullspace = M.nullspace()

    payload = {
        "candidate_count": len(labels),
        "constraint_rows": len(rows),
        "quasiprimary_dimension": len(nullspace),
        "basis": [
            [
                {"coeff": str(coeff), "label": label}
                for coeff, label in zip(vec, labels)
                if coeff != 0
            ]
            for vec in nullspace
        ],
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
