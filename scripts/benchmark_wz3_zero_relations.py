from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from pyope import (  # noqa: E402
    LocalOperatorBasis,
    RealizedGenerator,
    compute_backend,
    list_zero_relations,
    realize,
)
from tests.utils.w_z3_algebra import make_z3_free_field_data


def build_wz3_realized_generators(prefix: str):
    ops = make_z3_free_field_data(prefix=prefix)
    realized = {
        name: RealizedGenerator(name, realization=ops[name])
        for name in ["T", "J", "W", "Wbar", "G", "Gbar", "GW", "GbarWbar"]
    }
    return ops, realized


def run_once(weight: int, max_workers: int, repeats: int) -> dict[str, object]:
    os.environ["PYOPE_WL_MAX_WORKERS"] = str(max_workers)

    timings: list[float] = []
    relation_count = 0
    expression_count = 0
    basis_size = 0

    for index in range(repeats):
        ops, realized = build_wz3_realized_generators(
            prefix=f"bench_w{weight}_r{index}"
        )
        abstract_basis = LocalOperatorBasis(realized.values())
        basis = LocalOperatorBasis(
            [ops["b"], ops["c"], ops["beta"], ops["gamma"]],
            max_occurence=max(6, weight + 1),
        )

        abstract_expressions = abstract_basis.list(weight)
        expressions = [realize(expr) for expr in abstract_expressions]

        expression_count = len(expressions)
        basis_size = len(basis.list(weight, max_occurence=max(6, weight + 1)))

        started = time.perf_counter()
        with compute_backend("wolfram"):
            zero_relations = list_zero_relations(expressions, basis, weight=weight)
        timings.append(time.perf_counter() - started)
        relation_count = len(zero_relations)

    return {
        "abstract_expression_count": len(abstract_expressions),
        "free_field_generators": ["b", "c", "beta", "gamma"],
        "weight": weight,
        "max_workers": max_workers,
        "repeats": repeats,
        "expression_count": expression_count,
        "basis_size": basis_size,
        "relation_count": relation_count,
        "timings": timings,
        "best": min(timings),
        "median": sorted(timings)[len(timings) // 2],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", nargs="+", type=int, default=[6, 7])
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 4])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    results = []
    for weight in args.weights:
        for workers in args.workers:
            results.append(run_once(weight, workers, args.repeats))

    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
