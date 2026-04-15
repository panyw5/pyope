"""
Profile wolfram parallelism: measure per-chunk timing and total wall time.
Usage: python scripts/profile_wolfram_parallel.py [--weight 6] [--workers 1 2 4 8]
"""
from __future__ import annotations
import argparse, os, sys, time, threading, statistics
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_call_log = []
_call_lock = threading.Lock()

import pyope.wolfram_backend as _wbe
_orig_invoke = _wbe._invoke_wolfram

def _patched_invoke(script_path, payload, operator_names):
    t0 = time.perf_counter()
    result = _orig_invoke(script_path, payload, operator_names)
    elapsed = time.perf_counter() - t0
    with _call_lock:
        _call_log.append({"elapsed": elapsed, "op": payload.get("PYOPE_WL_OPERATION", "?")})
    return result

_wbe._invoke_wolfram = _patched_invoke

from pyope import BasicOperator, LocalOperatorBasis, MakeOPE, NO, One, OPE, clear_registry, d, make_realized
from pyope.operator_spaces import _canonicalize_and_cache
from pyope.wolfram_backend import chunk_exprs_for_wolfram, _get_chunk_max_items, simplify_with_wolfram

def setup_wz3():
    clear_registry()
    b = BasicOperator("b", fermionic=True, conformal_weight=Fraction(2))
    c = BasicOperator("c", fermionic=True, conformal_weight=Fraction(-1))
    beta = BasicOperator("beta", conformal_weight=Fraction(3, 2))
    gamma = BasicOperator("gamma", conformal_weight=Fraction(-1, 2))
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])
    J = 2*NO(b,c) + 3*NO(beta,gamma)
    G = NO(gamma, b)
    Gbar = 2*NO(d(beta),c) + 3*NO(beta,d(c))
    T = (-2*NO(b,d(c)) - Fraction(3,2)*NO(beta,d(gamma)) - NO(d(b),c) - Fraction(1,2)*NO(d(beta),gamma))
    W = beta; GW = b
    Wbar = (NO(beta,NO(beta,NO(gamma,NO(gamma,gamma)))) + 2*NO(beta,NO(gamma,NO(gamma,NO(b,c)))) - 4*NO(beta,NO(d(gamma),gamma)) - Fraction(4,3)*NO(gamma,NO(b,d(c))) + Fraction(2,3)*NO(gamma,NO(d(b),c)) + Fraction(2,3)*NO(d(beta),NO(gamma,gamma)) - Fraction(8,3)*NO(d(gamma),NO(b,c)) + Fraction(10,9)*d(d(gamma)))
    GbarWbar = (Fraction(8,3)*NO(b,NO(d(d(c)),c)) + 3*NO(beta,NO(beta,NO(gamma,NO(gamma,d(c))))) - 4*NO(beta,NO(gamma,NO(b,NO(d(c),c)))) - 4*NO(beta,NO(gamma,d(d(c)))) - 4*NO(beta,NO(d(gamma),d(c))) - Fraction(2,3)*NO(d(b),NO(d(c),c)) + 2*NO(d(beta),NO(beta,NO(gamma,NO(gamma,c)))) - Fraction(8,3)*NO(d(beta),NO(d(gamma),c)) + Fraction(2,3)*NO(d(d(beta)),NO(gamma,c)) + Fraction(10,9)*d(d(d(c))))
    return list(make_realized([T, J, W, Wbar, G, Gbar, GW, GbarWbar]))

def run_simplify(exprs, workers):
    os.environ["PYOPE_WL_MAX_WORKERS"] = str(workers)
    _canonicalize_and_cache.cache_clear()
    with _call_lock: _call_log.clear()
    chunks = chunk_exprs_for_wolfram(exprs, max_items=_get_chunk_max_items())
    actual_workers = min(workers, len(chunks))
    t0 = time.perf_counter()
    simplify_with_wolfram(exprs)
    wall = time.perf_counter() - t0
    with _call_lock: logs = list(_call_log)
    pts = [e["elapsed"] for e in logs]
    return {"workers": workers, "actual_workers": actual_workers, "n_chunks": len(chunks),
            "chunk_sizes": [len(c) for c in chunks], "wall_time": wall, "n_calls": len(pts),
            "min": min(pts) if pts else 0, "max": max(pts) if pts else 0,
            "mean": statistics.mean(pts) if pts else 0, "median": statistics.median(pts) if pts else 0,
            "sum": sum(pts), "all": pts}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--weight", type=int, default=6)
    p.add_argument("--workers", nargs="+", type=int, default=[1, 2, 4, 8])
    args = p.parse_args()
    print(f"\n=== W_Z3 profile  weight={args.weight} ===", flush=True)
    gens = setup_wz3()
    basis = LocalOperatorBasis(gens)
    exprs = [op.realize() for op in basis.list(weight=args.weight)]
    chunks0 = chunk_exprs_for_wolfram(exprs, max_items=_get_chunk_max_items())
    print(f"Expressions: {len(exprs)}  CHUNK_MAX_ITEMS={_get_chunk_max_items()}")
    print(f"Chunks: {len(chunks0)}  sizes={[len(c) for c in chunks0]}\n")
    hdr = f"{'W':>4}  {'act':>3}  {'wall':>7}  {'calls':>5}  {'min':>6}  {'med':>6}  {'max':>6}  {'sum':>8}  {'eff':>5}  {'spdup':>6}"
    print(hdr); print("-"*len(hdr))
    results = []; baseline = None
    for w in args.workers:
        print(f"  workers={w}...", end="", flush=True)
        r = run_simplify(exprs, w)
        results.append(r)
        sp = baseline / r["wall_time"] if baseline else 1.0
        if baseline is None: baseline = r["wall_time"]
        eff = sp / max(1, r["actual_workers"])
        print(f"\r{r['workers']:>4}  {r['actual_workers']:>3}  {r['wall_time']:>7.2f}  {r['n_calls']:>5}  {r['min']:>6.2f}  {r['median']:>6.2f}  {r['max']:>6.2f}  {r['sum']:>8.2f}  {eff:>4.0%}  {sp:>5.2f}x")
    r1 = results[0]
    print(f"\n=== Diagnosis ===")
    print(f"  Startup floor (min chunk time) : {r1['min']:.3f}s")
    print(f"  Mean chunk time               : {r1['mean']:.3f}s")
    print(f"  Startup fraction (lower bound): {r1['min']/r1['mean']:.0%}")
    print(f"  Ideal wall (perfect parallel) : {r1['max']:.2f}s")
    print(f"  Actual serial wall            : {r1['wall_time']:.2f}s")
    print(f"  Max achievable speedup        : {r1['wall_time']/r1['max']:.1f}x")
    print(f"\n  Per-chunk detail (workers=1):")
    for i, t in enumerate(r1["all"]):
        bar = "#" * int(t / r1["max"] * 25)
        print(f"    chunk {i+1:2d} ({r1['chunk_sizes'][i]:2d}): {t:6.3f}s  {bar}")

if __name__ == "__main__":
    main()
