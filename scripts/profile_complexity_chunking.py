"""
Test complexity-aware chunking vs naive chunking.
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
        _call_log.append({"elapsed": elapsed})
    return result

_wbe._invoke_wolfram = _patched_invoke

from pyope import BasicOperator, LocalOperatorBasis, MakeOPE, NO, One, OPE, clear_registry, d, make_realized
from pyope.operator_spaces import _canonicalize_and_cache
from pyope.wolfram_backend import op_to_wolfram_string, simplify_with_wolfram
from concurrent.futures import ThreadPoolExecutor

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

def complexity_aware_chunks(exprs, target_chunks=10):
    """Split by complexity (string length) instead of count."""
    complexities = [(i, len(op_to_wolfram_string(e))) for i, e in enumerate(exprs)]
    complexities.sort(key=lambda x: x[1])
    
    total_complexity = sum(c for _, c in complexities)
    target_per_chunk = total_complexity / target_chunks
    
    chunks = []
    current_chunk = []
    current_complexity = 0
    
    for idx, comp in complexities:
        current_chunk.append(exprs[idx])
        current_complexity += comp
        if current_complexity >= target_per_chunk and len(chunks) < target_chunks - 1:
            chunks.append(current_chunk)
            current_chunk = []
            current_complexity = 0
    
    if current_chunk:
        chunks.append(current_chunk)
    
    return chunks

def _eval_list_with_wolfram_direct(exprs):
    """Direct call to wolfram for a list."""
    from pyope.wolfram_backend import _eval_list_with_wolfram
    return _eval_list_with_wolfram(exprs, "CANONICALIZE_LIST")

def run_with_custom_chunks(exprs, chunks, workers):
    """Run simplify with custom chunking."""
    os.environ["PYOPE_WL_MAX_WORKERS"] = str(workers)
    _canonicalize_and_cache.cache_clear()
    with _call_lock: _call_log.clear()
    
    actual_workers = min(workers, len(chunks))
    
    t0 = time.perf_counter()
    if actual_workers == 1:
        results = []
        for chunk in chunks:
            results.extend(_eval_list_with_wolfram_direct(chunk))
    else:
        with ThreadPoolExecutor(max_workers=actual_workers) as executor:
            chunk_results = list(executor.map(_eval_list_with_wolfram_direct, chunks))
        results = []
        for cr in chunk_results:
            results.extend(cr)
    
    wall = time.perf_counter() - t0
    
    with _call_lock: logs = list(_call_log)
    pts = [e["elapsed"] for e in logs]
    
    return {
        "workers": workers, "actual_workers": actual_workers,
        "n_chunks": len(chunks), "chunk_sizes": [len(c) for c in chunks],
        "wall_time": wall, "n_calls": len(pts),
        "min": min(pts) if pts else 0, "max": max(pts) if pts else 0,
        "mean": statistics.mean(pts) if pts else 0,
        "sum": sum(pts), "all": pts
    }

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--weight", type=int, default=6)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    
    print(f"\n=== Complexity-aware chunking test  weight={args.weight} ===\n")
    gens = setup_wz3()
    basis = LocalOperatorBasis(gens)
    exprs = [op.realize() for op in basis.list(weight=args.weight)]
    
    print(f"Total expressions: {len(exprs)}")
    print(f"Workers: {args.workers}\n")
    
    # Analyze complexity distribution
    complexities = [len(op_to_wolfram_string(e)) for e in exprs]
    print(f"Complexity stats:")
    print(f"  min: {min(complexities)}")
    print(f"  max: {max(complexities)}")
    print(f"  mean: {statistics.mean(complexities):.1f}")
    print(f"  median: {statistics.median(complexities):.1f}")
    print(f"  total: {sum(complexities)}\n")
    
    # Test complexity-aware chunking
    chunks_complex = complexity_aware_chunks(exprs, target_chunks=10)
    chunk_complexities = [sum(len(op_to_wolfram_string(e)) for e in c) for c in chunks_complex]
    
    print(f"Complexity-aware chunks: {len(chunks_complex)}")
    print(f"  sizes: {[len(c) for c in chunks_complex]}")
    print(f"  complexities: {chunk_complexities}")
    print(f"  complexity std: {statistics.stdev(chunk_complexities):.1f}")
    print(f"  complexity range: {max(chunk_complexities) - min(chunk_complexities)}")
    
    print(f"\nRunning with complexity-aware chunking...", flush=True)
    r_complex = run_with_custom_chunks(exprs, chunks_complex, args.workers)
    
    print(f"\n=== Results ===")
    print(f"Wall time: {r_complex['wall_time']:.2f}s")
    print(f"Chunks: {r_complex['n_chunks']}")
    print(f"Per-chunk times:")
    for i, t in enumerate(r_complex['all']):
        bar = "#" * int(t / r_complex['max'] * 30)
        print(f"  chunk {i+1:2d} ({r_complex['chunk_sizes'][i]:2d} exprs, complexity {chunk_complexities[i]:5d}): {t:6.3f}s  {bar}")
    
    print(f"\nTime distribution:")
    print(f"  min: {r_complex['min']:.2f}s")
    print(f"  max: {r_complex['max']:.2f}s")
    print(f"  mean: {r_complex['mean']:.2f}s")
    print(f"  range: {r_complex['max'] - r_complex['min']:.2f}s")
    print(f"  std: {statistics.stdev(r_complex['all']):.2f}s")

if __name__ == "__main__":
    main()
