"""
Identify which expressions are slow to simplify.
"""
from __future__ import annotations
import sys, time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyope import BasicOperator, LocalOperatorBasis, MakeOPE, NO, One, OPE, clear_registry, d, make_realized
from pyope.wolfram_backend import op_to_wolfram_string, _eval_list_with_wolfram

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

print("Setting up...")
gens = setup_wz3()
basis = LocalOperatorBasis(gens)
exprs = [op.realize() for op in basis.list(weight=6)]

print(f"Total: {len(exprs)} expressions\n")

# Find top 10 most complex
complexities = [(i, e, len(op_to_wolfram_string(e))) for i, e in enumerate(exprs)]
complexities.sort(key=lambda x: x[2], reverse=True)

print("Top 10 most complex expressions:")
for rank, (idx, expr, comp) in enumerate(complexities[:10], 1):
    print(f"  {rank:2d}. expr[{idx:3d}]: complexity={comp:7d}  {str(expr)[:80]}...")

print("\nTesting individual simplification time for top 5...")
for rank, (idx, expr, comp) in enumerate(complexities[:5], 1):
    t0 = time.perf_counter()
    try:
        result = _eval_list_with_wolfram([expr], "CANONICALIZE_LIST")
        elapsed = time.perf_counter() - t0
        print(f"  {rank}. expr[{idx:3d}]: {elapsed:6.2f}s  (complexity={comp})")
    except Exception as e:
        elapsed = time.perf_counter() - t0
        print(f"  {rank}. expr[{idx:3d}]: FAILED after {elapsed:.2f}s - {e}")
