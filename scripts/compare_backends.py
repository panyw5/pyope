"""
Compare Wolfram vs SymPy backend for slow expressions.
"""
from __future__ import annotations
import sys, time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyope import (BasicOperator, LocalOperatorBasis, MakeOPE, NO, One, OPE, 
                   clear_registry, d, make_realized, compute_backend)
from pyope.simplify import simplify

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

print("Setting up W_Z3...")
gens = setup_wz3()
basis = LocalOperatorBasis(gens)
list_of_ops = [op for op in basis.list(weight=6)]
print(f"Generated {len(list_of_ops)} operators\n")

# Test the slowest expressions
test_indices = [304, 222, 196]  # Top 3 slowest

print("="*80)
print("BACKEND COMPARISON")
print("="*80)

for idx in test_indices:
    print(f"\n{'='*80}")
    print(f"Testing list_of_ops[{idx}]: {list_of_ops[idx]}")
    print(f"{'='*80}")
    
    expr = list_of_ops[idx].realize()
    
    # Test SymPy backend
    print("\n[SymPy backend]")
    try:
        t0 = time.perf_counter()
        with compute_backend("sympy"):
            result_sympy = simplify(expr)
        t_sympy = time.perf_counter() - t0
        print(f"  Time: {t_sympy:.2f}s")
        print(f"  Result length: {len(str(result_sympy))} chars")
    except Exception as e:
        t_sympy = None
        print(f"  FAILED: {e}")
    
    # Test Wolfram backend
    print("\n[Wolfram backend]")
    try:
        from pyope.wolfram_backend import _eval_list_with_wolfram
        t0 = time.perf_counter()
        result_wolfram = _eval_list_with_wolfram([expr], "CANONICALIZE_LIST")[0]
        t_wolfram = time.perf_counter() - t0
        print(f"  Time: {t_wolfram:.2f}s")
        print(f"  Result length: {len(str(result_wolfram))} chars")
    except Exception as e:
        t_wolfram = None
        print(f"  FAILED: {e}")
    
    # Compare
    if t_sympy and t_wolfram:
        speedup = t_wolfram / t_sympy
        print(f"\n  Speedup: SymPy is {speedup:.2f}x {'faster' if speedup > 1 else 'slower'} than Wolfram")
        
        # Check if results are equivalent
        try:
            diff = simplify(result_sympy - result_wolfram)
            if diff == 0:
                print(f"  ✓ Results are equivalent")
            else:
                print(f"  ✗ Results differ: {diff}")
        except:
            print(f"  ? Could not verify equivalence")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
