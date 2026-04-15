"""
Measure the overhead breakdown of Python-Wolfram communication.
"""
from __future__ import annotations
import sys, time, tempfile
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyope import BasicOperator, LocalOperatorBasis, MakeOPE, NO, One, OPE, clear_registry, d, make_realized
from pyope.wolfram_backend import op_to_wolfram_string

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
list_of_ops = [op for op in basis.list(weight=6)]
expr = list_of_ops[304].realize()

print(f"\nTesting expr[304]: NO(Wbar,NO(Wbar,NO(Wbar,Wbar)))")
print("="*70)

# 1. Measure serialization
print("\n1. Python serialization (op_to_wolfram_string)")
t0 = time.perf_counter()
wolfram_str = op_to_wolfram_string(expr)
t_serialize = time.perf_counter() - t0
print(f"   Time: {t_serialize:.3f}s")
print(f"   Length: {len(wolfram_str):,} characters")

# 2. Measure file I/O
print("\n2. File I/O (write + read)")
with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
    t0 = time.perf_counter()
    f.write(wolfram_str)
    temp_path = f.name
    t_write = time.perf_counter() - t0

t0 = time.perf_counter()
with open(temp_path, 'r') as f:
    content = f.read()
t_read = time.perf_counter() - t0
Path(temp_path).unlink()

print(f"   Write: {t_write:.3f}s")
print(f"   Read:  {t_read:.3f}s")
print(f"   Total: {t_write + t_read:.3f}s")

# 3. Measure total _invoke_wolfram time
print("\n3. Total _invoke_wolfram (includes all overhead + computation)")
from pyope.wolfram_backend import _eval_list_with_wolfram
t0 = time.perf_counter()
result = _eval_list_with_wolfram([expr], "CANONICALIZE_LIST")
t_total = time.perf_counter() - t0
print(f"   Time: {t_total:.3f}s")

# 4. Calculate overhead
print("\n" + "="*70)
print("OVERHEAD BREAKDOWN")
print("="*70)
print(f"Total time (measured):           {t_total:.2f}s")
print(f"User's direct Mathematica time:  27.00s (from screenshot)")
print(f"Estimated overhead:              {t_total - 27:.2f}s ({(t_total-27)/t_total*100:.1f}%)")
print()
print(f"Known overhead components:")
print(f"  - Python serialization:        {t_serialize:.2f}s")
print(f"  - File I/O (write+read):       {t_write + t_read:.2f}s")
print(f"  - Kernel startup + OPEdefs:    ~2.0s (estimated)")
print(f"  - ToExpression parsing:        ~2.0s (estimated)")
print(f"  - Result serialization:        ~1.0s (estimated)")
print(f"  - Result parsing in Python:    ~1.0s (estimated)")
print(f"  Total estimated overhead:      ~{t_serialize + t_write + t_read + 6:.1f}s")
print()
print(f"Actual computation (estimated):  {t_total - (t_serialize + t_write + t_read + 6):.1f}s")
print(f"User's measurement:              27.0s")
