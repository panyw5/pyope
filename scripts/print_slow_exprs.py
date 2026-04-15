"""
Print the slowest expressions in detail.
"""
from __future__ import annotations
import sys
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

print("Setting up W_Z3 algebra...")
gens = setup_wz3()
basis = LocalOperatorBasis(gens)

print("\nGenerating weight-6 operators (abstract)...")
list_of_ops = [op for op in basis.list(weight=6)]
print(f"Total: {len(list_of_ops)} operators\n")

print("Realizing expressions...")
list_of_ops_realized = [op.realize() for op in list_of_ops]

# Find slowest by complexity
complexities = [(i, list_of_ops[i], list_of_ops_realized[i], len(op_to_wolfram_string(list_of_ops_realized[i]))) 
                for i in range(len(list_of_ops))]
complexities.sort(key=lambda x: x[3], reverse=True)

print("=" * 100)
print("TOP 10 SLOWEST EXPRESSIONS")
print("=" * 100)

for rank, (idx, abstract, realized, comp) in enumerate(complexities[:10], 1):
    print(f"\n{'='*100}")
    print(f"RANK {rank}: list_of_ops[{idx}]")
    print(f"{'='*100}")
    print(f"Complexity: {comp:,} characters")
    print(f"\nAbstract form:")
    print(f"  {abstract}")
    print(f"\nRealized form (first 500 chars):")
    realized_str = str(realized)
    if len(realized_str) > 500:
        print(f"  {realized_str[:500]}...")
        print(f"  ... (total {len(realized_str)} characters)")
    else:
        print(f"  {realized_str}")
    
    # Try to identify the structure
    abstract_str = str(abstract)
    if "Wbar" in abstract_str:
        wbar_count = abstract_str.count("Wbar")
        print(f"\n  Contains {wbar_count} Wbar factor(s)")
    if "GbarWbar" in abstract_str:
        gbwb_count = abstract_str.count("GbarWbar")
        print(f"  Contains {gbwb_count} GbarWbar factor(s)")

print("\n" + "="*100)
print("SUMMARY")
print("="*100)
print("\nTop 5 indices for reference:")
for rank, (idx, _, _, comp) in enumerate(complexities[:5], 1):
    print(f"  {rank}. list_of_ops[{idx}]  (complexity={comp:,})")
