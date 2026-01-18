import sympy as sp
from pyope import BasisOperator, Fermionic, OPE, NO, d, One, MakeOPE

def check_system(name, lam, is_fermionic, ope_val):
    print(f"\n--- Checking {name} system (lambda={lam}, fermionic={is_fermionic}) ---")
    b = BasisOperator("b", conformal_weight=lam, fermionic=is_fermionic)
    c = BasisOperator("c", conformal_weight=1-lam, fermionic=is_fermionic)
    if is_fermionic:
        Fermionic(b, c)
    OPE[b, c] = MakeOPE([ope_val * One])
    
    # Try different T forms
    # Form 1: T = (lam-1) b dc + lam db c (VOA.wls form)
    T1 = (lam-1) * NO(b, d(c)) + lam * NO(d(b), c)
    # Form 2: T = -lam b dc + (1-lam) db c (Polchinski-like)
    T2 = -lam * NO(b, d(c)) + (1-lam) * NO(d(b), c)
    # Form 3: User's form
    T3 = -lam * NO(b, d(c)) - (1-lam) * NO(d(b), c)

    for i, T in enumerate([T1, T2, T3], 1):
        print(f"Testing Form {i}...")
        res_b = OPE(T, b)
        res_c = OPE(T, c)
        print(f"  T(z)b(w) pole 2: {res_b.pole(2)}, pole 1: {res_b.pole(1)}")
        print(f"  T(z)c(w) pole 2: {res_c.pole(2)}, pole 1: {res_c.pole(1)}")
        
        # Check central charge
        res_tt = OPE(T, T)
        print(f"  c/2 (pole 4): {res_tt.pole(4)}")

print("Reference: VOA.wls uses T = (lam-1) b dc + lam db c and gets h[b]=1-lam, h[c]=lam")
check_system("bc", 2, True, 1)
check_system("beta-gamma", 1.5, False, -1)
