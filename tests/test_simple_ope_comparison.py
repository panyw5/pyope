"""
简单 OPE 测试：对比 pyope 与 Mathematica
"""

import sys
sys.path.insert(0, '../src')

from pyope import (
    BasisOperator,
    OPE,
    MakeOPE,
    One,
    NO,
    Fermionic,
)

# 定义自由场
b = BasisOperator("b", bosonic=False, conformal_weight=-1)
c = BasisOperator("c", bosonic=False, conformal_weight=2)
beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)

# 声明费米子
Fermionic(b, c, beta, gamma)

# 定义自由场 OPE
OPE[b, c] = MakeOPE([One])
OPE[beta, gamma] = MakeOPE([-One])

print("=== pyope 基本 OPE 测试 ===\n")

# 测试 1: :bc: 与 :bc: 的 OPE
print("1. OPE[:bc:, :bc:]:")
bc = NO(b, c)
ope_bc_bc = OPE(bc, bc)
print(f"  max_pole = {ope_bc_bc.max_pole}")
print(f"  pole(2) = {ope_bc_bc.pole(2)}")
print(f"  pole(1) = {ope_bc_bc.pole(1)}")
print()

# 测试 2: :βγ: 与 :βγ: 的 OPE
print("2. OPE[:βγ:, :βγ:]:")
bg = NO(beta, gamma)
ope_bg_bg = OPE(bg, bg)
print(f"  max_pole = {ope_bg_bg.max_pole}")
print(f"  pole(2) = {ope_bg_bg.pole(2)}")
print(f"  pole(1) = {ope_bg_bg.pole(1)}")
print()

# 测试 3: (2:bc: + 3:βγ:) 与 (2:bc: + 3:βγ:) 的 OPE
print("3. OPE[2:bc: + 3:βγ:, 2:bc: + 3:βγ:]:")
j0 = 2*bc + 3*bg
ope_j0_j0 = OPE(j0, j0)
print(f"  max_pole = {ope_j0_j0.max_pole}")
print(f"  pole(2) = {ope_j0_j0.pole(2)}")
print(f"  pole(1) = {ope_j0_j0.pole(1)}")
print()

print("=== 对比 ===")
print("Mathematica: OPE[:bc:, :bc:] pole(2) = One")
print(f"pyope:       OPE[:bc:, :bc:] pole(2) = {ope_bc_bc.pole(2)}")
print()
print("Mathematica: OPE[:βγ:, :βγ:] pole(2) = -One")
print(f"pyope:       OPE[:βγ:, :βγ:] pole(2) = {ope_bg_bg.pole(2)}")
print()
print("Mathematica: OPE[j0, j0] pole(2) = -5*One")
print(f"pyope:       OPE[j0, j0] pole(2) = {ope_j0_j0.pole(2)}")
