"""
测试费米子正规序 OPE 的符号问题
"""
import sys
sys.path.insert(0, 'src')

from pyope import BasisOperator, OPE, NO, Bosonic, Fermionic, One
from fractions import Fraction

# 定义算符
b = BasisOperator('b', bosonic=True, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=True, conformal_weight=Fraction(-1))
beta = BasisOperator('β', bosonic=False, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=False, conformal_weight=Fraction(-1, 2))

Bosonic(b, c)
Fermionic(beta, gamma)

# 定义基本 OPE
OPE[b, c] = OPE.make([One])
OPE[beta, gamma] = OPE.make([-One])

print("=" * 60)
print("测试费米子正规序 OPE 符号")
print("=" * 60)
print()

# 测试 1: OPE(NO(b,c), NO(b,c))
print("测试 1: OPE(NO(b,c), NO(b,c))")
ope_bc_bc = OPE(NO(b, c), NO(b, c))
pole2_bc = ope_bc_bc.pole(2)
print(f"  极点(2) = {pole2_bc}")
print(f"  期望: One (玻色子)")
print(f"  状态: {'✓' if pole2_bc == One else '✗'}")
print()

# 测试 2: OPE(NO(β,γ), NO(β,γ))
print("测试 2: OPE(NO(β,γ), NO(β,γ))")
ope_betagamma_betagamma = OPE(NO(beta, gamma), NO(beta, gamma))
pole2_betagamma = ope_betagamma_betagamma.pole(2)
print(f"  极点(2) = {pole2_betagamma}")
print(f"  期望: -One (费米子)")
print(f"  状态: {'✓' if pole2_betagamma == -One else '✗'}")
print()

# 测试 3: OPE(j0, j0)
print("测试 3: OPE(j0, j0)")
j0 = 2 * NO(b, c) + 3 * NO(beta, gamma)
ope_j0_j0 = OPE(j0, j0)
pole2_j0 = ope_j0_j0.pole(2)
print(f"  极点(2) = {pole2_j0}")
print(f"  期望: -5*One")
print(f"  计算: 4*OPE(NO(b,c), NO(b,c))_pole(2) + 9*OPE(NO(β,γ), NO(β,γ))_pole(2)")
print(f"       = 4*1 + 9*(-1) = -5")
print(f"  状态: {'✓' if pole2_j0 == -5*One else '✗'}")
print()

print("=" * 60)
print("总结")
print("=" * 60)
if pole2_bc == One and pole2_betagamma == -One and pole2_j0 == -5*One:
    print("✅ 所有测试通过！")
else:
    print("❌ 存在符号错误")
    print()
    print("详细结果:")
    print(f"  OPE(NO(b,c), NO(b,c))_pole(2) = {pole2_bc} (期望 One)")
    print(f"  OPE(NO(β,γ), NO(β,γ))_pole(2) = {pole2_betagamma} (期望 -One)")
    print(f"  OPE(j0, j0)_pole(2) = {pole2_j0} (期望 -5*One)")
