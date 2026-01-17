# 测试文件统计性错误修复总结

## 修复的文件列表

### ✅ 已修复的文件（共 6 个）

1. **test_null_states_stage1.py**
   - 第 33 行：`Bosonic(b, c)` → `Fermionic(b, c)`
   - 第 98 行：`Fermionic(beta, gamma)` → `Bosonic(beta, gamma)`

2. **test_null_states_stage2.py**
   - 第 39-40 行：修正统计性声明

3. **test_null_states_verification.py**
   - 第 44-45 行：修正统计性声明

4. **test_null_states_grouped.py**
   - 第 38-39 行：修正统计性声明

5. **test_fermion_sign_bug.py**
   - 第 11-17 行：完全反了！修正为：
     - b, c: `bosonic=False` + `Fermionic(b, c)`
     - beta, gamma: `bosonic=True` + `Bosonic(beta, gamma)`

6. **test_level4_verification.py**
   - 第 34-35 行：修正统计性声明

## 统一的正确约定

```python
# bc 系统：费米子（鬼场）
b = BasisOperator('b', bosonic=False, conformal_weight=Fraction(2))
c = BasisOperator('c', bosonic=False, conformal_weight=Fraction(-1))
Fermionic(b, c)

# βγ 系统：玻色子
beta = BasisOperator('β', bosonic=True, conformal_weight=Fraction(3, 2))
gamma = BasisOperator('γ', bosonic=True, conformal_weight=Fraction(-1, 2))
Bosonic(beta, gamma)
```

## 验证状态

- ✅ test_null_states_stage1.py - 已运行，所有测试通过
- ⏳ 其他文件 - 待验证

## 其他发现的问题

### test_free_field_systems.py

该文件中的定义看起来不标准：
```python
# 第 36-37 行
self.beta = BasisOperator("β", bosonic=False, conformal_weight=-1/2)
self.gamma = BasisOperator("γ", bosonic=False, conformal_weight=3/2)

# 第 87-88 行
self.b = BasisOperator("b", bosonic=False, conformal_weight=-1)
self.c = BasisOperator("c", bosonic=False, conformal_weight=2)
```

**问题**：
1. 权重不标准（beta 应该是 3/2，gamma 应该是 -1/2）
2. b 和 c 的权重也不标准（b 应该是 2，c 应该是 -1）

**建议**：检查这个文件是否在测试特殊情况，或者需要修正。

## 下一步行动

1. ✅ 修复所有统计性声明错误（已完成）
2. ⏳ 运行所有测试验证修复
3. ⏳ 检查 test_free_field_systems.py 的意图
4. ⏳ 添加缺失的测试（费米子自身 NO 乘积等）
