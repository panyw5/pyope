# LocalOperator 架构重构总结

## 概述

成功实现了统一的 `LocalOperator` 类型系统，使得所有算符（包括线性组合）都是 `LocalOperator` 类型。

## 架构设计

```
LocalOperator (抽象基类)
├── Operator (基础算符，继承 sp.Symbol)
│   ├── BasisOperator (基础算符/生成元)
│   ├── DerivativeOperator (导数算符)
│   ├── NormalOrderedOperator (正规序算符)
│   └── ConstantOperator (常数算符)
└── CompositeOperator (复合算符)
    ├── OperatorSum (算符和: A + B)
    └── OperatorProduct (标量积: c*A)
```

## 主要特性

### 1. 统一类型系统
- 所有算符都是 `LocalOperator` 的实例
- `T` 是 `LocalOperator`
- `2*T` 是 `LocalOperator` (OperatorProduct)
- `2*T + J` 是 `LocalOperator` (OperatorSum)

### 2. 自动运算符重载
```python
T = BasisOperator("T", bosonic=True)
J = BasisOperator("J", bosonic=True)

# 加法返回 OperatorSum
sum_op = T + J  # OperatorSum(T, J)

# 标量乘法返回 OperatorProduct
prod_op = 2 * T  # OperatorProduct(2, T)

# 复合运算
complex_op = 2*T + 3*J  # OperatorSum(OperatorProduct(2, T), OperatorProduct(3, J))
```

### 3. OPE 和 NO 自动处理
```python
# OPE 自动展开线性组合
OPE(T, 2*T + J)  # = 2*OPE(T, T) + OPE(T, J)

# NO 自动分配
NO(2*T, J)  # = 2*NO(T, J)
```

### 4. 向后兼容
- 保持与 SymPy 的集成
- 支持旧的 `sp.Add` 和 `sp.Mul` 类型（向后兼容）
- 所有现有代码无需修改

## 实现细节

### 文件结构
- `local_operator.py`: 新增，定义 `LocalOperator` 基类和复合算符
- `operators.py`: 修改，让 `Operator` 继承 `LocalOperator`，重载运算符
- `api.py`: 修改，更新 `OPE` 和 `NO` 函数处理新类型
- `__init__.py`: 修改，导出新类型

### 关键实现
1. **LocalOperator 基类**
   - 继承 `sp.Expr` 以保持 SymPy 兼容性
   - 定义抽象 `parity` 属性
   - 重载 `__add__`, `__mul__` 等运算符

2. **OperatorSum**
   - 自动展平嵌套和
   - 单项自动简化
   - 支持 `terms` 属性访问各项

3. **OperatorProduct**
   - 自动合并系数
   - 分配律：`c*(A+B)` → `c*A + c*B`
   - 支持 `coeff` 和 `operator` 属性

4. **SymPy 兼容性**
   - 使用私有属性 `_terms`, `_coeff` 存储数据
   - 通过 `@property` 暴露 `args` 属性
   - 确保 `args` 中的元素是 SymPy 对象

## 测试结果

所有测试通过：
- ✓ 基础算符创建
- ✓ 算术运算（加法、标量乘法）
- ✓ 类型检查 (`is_local_operator`)
- ✓ OPE 计算（包括线性组合）
- ✓ NO 计算（包括线性组合）
- ✓ 宇称计算
- ✓ 向后兼容性

## 使用示例

```python
from pyope import BasisOperator, OPE, NO, is_local_operator

# 创建算符
T = BasisOperator("T", bosonic=True)
J = BasisOperator("J", bosonic=True)

# 算术运算
expr = 2*T + 3*J
print(is_local_operator(expr))  # True
print(type(expr))  # <class 'pyope.local_operator.OperatorSum'>

# OPE 自动处理
OPE(T, 2*T + J)  # 自动展开

# NO 自动处理
NO(2*T, J)  # 自动分配
```

## 优势

1. **数学直觉**：`2*T + J` 是一个算符，而不是 SymPy 表达式
2. **类型安全**：可以用 `isinstance(x, LocalOperator)` 检查
3. **自动简化**：运算符重载自动处理线性性
4. **向后兼容**：不破坏现有代码
5. **可扩展**：易于添加新的算符类型

## 未来改进

1. 实现更多算符类型（如张量积）
2. 优化性能（缓存计算结果）
3. 添加更多代数运算（对易子、反对易子）
4. 完善 LaTeX 显示
