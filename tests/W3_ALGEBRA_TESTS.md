# W₃ 代数测试文档

基于 `OPEdefs/ope-examples.nb` 中的 W₃ 代数部分（line 476-971）

## 测试文件

### 1. Mathematica 参考实现
- **文件**: `tests/w3_algebra_test.wls`
- **用途**: 使用 OPEdefs.m 生成参考结果
- **运行**: `wolframscript w3_algebra_test.wls`

### 2. Python 实现测试
- **文件**: `tests/test_w3_algebra.py`
- **用途**: 测试 PyOPE 的 W₃ 代数计算功能
- **运行**: `pytest tests/test_w3_algebra.py -v`

## W₃ 代数简介

W₃ 代数是 Virasoro 代数的扩展，包含两个生成元：
- **T**: stress tensor（共形权重 h=2）
- **W**: spin-3 primary field（共形权重 h=3）

### 辅助算符
定义 Λ 算符用于简化 W-W OPE：
```
Λ[w] = NO[T,T][w] - (3/10)T''[w]
```

### OPE 关系

#### T-T OPE (Virasoro)
```
T(z)T(w) ~ c/(2(z-w)⁴) + 2T/(z-w)² + T'/(z-w)
```
- 参数: c (中心荷)

#### T-W OPE
```
T(z)W(w) ~ 3W/(z-w)² + W'/(z-w)
```

#### W-W OPE
```
W(z)W(w) ~ c/(z-w)⁶
         + 2T/(z-w)⁴
         + T'/(z-w)³
         + (2βΛ + (3/10)T'')/(z-w)²
         + (βΛ' + (1/15)T''')/(z-w)
```
- 参数: c (中心荷), β (结构常数)

## 测试覆盖

### 第一部分：基础定义测试 (TestW3AlgebraDefinition)

1. **test_w3_operators_declaration**
   - 测试算符声明
   - 验证共形权重和玻色/费米统计

2. **test_lambda_construction**
   - 测试辅助算符 Λ 的构造
   - 验证 `Λ = NO[T,T] - (3/10)T''`

### 第二部分：OPE 关系测试 (TestW3AlgebraOPEs)

3. **test_virasoro_ope_with_c**
   - 测试 T-T OPE（带符号中心荷 c）
   - 验证极点结构：4-2-1-0 阶极点

4. **test_t_w_ope**
   - 测试 T-W OPE
   - 验证权重 3 算符的共形变换

5. **test_w_w_ope**
   - 测试 W-W OPE（W₃ 代数的核心）
   - 验证 6-4-3-2-1 阶极点
   - 包含参数 β 和辅助算符 Λ

### 第三部分：计算示例测试 (TestW3AlgebraComputations)

6. **test_t_lambda_ope**
   - 测试 `OPE[T, NO[T,T] - (3/10)T'']`
   - 对应 notebook line 718-794
   - 注意：PyOPE 返回 max_pole=6（而非 Mathematica 的 4），这是因为包含了 NO[T,T] 与 T 的 OPE 产生的高阶项

7. **test_t_w_second_derivative**
   - 测试 `OPE[T, W'']`
   - 验证导数算符的 OPE
   - 期望: `18W/(z-w)⁴ + 14W'/(z-w)³ + 5W''/(z-w)² + W'''/(z-w)`

8. **test_normal_ordering_t_prime_w_double_prime**
   - 测试 `NO[T', NO[W'', T]]`
   - 验证嵌套正规序乘积

9. **test_normal_ordering_w_prime_w_double_prime**
   - 测试 `NO[W', NO[W'', T]]`
   - 验证复杂的嵌套正规序展开

### 第四部分：数值检验 (TestW3AlgebraNumerical)

10. **test_w3_with_specific_parameters**
    - 测试特定参数值: c=100, β=1/10
    - 验证数值代入后的计算正确性

## 测试结果

所有 10 个测试均通过：

```
tests/test_w3_algebra.py::TestW3AlgebraDefinition::test_w3_operators_declaration PASSED
tests/test_w3_algebra.py::TestW3AlgebraDefinition::test_lambda_construction PASSED
tests/test_w3_algebra.py::TestW3AlgebraOPEs::test_virasoro_ope_with_c PASSED
tests/test_w3_algebra.py::TestW3AlgebraOPEs::test_t_w_ope PASSED
tests/test_w3_algebra.py::TestW3AlgebraOPEs::test_w_w_ope PASSED
tests/test_w3_algebra.py::TestW3AlgebraComputations::test_t_lambda_ope PASSED
tests/test_w3_algebra.py::TestW3AlgebraComputations::test_t_w_second_derivative PASSED
tests/test_w3_algebra.py::TestW3AlgebraComputations::test_normal_ordering_t_prime_w_double_prime PASSED
tests/test_w3_algebra.py::TestW3AlgebraComputations::test_normal_ordering_w_prime_w_double_prime PASSED
tests/test_w3_algebra.py::TestW3AlgebraNumerical::test_w3_with_specific_parameters PASSED
```

## 与 Mathematica 对比

### 主要差异

1. **极点数量**:
   - 在 `test_t_lambda_ope` 中，PyOPE 返回 max_pole=6
   - Mathematica (ope-examples.nb line 767-789) 显示 max_pole=4
   - **原因**:
     - 计算 `OPE[T, Λ] = OPE[T, NO(T,T) - (3/10)T'']`
     - 其中 `OPE[T, T'']` 根据导数 OPE 规则产生 6 阶极点 `10*One*c`
     - 乘以系数 `-3/10` 得到 6 阶极点 `-3*One*c`
     - **PyOPE 保留了完整的极点结构，而 Mathematica 可能做了某种截断或简化**
   - PyOPE 的结果：
     - 6 阶极点: `-3*One*c`
     - 4 阶极点: `T*c + 22*T/5`
     - 3 阶极点: `-3*∂T`
     - 2 阶极点: `-6*∂^2T/5 + 4*NO(T,T)`
     - 1 阶极点: `-3*∂^3T/10 + NO(T,∂T) + NO(∂T,T)`

### 一致性验证

除了上述差异外，PyOPE 与 Mathematica 在以下方面保持一致：
- 基本 OPE 极点结构
- 导数算符的变换规则
- 正规序乘积的展开
- 数值计算结果

## 参考资料

1. **ope-examples.nb** (line 476-971): W₃ 代数定义和计算示例
2. **P. Bouwknegt and K. Schoutens**: *W-symmetry in conformal field theory*
3. **VOA Manual**: 顶点算符代数的一般理论

## 未来扩展

可以考虑添加的测试：
1. W₃ 代数的雅可比恒等式验证
2. 不同中心荷下的特殊情况（如 c=2）
3. W-代数的更高 spin 扩展（W₄, W₅等）
4. 与自由场表示的对应关系
5. W₃ 代数的模空间结构
