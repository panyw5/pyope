# Mathematica 测试文件使用指南

## 概述

本目录包含了为不同代数类型生成的 Mathematica 测试脚本，用于计算嵌套正规序（nested NO）乘积展开。这些测试结果将作为 Python pyope 实现的参考基准。

## 文件列表

### Mathematica 测试脚本

1. **test_virasoro_mathematica.wls** - Virasoro 代数测试
2. **test_sl2_mathematica.wls** - sl(2) Kac-Moody 代数测试
3. **test_u1_mathematica.wls** - U(1) 流代数测试
4. **test_bcbetagamma_mathematica.wls** - bc-βγ 自由场系统测试
5. **test_w3_mathematica.wls** - W₃ 代数测试
6. **test_n4_mathematica.wls** - N=4 小超共形代数测试

### 辅助脚本

- **run_all_mathematica_tests.sh** - 批量运行所有 Mathematica 测试
- **parse_mathematica_results.py** - 解析 Mathematica 输出并生成 Python 测试框架

### 文档

- **README_MATHEMATICA_TESTS.md** - 详细的测试文件说明文档

## 使用步骤

### 第一步：运行 Mathematica 测试

确保已安装 Wolfram Engine 或 Mathematica，然后运行：

```bash
cd /Users/lelouch/pyope/tests
./run_all_mathematica_tests.sh
```

这将：
- 执行所有 6 个代数类型的测试
- 将结果保存到 `mathematica_results/` 目录
- 每个代数生成一个独立的结果文件

### 第二步：检查结果

查看生成的结果文件：

```bash
ls -lh mathematica_results/
cat mathematica_results/virasoro_results.txt
```

结果文件格式：
```
========================================
Test: VIR-1.1: NO(NO(T, T), T)
Output: <展开后的表达式>

========================================
Test: VIR-1.2: NO(T, NO(T, T))
Output: <展开后的表达式>
...
```

### 第三步：解析结果并生成 Python 测试（可选）

运行解析脚本：

```bash
python3 parse_mathematica_results.py
```

这将：
- 解析所有 Mathematica 结果文件
- 生成 Python 测试框架到 `generated_tests/` 目录
- 每个代数生成一个 Python 测试文件

## 测试覆盖范围

### 1. Virasoro 代数 (15 个测试)
- 简单嵌套 NO
- 带导数的嵌套 NO
- 三重嵌套 NO
- 高阶导数
- 复杂嵌套结构

### 2. sl(2) Kac-Moody 代数 (22 个测试)
- 简单嵌套 NO（所有生成元组合）
- 三重嵌套 NO
- 带导数的嵌套 NO
- Sugawara 构造
- 复杂嵌套结构

### 3. U(1) 流代数 (26 个测试)
- J 的嵌套 NO
- T 和 J 的混合
- 带导数的嵌套 NO
- 混合 T、J 与导数
- 复杂嵌套结构

### 4. bc-βγ 系统 (30 个测试)
- bc 系统的嵌套 NO
- βγ 系统的嵌套 NO
- 带导数的嵌套 NO
- bc 和 βγ 的混合
- 高阶导数
- 复杂嵌套结构

### 5. W₃ 代数 (27 个测试)
- T 和 W 的嵌套 NO
- 仅 T 的嵌套 NO
- 带导数的嵌套 NO
- 三重嵌套 NO
- 高阶导数
- 复杂嵌套结构

### 6. N=4 小超共形代数 (28 个测试)
- J 生成元的嵌套 NO
- G 生成元的嵌套 NO
- J 和 G 的混合
- 带 T 的嵌套 NO
- 三重嵌套结构
- 带导数的复杂嵌套
- Sugawara 类构造

**总计：148 个测试用例**

## 代数定义

### Virasoro
```mathematica
OPE[T, T] = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
```

### sl(2)
```mathematica
OPE[J+, J0] = J+/(z-w)
OPE[J0, J0] = k/2/(z-w)^2
OPE[J0, J-] = -J-/(z-w)
OPE[J+, J-] = k/(z-w)^2 + 2J0/(z-w)
```

### U(1)
```mathematica
OPE[J, J] = 1/(z-w)^2
OPE[T, J] = J/(z-w)^2 + ∂J/(z-w)
OPE[T, T] = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
```

### bc-βγ
```mathematica
OPE[b, c] = 1/(z-w)
OPE[β, γ] = -1/(z-w)
```

### W₃
```mathematica
OPE[T, T] = c/2/(z-w)^4 + 2T/(z-w)^2 + ∂T/(z-w)
OPE[T, W] = 3W/(z-w)^2 + ∂W/(z-w)
OPE[W, W] = c/(z-w)^6 + 2T/(z-w)^4 + ... (复杂)
```

### N=4
从 bc-βγ 系统构造：
```mathematica
J+ = β
J0 = NO[b,c] + 2NO[β,γ]
J- = -NO[β,NO[γ,γ]] - NO[γ,NO[b,c]] + 3/2∂γ
G+ = b
G- = NO[b,γ]
T = -3/2NO[b,∂c] - NO[β,∂γ] - 1/2NO[∂b,c]
```

## 注意事项

1. **路径依赖**：所有测试脚本使用绝对路径 `/Users/lelouch/pyope/OPEdefs/OPEdefs.m`
   - 如果 OPEdefs.m 位置不同，需要修改每个 .wls 文件的第 8 行

2. **符号参数**：测试保留符号参数（c, k, β 等），不进行数值替换

3. **输出格式**：Mathematica 输出使用其标准格式
   - `Derivative[n][X]` 表示 n 阶导数
   - `NO[A, B]` 表示正规序乘积

4. **执行时间**：某些复杂测试可能需要较长时间
   - W₃ 和 N=4 测试尤其耗时
   - 建议单独运行或增加超时时间

## 故障排除

### 问题：wolframscript 命令未找到
**解决方案**：
```bash
# 安装 Wolfram Engine (免费)
# 或确保 Mathematica 已安装并添加到 PATH
export PATH="/Applications/Mathematica.app/Contents/MacOS:$PATH"
```

### 问题：OPEdefs.m 加载失败
**解决方案**：
```bash
# 检查文件是否存在
ls -l /Users/lelouch/pyope/OPEdefs/OPEdefs.m

# 如果路径不同，修改所有 .wls 文件的第 8 行
sed -i '' 's|/Users/lelouch/pyope/OPEdefs/OPEdefs.m|/your/path/OPEdefs.m|g' test_*_mathematica.wls
```

### 问题：测试运行时间过长
**解决方案**：
```bash
# 单独运行特定测试
wolframscript -file test_virasoro_mathematica.wls

# 或减少测试用例数量（编辑 .wls 文件）
```

## 下一步

1. 运行所有 Mathematica 测试生成参考结果
2. 在 Python 中实现对应的 `expand_nested_no` 功能
3. 编写 Python 测试用例，与 Mathematica 结果对比
4. 确保所有测试通过

## 贡献

如需添加新的代数类型测试：
1. 创建新的 `test_<algebra>_mathematica.wls` 文件
2. 定义代数的生成元和 OPE
3. 添加各种嵌套 NO 测试用例
4. 更新 `run_all_mathematica_tests.sh`
5. 更新本文档

## 参考

- OPEdefs 文档：包含在 `OPEdefs/OPEdefs.m` 头部注释
- VOA 计算示例：`.claude/skills/voa/computations/VOA.wls`
- 原始测试：`test_expand_nested_no_mathematica.wls`
