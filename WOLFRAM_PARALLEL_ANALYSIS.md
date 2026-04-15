# Wolfram Backend 并行化性能分析报告

## 执行摘要

对 W_Z3 代数 weight=6 (305 个表达式) 的 `simplify_with_wolfram` 进行了详细 profiling，发现：

- **并行化已接近理论极限**：8 workers 达到 50s，理论下限 38s（仅差 12s）
- **瓶颈不是启动开销**：wolframscript 启动仅占 9-11% 的时间
- **真正瓶颈**：少数超慢表达式（特别是 `NO(Wbar,NO(Wbar,NO(Wbar,Wbar)))`）

## 性能数据

### 不同 worker 配置的性能

| Workers | Wall Time | Speedup | Efficiency |
|---------|-----------|---------|------------|
| 1       | 146.02s   | 1.00x   | 100%       |
| 2       | 88.17s    | 1.66x   | 83%        |
| 4       | 55.30s    | 2.64x   | 66%        |
| 8       | 50.08s    | 2.92x   | 36%        |

**理论极限**：38s（最慢单个表达式的时间）

### 最慢的表达式

| 排名 | 索引 | 抽象形式 | 复杂度 | 时间 |
|------|------|----------|--------|------|
| 1 | 304 | `NO(Wbar,NO(Wbar,NO(Wbar,Wbar)))` | 733,966 字符 | 38.09s |
| 2 | 222 | `NO(GbarWbar,NO(J,NO(Wbar,Wbar)))` | 204,667 字符 | 17.28s |
| 3 | 196 | `NO(Gbar,NO(Wbar,NO(Wbar,Wbar)))` | 171,431 字符 | 6.28s |
| 4 | 211 | `NO(∂GbarWbar,NO(Wbar,Wbar))` | 100,040 字符 | 16.64s |
| 5 | 228 | `NO(GbarWbar,NO(∂Wbar,Wbar))` | 99,800 字符 | 15.89s |

**Top 5 表达式占总时间的 64%**（94s / 146s）

## 根本原因

### Wbar 的指数级展开

Wbar 定义包含 8 项，每项都是嵌套的 NO 乘积。计算 NO 乘积时：

```
NO(Wbar, Wbar)                      → ~8×8 = 64 项
NO(Wbar, NO(Wbar, Wbar))            → ~8×64 = 512 项  
NO(Wbar, NO(Wbar, NO(Wbar, Wbar))) → ~8×512 = 4096 项
```

每一项都需要计算 Wick contraction（b·c, β·γ 的配对），导致：
- 表达式长度爆炸（733,966 字符）
- Wolfram 简化算法接近最坏情况

### Chunking 分析

**原始 chunking（CHUNK_MAX_ITEMS=32）**：
- 10 chunks，大小 [32, 32, 32, 32, 32, 32, 32, 32, 32, 17]
- chunk 7: 47.3s（包含 expr[222]）
- chunk 10: 42.4s（包含 expr[304]）
- **偶然地**将两个最慢表达式分散到不同 chunk，可以并行

**更小的 chunking（CHUNK_MAX_ITEMS=16）**：
- 20 chunks
- 8 workers: 49.17s（vs 50.08s，仅快 0.9s）
- 1 worker: 158.73s（vs 146.02s，慢 12.7s）
- **结论**：更小的 chunk 增加启动开销，没有实质提升

## 改进建议

### 1. 短期（无需代码改动）

**保持当前配置**：
- `PYOPE_WL_MAX_WORKERS=4` 或 `8`
- `PYOPE_WL_CHUNK_MAX_ITEMS=32`（默认）
- 当前性能已接近理论极限

### 2. 中期（需要代码改动）

**智能 chunking**：
```python
def smart_chunk_exprs(exprs, target_chunks=10):
    # 1. 预计算复杂度
    complexities = [(i, e, len(op_to_wolfram_string(e))) for i, e in enumerate(exprs)]
    
    # 2. 识别超慢表达式（complexity > 50000）
    slow_exprs = [x for x in complexities if x[2] > 50000]
    fast_exprs = [x for x in complexities if x[2] <= 50000]
    
    # 3. 每个慢表达式单独一个 chunk
    chunks = [[exprs[idx]] for idx, _, _ in slow_exprs]
    
    # 4. 快表达式均分
    # ...
    
    return chunks
```

**预期收益**：5-10s

### 3. 长期（架构优化）

**a) 测试 SymPy backend**：
```python
with compute_backend("sympy"):
    result = simplify(expr)
```
SymPy 可能对深度嵌套的 NO 乘积有更好的性能。

**b) 预简化策略**：
在构造高阶 NO 乘积前，先简化中间结果。

**c) 持久化 Wolfram kernel**：
使用 `wolframclient` Python 库保持长连接，消除启动开销（~13s 总计）。

**d) 过滤不需要的项**：
如果研究不需要 Wbar^4 等高阶项，在生成 basis 后过滤：
```python
list_of_ops = [op for op in basis.list(weight=6) 
               if op not in [list_of_ops[304], ...]]
```

## 结论

当前的并行化实现是**高效的**，性能瓶颈在于：
1. 少数表达式的固有复杂度（Wbar^4 → 指数级展开）
2. Wolfram 简化算法在这些情况下的性能特征

进一步优化应该：
- 从算法层面解决（SymPy backend、预简化）
- 或从问题层面规避（过滤不需要的高阶项）

而不是继续调整并行化参数（已接近极限）。

---

生成时间：2025-01-XX
测试环境：macOS, W_Z3 algebra weight=6 (305 expressions)
