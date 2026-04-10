# Jacobi 恒等式教程 Notebook

## 文件信息

- **文件名**: `jacobi_identity_demo.ipynb`
- **位置**: `/Users/lelouch/pyope/demo/jacobi_identity_demo.ipynb`
- **创建日期**: 2026-01-07
- **状态**: ✅ 已完成并验证

## 内容概览

这是一个全面的教学性 Jupyter notebook，深入讲解顶点算符代数中的 Jacobi 恒等式。

### 主要章节

1. **数学背景**：Jacobi 恒等式的定义和物理意义
2. **环境设置**：导入必要的 pyope 模块
3. **Virasoro 代数示例**：
   - 定义能动张量 T 和其 OPE
   - 使用 `check_jacobi_identity` 检查恒等式
   - 可视化结果矩阵
4. **流代数示例**：U(1) 流代数的 Jacobi 恒等式验证
5. **交互式探索**：
   - 改变中心荷参数
   - 展示错误 OPE 定义的失败案例
6. **API 使用指南**：
   - `check_jacobi_identity` vs `verify_jacobi_identity`
   - `simplify_func` 参数的使用
   - 结果矩阵的解读
7. **理论背景**：Jacobi 恒等式的推导和物理意义
8. **练习与探索**：混合算符和自定义算符的练习
9. **总结**：关键要点和进一步学习资源
10. **技术附录**：实现细节、性能考虑和已知限制

## 学习目标

完成本 notebook 后，你将能够：

- ✅ 理解 Jacobi 恒等式的数学定义和物理意义
- ✅ 使用 pyope 验证任意算符的 Jacobi 恒等式
- ✅ 解读 Jacobi 恒等式的验证结果
- ✅ 识别错误的 OPE 定义
- ✅ 应用 Jacobi 恒等式到实际的 VOA 计算中

## 前置知识

- 基本的顶点算符代数概念
- 算符积展开 (OPE) 的定义
- Python 和 SymPy 基础

## 使用方法

```bash
# 启动 Jupyter Notebook
cd /Users/lelouch/pyope/demo
jupyter notebook jacobi_identity_demo.ipynb
```

或者使用 JupyterLab:

```bash
jupyter lab jacobi_identity_demo.ipynb
```

## 验证状态

✅ **所有代码已验证**：
- Virasoro 代数 Jacobi 恒等式：通过
- 流代数 Jacobi 恒等式：通过
- 所有示例代码可执行

## 参考资料

- `tests/test_jacobi_virasoro.py` - 完整测试用例
- `tests/JACOBI_README.md` - 快速开始指南
- `src/pyope/jacobi.py` - 实现代码
- V. Kac, "Vertex Algebras for Beginners", Chapter 3

## 特色功能

- 📊 **可视化**：使用 matplotlib 展示结果矩阵
- 🔄 **交互式**：可修改参数观察结果变化
- 📝 **详细注释**：每个步骤都有清晰的中文解释
- ✨ **实用示例**：Virasoro 代数和流代数的完整演示
- 🎯 **练习题**：帮助巩固理解的实践练习

## 适用人群

- 顶点算符代数研究者
- 共形场论学习者
- pyope 库用户
- 数学物理专业学生

## 技术细节

- **Cells 数量**: 37
- **代码 cells**: 约 15 个
- **Markdown cells**: 约 22 个
- **预计完成时间**: 45-60 分钟
- **难度级别**: 中级

## 相关文件

- `demo/pyope_ope_demo.ipynb` - OPE 计算基础教程
- `demo/pyope_basic_demo.ipynb` - pyope 基础功能演示
- `tests/test_jacobi_virasoro.py` - 自动化测试

---

**创建者**: Claude (Anthropic)  
**版本**: 1.0  
**最后更新**: 2026-01-07
