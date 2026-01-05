# 任务

将 `OPEdefs.m` 以及相关的 Mathematica 代码包重构成 `python` 库 `pyope`，用于进行顶点算符代数符号计算。

# 库结构
- `README.md` 包含我个人对 `pyope` 的功能初步构想，请你继续帮我完善它
- **关键参考资料**: `OPEdefs` 文件夹包含 `OPEdefs.m` 以及相关的 Mathematica 代码包
- **关键参考资料**: `voa` SKILL 中的 `voa-manual.md` 的 section `3.3 Implementation`
- `src` 文件夹包含开发中的 `python` 库源文件
- `tests` 文件夹包含用于测试库功能的 `python` 脚本
- `papers` 文件夹包含顶点算符代数计算相关的论文，可以作为构建测试框架的参考
- `demo` 文件夹包含给用户展示与测试功能的 Jupyter Notebooks

# subagent

- 当需要阅读比较长的论文以及较大的 `.nb` 文件时，请你用 `collaborating-with-gemini` 调用 `gemini` 完成，因为 `gemini` 拥有较大的上下文窗口
- 对于其它较为复杂的任务，根据实际情况调用 `subagent`

# 文件读取与写入
- 分次分批分段读取，每次只读取一小部分，避免出错
- 分次分批分段写入，每次只写入一小部分，避免出错

# 开发与测试
- 功能开发请逐一完成，每个功能都应当进行充分测试
- 测试时，请你准备两个文件：一个用 `wolfram-engine` 或者 `wolfram-script` 进行计算，另一个用 `pyope` 计算，对比两个文件的结果
- 构思测试时，请认真参考 `voa` SKILL 中的 `voa-manual`、`comutations` 中的 `.nb` 文件