# 发布到 PyPI

本文档对应当前仓库 `pyope` 的发布流程。

## 1. 首次发布前检查

- 确认你已经注册 PyPI 账号：<https://pypi.org/account/register/>
- 确认 `pyproject.toml` 里的版本号已经更新
- 确认包名 `pyope-voa` 在 PyPI 上可用，或者该名字已经归你所有
- 建议先在 TestPyPI 试发一次

当前项目的关键元数据来自：

- 包名：`pyope-voa`
- 版本：`0.1.0.post1`
- 构建后端：`setuptools.build_meta`

## 2. 安装发布工具

```bash
python -m pip install --upgrade build twine
```

## 3. 本地构建

在仓库根目录运行：

```bash
python -m build
```

成功后会生成：

- `dist/pyope_voa-<version>.tar.gz`
- `dist/pyope_voa-<version>-py3-none-any.whl`

## 4. 检查分发包

```bash
python -m twine check dist/*
```

## 5. 上传到 TestPyPI（推荐先做）

```bash
python -m twine upload --repository testpypi dist/*
```

安装测试：

```bash
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pyope-voa
```

## 6. 上传到正式 PyPI

```bash
python -m twine upload dist/*
```

## 7. 推荐使用 API Token

上传时，用户名使用：

```text
__token__
```

密码使用你在 PyPI 里创建的 API Token。

也可以写入 `~/.pypirc`：

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = <your-pypi-token>

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = <your-testpypi-token>
```

## 8. 以后每次发版

典型流程：

1. 修改 `pyproject.toml` 中的 `version`
2. 同步修改 `src/pyope/__init__.py` 中的 `__version__`
3. 清理旧产物
4. 重新构建并检查
5. 上传到 PyPI

建议命令：

```bash
rm -rf dist build *.egg-info
python -m build
python -m twine check dist/*
python -m twine upload dist/*
```

## 9. 当前仓库还需要注意的点

- `pyproject.toml` 里的项目链接已经指向真实仓库 `https://github.com/panyw5/pyope`
- 当前仓库已补充 `LICENSE`
- 发版前请确认版本号不是已经在 PyPI 存在的旧版本
