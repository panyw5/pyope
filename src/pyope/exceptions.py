"""
异常类定义模块

本模块定义了 pyope 库中使用的自定义异常类型。
"""


class IllegalOperatorProductError(TypeError):
    """
    非法算符乘积错误

    当检测到直接的算符乘法（如 T * J）时抛出此异常。
    在 VOA 计算中，算符的乘积必须通过 NO(A, B) 显式表示正规序乘积。

    Attributes:
        expr: 导致错误的表达式
        context: 错误发生的上下文（函数名等）
    """

    def __init__(self, expr, context: str, hint: str = ""):
        """
        创建非法算符乘积错误

        Args:
            expr: 导致错误的表达式
            context: 错误发生的上下文（如 "NO(left)", "OPE(right)" 等）
            hint: 额外的提示信息（如何修复）
        """
        msg = (
            f"Illegal operator product detected in {context}: {expr}\n"
            f"Direct multiplication between operator-containing factors is forbidden.\n"
        )
        if hint:
            msg += f"Hint: {hint}\n"
        else:
            msg += (
                "Hint: Use NO(A, B) to represent operator products (normal ordering).\n"
            )

        super().__init__(msg)
        self.expr = expr
        self.context = context
