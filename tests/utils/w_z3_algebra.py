"""W_{Z3} algebra 测试辅助模块。

测试对象:
- `pyope` 对 `W_{Z3}` 代数的 OPE 恢复、Jacobi 检查和 weight-8 null-state 化零。
- 与 `pyope` 并行的 Wolfram 参考链路，确保抽象代数关系和 free-field realization
    在 `OPEdefs.m` / `OPEdefs.wls` 中也能被一致复现。

测试方法与验证渠道:
- Python 渠道: 直接调用 `pyope` 的 `OPE`、`check_jacobi_identity`、`simplify` 等接口。
- Wolfram 渠道: 通过 `wolframscript` 调用 `OPEdefs.m` / `OPEdefs.wls` 做参考验证。
- 文档渠道: 从 `tests/W_Z3-algebra.md` 提取 strong generator 顺序、OPE、Jacobi
    residual/null-expression 以及选定的 weight-8 null states，作为 ground truth。

测试数据来源:
- `tests/W_Z3-algebra.md`: 用户提供的 `W_{Z3}` 文档，是主要的 OPE 与 null-state 数据源。
- `OPEdefs/OPEdefs.m` 与 `OPEdefs/OPEdefs.wls`: Mathematica / Wolfram 参考实现。
- 本文件中的 ground-truth 常量: 对文档和 Wolfram 输入的结构化转录，供多个 pytest
    文件共享，避免测试里重复声明同一份代数数据。
"""

import itertools
import json
import re
import shutil
import subprocess
import tempfile
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy import Add, Mul

from pyope import (
    BasicOperator,
    Bosonic,
    Fermionic,
    MakeOPE,
    NO,
    OPE,
    One,
    Zero,
    bracket,
    check_jacobi_identity,
    clear_registry,
    d,
    simplify,
    simplify_with_wolfram,
)
from pyope.operators import DerivativeOperator, NormalOrderedOperator


ROOT = Path(__file__).resolve().parents[2]
W_Z3_MD = ROOT / "tests" / "W_Z3-algebra.md"

# 这份 OPE ground truth 是 `tests/W_Z3-algebra.md` 中文档数据的结构化转录。
# `OPEdefs.m` 里 OPE 声明要严格按照算符声明顺序。
# 我们选择的算符声明顺序是 T, J, W, Wbar, G, Gbar, GW, GbarWbar。
# 所以 OPE 用的是 `OPE[W,G] = - GW` 而不是 `OPE[G, W] = GW`。
W_Z3_OPE_GROUND_TRUTH = {
    "T_T": {"4": "-15/2 * One", "2": "2 * T", "1": "Derivative[1][T]"},
    "T_J": {"2": "J", "1": "Derivative[1][J]"},
    "J_J": {"2": "-5 * One"},
    "T_G": {"2": "3/2 * G", "1": "Derivative[1][G]"},
    "T_Gbar": {"2": "3/2 * Gbar", "1": "Derivative[1][Gbar]"},
    "T_W": {"2": "3/2 * W", "1": "Derivative[1][W]"},
    "T_Wbar": {"2": "3/2 * Wbar", "1": "Derivative[1][Wbar]"},
    "T_GW": {"2": "2 * GW", "1": "Derivative[1][GW]"},
    "T_GbarWbar": {"2": "2 * GbarWbar", "1": "Derivative[1][GbarWbar]"},
    "J_G": {"1": "-G"},
    "J_Gbar": {"1": "Gbar"},
    "J_W": {"1": "3 * W"},
    "J_Wbar": {"1": "-3 * Wbar"},
    "J_GW": {"1": "2 * GW"},
    "J_GbarWbar": {"1": "-2 * GbarWbar"},
    "W_Wbar": {
        "3": "-20/9 * One",
        "2": "4/3 * J",
        "1": "2/3 * T - 1/3 * NO[J, J] + 2/3 * Derivative[1][J]",
    },
    "W_G": {"1": "-GW"},
    "Wbar_Gbar": {"1": "-GbarWbar"},
    "G_W": {"1": "GW"},
    "Gbar_W": {},
    "Gbar_Wbar": {"1": "GbarWbar"},
    "G_Wbar": {},
    "G_Gbar": {
        "3": "5 * One",
        "2": "J",
        "1": "-T + 1/2 * Derivative[1][J]",
    },
    "Gbar_GW": {"2": "-3 * W", "1": "-Derivative[1][W]"},
    "G_GbarWbar": {"2": "-3 * Wbar", "1": "-Derivative[1][Wbar]"},
    "Wbar_GW": {
        "2": "4/3 * G",
        "1": "2/3 * NO[J, G] + 2/3 * Derivative[1][G]",
    },
    "W_GbarWbar": {
        "2": "-4/3 * Gbar",
        "1": "2/3 * NO[J, Gbar] - 2/3 * Derivative[1][Gbar]",
    },
    "GW_GbarWbar": {
        "4": "-20/3 * One",
        "3": "8/3 * J",
        "2": "2 * T - 1/3 * NO[J, J] + 4/3 * Derivative[1][J]",
        "1": "2/3 * NO[G, Gbar] - 2/3 * NO[J, T] - 1/3 * NO[J, Derivative[1][J]] + 4/3 * Derivative[1][T] + 1/3 * Derivative[2][J]",
    },
}


# 选取文档中已经人工筛出的 weight-8 null-state 条目，供 Python/Wolfram 双渠道复用。
W_Z3_NULL_GROUND_TRUTH = {
    "basis_ids": [0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
    "special_relations": ["Particular T4 Relation"],
    "reference": "tests/W_Z3-algebra.md",
}


# 记录抽象 Jacobi 检查时的生成元顺序和可直接比较的三元组集合。
_JACOBI_GENERATOR_ORDER = ["T", "J", "W", "Wbar", "G", "Gbar", "GW", "GbarWbar"]


def _directly_supported_jacobi_triples():
    triples = []
    for triple in itertools.combinations_with_replacement(_JACOBI_GENERATOR_ORDER, 3):
        a, b, c = triple
        if (
            f"{a}_{b}" in W_Z3_OPE_GROUND_TRUTH
            and f"{a}_{c}" in W_Z3_OPE_GROUND_TRUTH
            and f"{b}_{c}" in W_Z3_OPE_GROUND_TRUTH
        ):
            triples.append(triple)
    return triples


W_Z3_JACOBI_GROUND_TRUTH = {
    "generator_order": _JACOBI_GENERATOR_ORDER,
    "triples": _directly_supported_jacobi_triples(),
    "all_zero": False,
}


W_Z3_JACOBI_NULL_EXPRESSIONS = [
    "0",
    "2/3 * NO[J, GW] - 2 * NO[W, G] - 4/3 * Derivative[1][GW]",
    "-2/3 * NO[J, GbarWbar] - 2 * NO[Wbar, Gbar] - 4/3 * Derivative[1][GbarWbar]",
    "-2/3 * NO[J, GW] + 2 * NO[W, G] + 4/3 * Derivative[1][GW]",
    "-2/3 * NO[Gbar, GW] - 2/3 * NO[J, Derivative[1][W]] + 2 * NO[T, W] + NO[Derivative[1][J], W] - 2/3 * Derivative[2][W]",
    "2/3 * NO[J, GbarWbar] + 2 * NO[Wbar, Gbar] + 4/3 * Derivative[1][GbarWbar]",
    "8/3 * NO[Gbar, GbarWbar]",
    "-8/3 * NO[G, GW]",
    "2/3 * NO[G, GbarWbar] - 2/3 * NO[J, Derivative[1][Wbar]] - 2 * NO[T, Wbar] + NO[Derivative[1][J], Wbar] + 2/3 * Derivative[2][Wbar]",
    "2/3 * NO[Gbar, GW] + 2/3 * NO[J, Derivative[1][W]] - 2 * NO[T, W] - NO[Derivative[1][J], W] + 2/3 * Derivative[2][W]",
    "8/3 * NO[G, GW]",
    "-2/3 * NO[G, GbarWbar] + 2/3 * NO[J, Derivative[1][Wbar]] + 2 * NO[T, Wbar] - NO[Derivative[1][J], Wbar] - 2/3 * Derivative[2][Wbar]",
    "8/3 * NO[T, GW] + 4/3 * NO[Derivative[1][J], GW] - 8/3 * NO[Derivative[1][W], G] - 4/3 * Derivative[2][GW]",
    "-2/3 * NO[J, Derivative[1][GW]] + 8/3 * NO[T, GW] + 2 * NO[W, Derivative[1][G]] + 2/3 * NO[Derivative[1][J], GW] - 2/3 * NO[Derivative[1][W], G]",
    "-2/3 * NO[J, Derivative[1][GbarWbar]] - 8/3 * NO[T, GbarWbar] - 2 * NO[Wbar, Derivative[1][Gbar]] + 2/3 * NO[Derivative[1][J], GbarWbar] + 2/3 * NO[Derivative[1][Wbar], Gbar]",
    "-8/3 * NO[Gbar, GbarWbar]",
    "-8/3 * NO[T, GbarWbar] + 4/3 * NO[Derivative[1][J], GbarWbar] + 8/3 * NO[Derivative[1][Wbar], Gbar] + 4/3 * Derivative[2][GbarWbar]",
]


W_Z3_JACOBI_NULL_EXPRESSIONS_CANONICAL = [
    "2/3 * NO[J, GW] - 2 * NO[W, G] - 4/3 * Derivative[1][GW]",
    "-2/3 * NO[J, GbarWbar] - 2 * NO[Wbar, Gbar] - 4/3 * Derivative[1][GbarWbar]",
    "-2/3 * NO[Gbar, GW] - 2/3 * NO[J, Derivative[1][W]] + 2 * NO[T, W] + NO[Derivative[1][J], W] - 2/3 * Derivative[2][W]",
    "2/3 * NO[G, GbarWbar] - 2/3 * NO[J, Derivative[1][Wbar]] - 2 * NO[T, Wbar] + NO[Derivative[1][J], Wbar] + 2/3 * Derivative[2][Wbar]",
    "8/3 * NO[G, GW]",
    "8/3 * NO[Gbar, GbarWbar]",
    "8/3 * NO[T, GW] + 4/3 * NO[Derivative[1][J], GW] - 8/3 * NO[Derivative[1][W], G] - 4/3 * Derivative[2][GW]",
    "-2/3 * NO[J, Derivative[1][GbarWbar]] - 8/3 * NO[T, GbarWbar] - 2 * NO[Wbar, Derivative[1][Gbar]] + 2/3 * NO[Derivative[1][J], GbarWbar] + 2/3 * NO[Derivative[1][Wbar], Gbar]",
]


class _WLParser:
    """Parse the limited Wolfram syntax used by the documented ground-truth strings."""

    def __init__(self, text, symbols):
        self.tokens = re.findall(
            r"Derivative|NO|One|[A-Za-z][A-Za-z0-9]*|\d+/\d+|\d+|\[|\]|\(|\)|,|\+|-|\*",
            text,
        )
        self.pos = 0
        self.symbols = symbols

    def peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def pop(self, expected=None):
        token = self.peek()
        if expected is not None and token != expected:
            raise ValueError(f"Expected {expected!r}, got {token!r}")
        self.pos += 1
        return token

    def parse(self):
        expr = self.parse_expr()
        if self.peek() is not None:
            raise ValueError(f"Unexpected trailing token {self.peek()!r}")
        return expr

    def parse_expr(self):
        value = self.parse_term()
        while self.peek() in {"+", "-"}:
            op = self.pop()
            right = self.parse_term()
            value = value + right if op == "+" else value - right
        return value

    def parse_term(self):
        value = self.parse_factor()
        while self.peek() == "*":
            self.pop("*")
            value = value * self.parse_factor()
        return value

    def parse_factor(self):
        token = self.peek()
        if token == "-":
            self.pop("-")
            return -self.parse_factor()
        if token == "(":
            self.pop("(")
            expr = self.parse_expr()
            self.pop(")")
            return expr
        if token == "NO":
            self.pop("NO")
            self.pop("[")
            args = [self.parse_expr()]
            while self.peek() == ",":
                self.pop(",")
                args.append(self.parse_expr())
            self.pop("]")
            return NO(*args)
        if token == "Derivative":
            self.pop("Derivative")
            self.pop("[")
            order = int(self.pop())
            self.pop("]")
            self.pop("[")
            base = self.parse_expr()
            self.pop("]")
            return d(base, order)
        self.pop()
        if token == "One":
            return One
        if re.fullmatch(r"\d+/\d+", token):
            p, q = token.split("/")
            return sp.Rational(int(p), int(q))
        if re.fullmatch(r"\d+", token):
            return sp.Integer(token)
        if token not in self.symbols:
            raise KeyError(f"Unknown symbol {token!r}")
        return self.symbols[token]


def parse_wl_expr(text, symbols):
    """Convert a Wolfram-style reference expression into the SymPy/pyope form used in tests."""

    return _WLParser(text, symbols).parse()


def _safe_op_name(base, prefix):
    if not prefix:
        return base
    sanitized = re.sub(r"[^A-Za-z0-9]", "", prefix)
    return f"{base}{sanitized}"


def make_z3_free_field_data(prefix="wz3ff"):
    """Build the `bc beta-gamma` free-field realization used by the computation tests.

    This is the concrete test object for the OPE recovery and null-state checks in both the
    pure-`pyope` and Wolfram-backed test files.
    """

    clear_registry()
    b = BasicOperator(
        _safe_op_name("b", prefix), fermionic=True, conformal_weight=Fraction(2)
    )
    c = BasicOperator(
        _safe_op_name("c", prefix), fermionic=True, conformal_weight=Fraction(-1)
    )
    beta = BasicOperator(_safe_op_name("beta", prefix), conformal_weight=Fraction(3, 2))
    gamma = BasicOperator(
        _safe_op_name("gamma", prefix), conformal_weight=Fraction(-1, 2)
    )

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    p = 3
    J = p * NO(beta, gamma) + (p - 1) * NO(b, c)
    G = NO(b, gamma)
    Gbar = p * NO(beta, d(c)) + (p - 1) * NO(d(beta), c)
    T = (
        -Fraction(1, 2) * p * NO(beta, d(gamma))
        + (1 - Fraction(1, 2) * p) * NO(d(beta), gamma)
        - Fraction(1, 2) * (p + 1) * NO(b, d(c))
        + Fraction(1, 2) * (1 - p) * NO(d(b), c)
    )
    W = beta
    Wbar = (
        NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
        + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c))))
        - 4 * NO(beta, NO(d(gamma), gamma))
        - Fraction(4, 3) * NO(gamma, NO(b, d(c)))
        + Fraction(2, 3) * NO(gamma, NO(d(b), c))
        + Fraction(2, 3) * NO(d(beta), NO(gamma, gamma))
        - Fraction(8, 3) * NO(d(gamma), NO(b, c))
        + Fraction(10, 9) * d(d(gamma))
    )
    GW = bracket(G, W, 1)
    GbarWbar = bracket(Gbar, Wbar, 1)
    return {
        "b": b,
        "c": c,
        "beta": beta,
        "gamma": gamma,
        "T": T,
        "J": J,
        "W": W,
        "Wbar": Wbar,
        "G": G,
        "Gbar": Gbar,
        "GW": GW,
        "GbarWbar": GbarWbar,
    }


def make_z3_realization_target_data(prefix="wz3real"):
    """Create abstract generators together with their free-field realization map.

    The separate target symbols let tests realize abstract Jacobi residuals into concrete
    free fields without reusing the same operator objects from the source algebra.
    """

    clear_registry()

    b = BasicOperator(
        _safe_op_name("b", prefix), fermionic=True, conformal_weight=Fraction(2)
    )
    c = BasicOperator(
        _safe_op_name("c", prefix), fermionic=True, conformal_weight=Fraction(-1)
    )
    beta = BasicOperator(_safe_op_name("beta", prefix), conformal_weight=Fraction(3, 2))
    gamma = BasicOperator(
        _safe_op_name("gamma", prefix), conformal_weight=Fraction(-1, 2)
    )

    Bosonic(beta, gamma)
    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    p = 3
    Jnew = BasicOperator(_safe_op_name("Jnew", prefix), conformal_weight=1)
    Tnew = BasicOperator(_safe_op_name("Tnew", prefix), conformal_weight=2)
    Wnew = BasicOperator(_safe_op_name("Wnew", prefix), conformal_weight=Fraction(3, 2))
    Wbarnew = BasicOperator(
        _safe_op_name("Wbarnew", prefix), conformal_weight=Fraction(3, 2)
    )
    Gnew = BasicOperator(
        _safe_op_name("Gnew", prefix), fermionic=True, conformal_weight=Fraction(3, 2)
    )
    Gbarnew = BasicOperator(
        _safe_op_name("Gbarnew", prefix),
        fermionic=True,
        conformal_weight=Fraction(3, 2),
    )
    GWnew = BasicOperator(
        _safe_op_name("GWnew", prefix), fermionic=True, conformal_weight=2
    )
    GbarWbarnew = BasicOperator(
        _safe_op_name("GbarWbarnew", prefix), fermionic=True, conformal_weight=2
    )

    Bosonic(Jnew, Tnew, Wnew, Wbarnew)
    Fermionic(Gnew, Gbarnew, GWnew, GbarWbarnew)

    realizations = {
        Jnew: p * NO(beta, gamma) + (p - 1) * NO(b, c),
        Gnew: NO(b, gamma),
        Gbarnew: p * NO(beta, d(c)) + (p - 1) * NO(d(beta), c),
        Tnew: (
            -Fraction(1, 2) * p * NO(beta, d(gamma))
            + (1 - Fraction(1, 2) * p) * NO(d(beta), gamma)
            - Fraction(1, 2) * (p + 1) * NO(b, d(c))
            + Fraction(1, 2) * (1 - p) * NO(d(b), c)
        ),
        Wnew: beta,
        Wbarnew: (
            NO(beta, NO(beta, NO(gamma, NO(gamma, gamma))))
            + 2 * NO(beta, NO(gamma, NO(gamma, NO(b, c))))
            - 4 * NO(beta, NO(d(gamma), gamma))
            - Fraction(4, 3) * NO(gamma, NO(b, d(c)))
            + Fraction(2, 3) * NO(gamma, NO(d(b), c))
            + Fraction(2, 3) * NO(d(beta), NO(gamma, gamma))
            - Fraction(8, 3) * NO(d(gamma), NO(b, c))
            + Fraction(10, 9) * d(d(gamma))
        ),
    }
    realizations[GWnew] = bracket(realizations[Gnew], realizations[Wnew], 1)
    realizations[GbarWbarnew] = bracket(realizations[Gbarnew], realizations[Wbarnew], 1)

    return {
        "b": b,
        "c": c,
        "beta": beta,
        "gamma": gamma,
        "T": Tnew,
        "J": Jnew,
        "W": Wnew,
        "Wbar": Wbarnew,
        "G": Gnew,
        "Gbar": Gbarnew,
        "GW": GWnew,
        "GbarWbar": GbarWbarnew,
        "realizations": realizations,
    }


def _pair_to_expr(left, right, pole_map, ops):
    max_pole = max((int(k) for k in pole_map), default=0)
    data = []
    for pole in range(max_pole, 0, -1):
        expr = parse_wl_expr(pole_map.get(str(pole), "0"), ops)
        data.append(expr)
    OPE[left, right] = MakeOPE(data)


def make_z3_abstract_data(prefix="wz3abs"):
    """Declare the abstract strong generators and load the documented OPE ground truth.

    This is the abstract test object used when the tests intentionally forget the free-field
    realization and only keep the generator list plus the OPE table.
    """

    clear_registry()
    T = BasicOperator(_safe_op_name("T", prefix), conformal_weight=2)
    J = BasicOperator(_safe_op_name("J", prefix), conformal_weight=1)
    W = BasicOperator(_safe_op_name("W", prefix), conformal_weight=Fraction(3, 2))
    Wbar = BasicOperator(_safe_op_name("Wbar", prefix), conformal_weight=Fraction(3, 2))
    G = BasicOperator(
        _safe_op_name("G", prefix), fermionic=True, conformal_weight=Fraction(3, 2)
    )
    Gbar = BasicOperator(
        _safe_op_name("Gbar", prefix), fermionic=True, conformal_weight=Fraction(3, 2)
    )
    GW = BasicOperator(_safe_op_name("GW", prefix), fermionic=True, conformal_weight=2)
    GbarWbar = BasicOperator(
        _safe_op_name("GbarWbar", prefix), fermionic=True, conformal_weight=2
    )
    Bosonic(T, J, W, Wbar)
    Fermionic(G, Gbar, GW, GbarWbar)
    ops = {
        "T": T,
        "J": J,
        "W": W,
        "Wbar": Wbar,
        "G": G,
        "Gbar": Gbar,
        "GW": GW,
        "GbarWbar": GbarWbar,
        "One": One,
    }
    for key, pole_map in W_Z3_OPE_GROUND_TRUTH.items():
        left_name, right_name = key.split("_", 1)
        _pair_to_expr(ops[left_name], ops[right_name], pole_map, ops)
    return ops


def expected_ope_expr_map(ops):
    """Parse the documented OPE table into exact `pyope` expressions for assertion use."""

    parsed = {}
    symbols = {**ops, "One": One}
    for key, pole_map in W_Z3_OPE_GROUND_TRUTH.items():
        parsed[key] = {
            int(pole): simplify(parse_wl_expr(expr, symbols))
            for pole, expr in pole_map.items()
        }
    return parsed


def load_selected_null_relation_sources():
    """Extract the selected weight-8 null-state source text from `tests/W_Z3-algebra.md`."""

    text = W_Z3_MD.read_text(encoding="utf-8")
    results = {}
    for basis_id in W_Z3_NULL_GROUND_TRUTH["basis_ids"]:
        pattern = rf"## Basis {basis_id}.*?```wl\n(.*?)\n```"
        match = re.search(pattern, text, re.S)
        if not match:
            raise ValueError(f"Missing Basis {basis_id} block")
        body = match.group(1).strip()
        body = body.replace("== 0", "").strip()
        if body.startswith("(") and body.endswith(")"):
            body = body[1:-1].strip()
        results[f"Basis {basis_id}"] = body
    special = re.search(r"## Particular T4 Relation.*?```wl\n(.*?)\n```", text, re.S)
    if not special:
        raise ValueError("Missing Particular T4 Relation block")
    body = special.group(1).strip().replace("== 0", "").strip()
    if body.startswith("(") and body.endswith(")"):
        body = body[1:-1].strip()
    results["Particular T4 Relation"] = body
    return results


def build_null_relations(ops):
    """Parse and normalize the selected null relations for Python and Wolfram checks.

    The relation bodies come from the markdown document, then pass through Wolfram-backed
    simplification so downstream tests compare canonicalized expressions instead of raw text.
    """

    symbols = {**ops, "One": One}
    parsed = {
        name: parse_wl_expr(source, symbols)
        for name, source in load_selected_null_relation_sources().items()
    }
    simplified_values = simplify_with_wolfram(list(parsed.values()))
    return {
        name: value
        for name, value in zip(parsed.keys(), simplified_values, strict=False)
    }


def assert_zero_jacobi_matrix(matrix):
    for row in matrix:
        for value in row:
            assert value == Zero


def compute_python_jacobi_summary():
    """Summarize directly supported abstract Jacobi triples that vanish inside `pyope`."""

    ops = make_z3_abstract_data()
    summary = {}
    for triple in W_Z3_JACOBI_GROUND_TRUTH["triples"]:
        matrix = check_jacobi_identity(*(ops[name] for name in triple))
        assert_zero_jacobi_matrix(matrix)
        summary["|".join(triple)] = {
            "rows": len(matrix),
            "cols": len(matrix[0]) if matrix else 0,
            "all_zero": True,
        }
    return summary


def collect_python_jacobi_residuals(ops=None, include_zero=False):
    """Collect nonzero abstract Jacobi residuals for comparison with documented null expressions."""

    if ops is None:
        ops = make_z3_abstract_data()
    residuals = []
    for triple in itertools.combinations_with_replacement(_JACOBI_GENERATOR_ORDER, 3):
        matrix = check_jacobi_identity(*(ops[name] for name in triple))
        for row in matrix:
            for value in row:
                simplified = simplify(value)
                if include_zero or (simplified != 0 and simplified != Zero):
                    residuals.append(simplified)
    return residuals


def expected_jacobi_null_expressions(ops):
    """Parse the full Jacobi null-expression list transcribed from the reference document."""

    symbols = {**ops, "One": One}
    return [
        simplify(parse_wl_expr(expr, symbols)) for expr in W_Z3_JACOBI_NULL_EXPRESSIONS
    ]


def expected_jacobi_null_expressions_canonical(ops):
    """Parse the deduplicated Jacobi null-expression set used for sign-insensitive comparison."""

    symbols = {**ops, "One": One}
    return [
        simplify(parse_wl_expr(expr, symbols))
        for expr in W_Z3_JACOBI_NULL_EXPRESSIONS_CANONICAL
    ]


def canonicalize_up_to_sign(expr):
    simplified = simplify(expr)
    if simplified == 0 or simplified == Zero:
        return "0"
    rep = repr(simplified)
    neg = repr(simplify(-simplified))
    return min(rep, neg)


def realize_exprs_with_free_fields(exprs, abstract_ops=None):
    """Realize abstract expressions in the free-field model before Wolfram/Python simplification."""

    if abstract_ops is None:
        abstract_ops = make_z3_abstract_data("abstractrealize")
    target_ops = make_z3_realization_target_data("realizetarget")
    direct_map = {
        abstract_ops[name]: target_ops["realizations"][target_ops[name]]
        for name in ["T", "J", "W", "Wbar", "G", "Gbar", "GW", "GbarWbar"]
    }
    return [simplify(recursive_realize(expr, direct_map)) for expr in exprs]


def recursive_realize(expr, realization_map):
    replaced = _recursive_realize_once(expr, realization_map)
    prev = None
    current = simplify(replaced)
    for _ in range(10):
        if current == prev:
            break
        prev = current
        current = simplify(_recursive_realize_once(current, realization_map))
    return current


def _recursive_realize_once(expr, realization_map):
    if expr in realization_map:
        return _recursive_realize_once(realization_map[expr], realization_map)
    if isinstance(expr, DerivativeOperator):
        base = _recursive_realize_once(expr.base, realization_map)
        return d(base, expr.order)
    if isinstance(expr, NormalOrderedOperator):
        left = _recursive_realize_once(expr.left, realization_map)
        right = _recursive_realize_once(expr.right, realization_map)
        return NO(left, right)
    if isinstance(expr, Add):
        return Add(
            *(_recursive_realize_once(arg, realization_map) for arg in expr.args)
        )
    if isinstance(expr, Mul):
        return Mul(
            *(_recursive_realize_once(arg, realization_map) for arg in expr.args)
        )
    return expr


def wolframscript_available():
    """Return whether the external Wolfram validation channel is available in the environment."""

    return shutil.which("wolframscript") is not None


def run_wolfram_jacobi_summary():
    """Run an `OPEdefs.m` Jacobi summary script through `wolframscript` and decode its JSON output."""

    triples = W_Z3_JACOBI_GROUND_TRUTH["triples"]
    triple_code = ", ".join(
        "{" + ", ".join(f'"{name}"' for name in triple) + "}" for triple in triples
    )
    script = f'''
Get["{(ROOT / "OPEdefs" / "OPEdefs.m").as_posix()}"];
Bosonic[T, J, W, Wbar];
Fermionic[G, Gbar, GW, GbarWbar];
OPE[T, T] = MakeOPE[{{-15/2 One, 0, 2 T, Derivative[1][T]}}];
OPE[T, J] = MakeOPE[{{0, J, Derivative[1][J]}}];
OPE[J, J] = MakeOPE[{{-5 One, 0}}];
OPE[T, G] = MakeOPE[{{0, 3/2 G, Derivative[1][G]}}];
OPE[T, Gbar] = MakeOPE[{{0, 3/2 Gbar, Derivative[1][Gbar]}}];
OPE[T, W] = MakeOPE[{{0, 3/2 W, Derivative[1][W]}}];
OPE[T, Wbar] = MakeOPE[{{0, 3/2 Wbar, Derivative[1][Wbar]}}];
OPE[T, GW] = MakeOPE[{{0, 2 GW, Derivative[1][GW]}}];
OPE[T, GbarWbar] = MakeOPE[{{0, 2 GbarWbar, Derivative[1][GbarWbar]}}];
OPE[J, G] = MakeOPE[{{-G}}];
OPE[J, Gbar] = MakeOPE[{{Gbar}}];
OPE[J, W] = MakeOPE[{{3 W}}];
OPE[J, Wbar] = MakeOPE[{{-3 Wbar}}];
OPE[J, GW] = MakeOPE[{{2 GW}}];
OPE[J, GbarWbar] = MakeOPE[{{-2 GbarWbar}}];
OPE[W, Wbar] = MakeOPE[{{-20/9 One, 4/3 J, 2/3 T - 1/3 NO[J, J] + 2/3 Derivative[1][J]}}];
OPE[G, W] = MakeOPE[{{GW}}];
OPE[Gbar, W] = MakeOPE[{{0}}];
OPE[Gbar, Wbar] = MakeOPE[{{GbarWbar}}];
OPE[G, Wbar] = MakeOPE[{{0}}];
OPE[G, Gbar] = MakeOPE[{{5 One, J, -T + 1/2 Derivative[1][J]}}];
OPE[Gbar, GW] = MakeOPE[{{-3 W, -Derivative[1][W]}}];
OPE[G, GbarWbar] = MakeOPE[{{-3 Wbar, -Derivative[1][Wbar]}}];
OPE[Wbar, GW] = MakeOPE[{{4/3 G, 2/3 NO[J, G] + 2/3 Derivative[1][G]}}];
OPE[W, GbarWbar] = MakeOPE[{{-4/3 Gbar, 2/3 NO[J, Gbar] - 2/3 Derivative[1][Gbar]}}];
OPE[GW, GbarWbar] = MakeOPE[{{-20/3 One, 8/3 J, 2 T - 1/3 NO[J, J] + 4/3 Derivative[1][J], 2/3 NO[G, Gbar] - 2/3 NO[J, T] - 1/3 NO[J, Derivative[1][J]] + 4/3 Derivative[1][T] + 1/3 Derivative[2][J]}}];
triples = {{{triple_code}}};
summary = Association@Table[
  Module[{{res, key}},
    key = StringRiffle[triple, "|"];
    res = OPEJacobi @@ (ToExpression /@ triple);
    key -> <|"rows" -> Length[res], "cols" -> If[Length[res] == 0, 0, Length[res[[1]]]], "all_zero" -> And @@ Flatten[Map[# === 0 &, res, {{2}}]]|>
  ],
  {{triple, triples}}
];
Print[ExportString[Normal[summary], "JSON"]];
'''
    with tempfile.TemporaryDirectory(prefix="wz3-jacobi-") as tmpdir:
        path = Path(tmpdir) / "jacobi.wls"
        path.write_text(script, encoding="utf-8")
        result = subprocess.run(
            ["wolframscript", "-file", str(path)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout)
    payload_lines = []
    recording = False
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if not recording and (stripped.startswith("[") or stripped.startswith("{")):
            recording = True
        if recording:
            payload_lines.append(line)
    if not payload_lines:
        raise RuntimeError(f"No JSON payload found in output: {result.stdout}")
    return json.loads("\n".join(payload_lines))
