"""Wolfram-backed compute helpers for pyope."""

from __future__ import annotations

import ast
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence, Set

import sympy as sp
from sympy import Add, Mul

from .constants import One, Zero
from .local_operator import assert_no_illegal_operator_mul, extract_scalar_operator
from .ope_data import OPEData
from .operators import (
    BasisOperator,
    DerivativeOperator,
    NormalOrderedOperator,
    Operator,
)
from .registry import ope_registry

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WOLFRAM_WRAPPER = _REPO_ROOT / "OPEdefs" / "OPEdefs.wls"


class WolframBackendError(RuntimeError):
    """Raised when the Wolfram backend cannot complete a request."""


def compute_ope(left: Any, right: Any) -> OPEData:
    return _run_operation("OPE", left, right)


def compute_no(left: Any, right: Any) -> Any:
    if not isinstance(left, BasisOperator) or not isinstance(right, BasisOperator):
        from .api import NO
        from .backend import compute_backend

        with compute_backend("sympy"):
            return NO(left, right)
    return _run_operation("NO", left, right)


def simplify_with_wolfram(expr: Any) -> Any:
    return _run_eval(expr)


def _run_operation(operation: str, left: Any, right: Any) -> Any:
    script_path = str(_WOLFRAM_WRAPPER)
    operator_names = _collect_protocol_operator_names((left, right))
    env = {
        "PYOPE_WL_OPERATION": operation,
        "PYOPE_WL_REGISTRY": _encode_registry_state(),
        "PYOPE_WL_LEFT": _encode_expr(left),
        "PYOPE_WL_RIGHT": _encode_expr(right),
        "PYOPE_WL_OPERATORS": ",".join(sorted(operator_names)),
    }
    return _invoke_wolfram(script_path, env, operator_names)


def _run_eval(expr: Any) -> Any:
    script_path = str(_WOLFRAM_WRAPPER)
    operator_names = _collect_protocol_operator_names((expr,))
    env = {
        "PYOPE_WL_OPERATION": "EVAL",
        "PYOPE_WL_REGISTRY": _encode_registry_state(),
        "PYOPE_WL_EXPR": _encode_value(expr),
        "PYOPE_WL_OPERATORS": ",".join(sorted(operator_names)),
    }
    return _invoke_wolfram(script_path, env, operator_names)


def _run_eval_list(exprs: Sequence[Any], operation: str) -> list[Any]:
    script_path = str(_WOLFRAM_WRAPPER)
    expr_list = list(exprs)
    operator_names = _collect_protocol_operator_names(expr_list)
    env = {
        "PYOPE_WL_OPERATION": operation,
        "PYOPE_WL_REGISTRY": _encode_registry_state(),
        "PYOPE_WL_EXPR_LIST": _encode_wolfram_expr_list(expr_list),
        "PYOPE_WL_OPERATORS": ",".join(sorted(operator_names)),
    }
    result = _invoke_wolfram(script_path, env, operator_names)
    if not isinstance(result, list):
        raise WolframBackendError(
            f"Expected list result from Wolfram operation {operation}, got {type(result)!r}"
        )
    return result


def _invoke_wolfram(
    script_path: str, env: Dict[str, str], operator_names: Set[str]
) -> Any:

    result = None
    stdout = ""
    stderr = ""
    for _ in range(2):
        result = subprocess.run(
            ["wolframscript", "-file", script_path],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, **env},
        )
        stdout = result.stdout or ""
        stderr = result.stderr or ""
        if result.returncode == 0:
            break

    if result is None or result.returncode != 0:
        raise WolframBackendError(
            "wolframscript execution failed"
            f" (exit={result.returncode if result is not None else 'unknown'})."
            f" wrapper: {script_path}."
            f" stdout: {stdout.strip() or '<empty>'}"
            f" stderr: {stderr.strip() or '<empty>'}"
        )

    marker = "PYOPE_RESULT:"
    payload = None
    for line in stdout.splitlines():
        if line.startswith(marker):
            payload = line[len(marker) :].strip()

    if payload is None:
        raise WolframBackendError(
            "wolframscript did not emit a parseable result marker. "
            f"stdout: {stdout.strip() or '<empty>'} stderr: {stderr.strip() or '<empty>'}"
        )

    return _decode_expr(payload, operator_names)


def _encode_registry_state() -> str:
    declarations = []
    for operator, parity in ope_registry._parities.items():
        if not isinstance(operator, BasisOperator):
            continue
        decl = "Bosonic" if parity == 0 else "Fermionic"
        declarations.append(f"{decl}[{_encode_symbol_name(operator.name)}]")

    definitions = []
    for (left_name, right_name), ope_data in ope_registry._opes.items():
        encoded_data = _encode_ope_data(ope_data)
        definitions.append(
            f"OPE[{_encode_symbol_name(left_name)}, {_encode_symbol_name(right_name)}] = {encoded_data}"
        )

    return "; ".join(declarations + definitions)


def _encode_ope_data(ope_data: OPEData) -> str:
    poles = []
    max_pole = ope_data.max_pole
    for order in range(max_pole, 0, -1):
        poles.append(_encode_value(ope_data.pole(order)))
    return "MakeOPE[{" + ", ".join(poles) + "}]"


def _encode_value(value: Any) -> str:
    if isinstance(value, list):
        return "{" + ", ".join(_encode_value(item) for item in value) + "}"
    if isinstance(value, tuple):
        return "PyTuple[{" + ", ".join(_encode_value(item) for item in value) + "}]"
    if isinstance(value, Mapping):
        entries = []
        for key, item in value.items():
            if not isinstance(key, str):
                raise WolframBackendError(
                    "Wolfram backend only supports dict containers with string keys"
                )
            entries.append(
                "{" + _encode_string_literal(key) + ", " + _encode_value(item) + "}"
            )
        return "PyDict[{" + ", ".join(entries) + "}]"
    if isinstance(value, OPEData):
        return _encode_ope_data(value)
    if isinstance(value, str):
        return _encode_string_literal(value)
    return _encode_expr(value)


def _encode_expr(expr: Any) -> str:
    if expr == 0 or expr is Zero:
        return "0"
    if expr == 1 or expr is One:
        return "One"
    if isinstance(expr, sp.Float):
        rational = sp.nsimplify(expr)
        if isinstance(rational, (sp.Rational, sp.Integer)):
            return _encode_expr(rational)
        return str(expr)
    if isinstance(expr, sp.Integer):
        return str(int(expr))
    if isinstance(expr, sp.Rational):
        return f"({expr.p}/{expr.q})"
    if isinstance(expr, (int, float)):
        return repr(expr)
    if isinstance(expr, BasisOperator):
        return _encode_symbol_name(expr.name)
    if isinstance(expr, DerivativeOperator):
        return f"Derivative[{expr.order}][{_encode_value(expr.base)}]"
    if isinstance(expr, NormalOrderedOperator):
        return f"NO[{_encode_value(expr.left)}, {_encode_value(expr.right)}]"
    if isinstance(expr, sp.Symbol):
        return _encode_symbol_name(expr.name)
    if isinstance(expr, Add):
        return "(" + " + ".join(_encode_expr(arg) for arg in expr.args) + ")"
    if isinstance(expr, Mul):
        return _encode_mul(expr)

    raise WolframBackendError(f"Unsupported expression for Wolfram backend: {expr!r}")


def _encode_expr_list(exprs: Sequence[Any]) -> str:
    return "[" + ", ".join(_encode_value(expr) for expr in exprs) + "]"


def _encode_wolfram_expr_list(exprs: Sequence[Any]) -> str:
    return "{" + ", ".join(_encode_value(expr) for expr in exprs) + "}"


def chunk_exprs_for_wolfram(
    exprs: Sequence[Any], max_items: int = 32, max_chars: int = 100_000
) -> list[list[Any]]:
    if max_items <= 0:
        raise ValueError("max_items must be positive")
    if max_chars <= 0:
        raise ValueError("max_chars must be positive")

    chunks: list[list[Any]] = []
    current_chunk: list[Any] = []
    current_chars = 2

    for expr in exprs:
        encoded = _encode_value(expr)
        encoded_len = len(encoded)
        if encoded_len + 2 > max_chars:
            raise WolframBackendError(
                "Single expression exceeds Wolfram transport size limit "
                f"({encoded_len} chars > {max_chars})"
            )

        separator_len = 2 if current_chunk else 0
        would_exceed_items = len(current_chunk) >= max_items
        would_exceed_chars = current_chars + separator_len + encoded_len > max_chars
        if current_chunk and (would_exceed_items or would_exceed_chars):
            chunks.append(current_chunk)
            current_chunk = []
            current_chars = 2
            separator_len = 0

        current_chunk.append(expr)
        current_chars += separator_len + encoded_len

    if current_chunk:
        chunks.append(current_chunk)

    return chunks


def _encode_mul(expr: Mul) -> str:
    if not expr.has(Operator):
        return "(" + " * ".join(_encode_value(arg) for arg in expr.args) + ")"

    assert_no_illegal_operator_mul(expr, context="wolfram_backend.encode")
    coeff, op = extract_scalar_operator(expr)

    if op == 1:
        return _encode_value(coeff)
    if coeff == 1:
        return _encode_value(op)
    return f"({_encode_value(coeff)} * {_encode_value(op)})"


def _encode_symbol_name(name: str) -> str:
    return name


def _encode_string_literal(value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return '"' + escaped + '"'


def _decode_expr(payload: str, operator_names: Iterable[str] = ()) -> Any:
    from .api import MakeOPE, NO
    from .backend import compute_backend
    from .operators import dn

    protocol_operator_names: Set[str] = set(operator_names)

    def normalize_protocol_value(value: Any) -> Any:
        if isinstance(value, sp.Symbol) and value.name in protocol_operator_names:
            return BasisOperator(value.name)
        if isinstance(value, NormalOrderedOperator):
            left = normalize_protocol_value(value.left)
            right = normalize_protocol_value(value.right)
            if left is not value.left or right is not value.right:
                return NormalOrderedOperator(left, right)
        if isinstance(value, DerivativeOperator):
            base = normalize_protocol_value(value.base)
            if base is not value.base:
                return DerivativeOperator(base, value.order)
        if isinstance(value, Add):
            return sp.Add(*[normalize_protocol_value(arg) for arg in value.args])
        if isinstance(value, Mul):
            return sp.Mul(*[normalize_protocol_value(arg) for arg in value.args])
        return value

    def decoded_dn(order: Any, value: Any) -> Any:
        return dn(order, normalize_protocol_value(value))

    def decoded_no(left: Any, right: Any) -> Any:
        with compute_backend("sympy"):
            return NO(normalize_protocol_value(left), normalize_protocol_value(right))

    names: Dict[str, Any] = {
        "MakeOPE": MakeOPE,
        "NO": decoded_no,
        "One": One,
        "PyDict": lambda items: dict(items),
        "PyTuple": lambda items: tuple(items),
        "Zero": Zero,
        "dn": decoded_dn,
    }

    for operator in ope_registry._parities.keys():
        if isinstance(operator, BasisOperator):
            names[operator.name] = operator
            protocol_operator_names.add(operator.name)

    for left_name, right_name in ope_registry._opes.keys():
        if left_name not in names:
            names[left_name] = BasisOperator(left_name)
        if right_name not in names:
            names[right_name] = BasisOperator(right_name)
        protocol_operator_names.add(left_name)
        protocol_operator_names.add(right_name)

    for name in protocol_operator_names:
        names.setdefault(name, BasisOperator(name))

    reserved = {"MakeOPE", "NO", "One", "PyDict", "PyTuple", "Zero", "dn"}
    for token in set(re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", payload)):
        if token not in names and token not in reserved:
            names[token] = sp.Symbol(token)

    try:
        node = ast.parse(payload, mode="eval")
        return _eval_protocol_ast(node, names, payload)
    except WolframBackendError:
        raise
    except Exception as exc:  # pragma: no cover - defensive
        raise WolframBackendError(
            f"Failed to decode Wolfram payload: {payload}"
        ) from exc


def _collect_protocol_operator_names(exprs: Iterable[Any]) -> Set[str]:
    names: Set[str] = set()

    def visit(expr: Any) -> None:
        if isinstance(expr, list | tuple):
            for item in expr:
                visit(item)
            return
        if isinstance(expr, Mapping):
            for item in expr.values():
                visit(item)
            return
        if isinstance(expr, BasisOperator):
            names.add(expr.name)
            return
        if isinstance(expr, DerivativeOperator):
            visit(expr.base)
            return
        if isinstance(expr, NormalOrderedOperator):
            visit(expr.left)
            visit(expr.right)
            return
        if isinstance(expr, sp.Expr):
            for operator in expr.atoms(Operator):
                if isinstance(operator, BasisOperator):
                    names.add(operator.name)
                elif isinstance(operator, DerivativeOperator):
                    visit(operator.base)
                elif isinstance(operator, NormalOrderedOperator):
                    visit(operator.left)
                    visit(operator.right)

    for expr in exprs:
        visit(expr)

    for operator in ope_registry._parities.keys():
        if isinstance(operator, BasisOperator):
            names.add(operator.name)

    for left_name, right_name in ope_registry._opes.keys():
        names.add(left_name)
        names.add(right_name)

    return names


def _eval_protocol_ast(node: ast.AST, names: Dict[str, Any], payload: str) -> Any:
    if isinstance(node, ast.Expression):
        return _eval_protocol_ast(node.body, names, payload)

    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in names:
            raise WolframBackendError(
                f"Unsupported Wolfram function in payload: {payload}"
            )
        if node.keywords:
            raise WolframBackendError(
                f"Keyword arguments are not allowed in payload: {payload}"
            )
        func = names[node.func.id]
        args = [_eval_protocol_ast(arg, names, payload) for arg in node.args]
        return func(*args)

    if isinstance(node, ast.Name):
        if node.id not in names:
            raise WolframBackendError(f"Unknown token in Wolfram payload: {node.id}")
        return names[node.id]

    if isinstance(node, ast.List):
        return [_eval_protocol_ast(elt, names, payload) for elt in node.elts]

    if isinstance(node, ast.Tuple):
        return tuple(_eval_protocol_ast(elt, names, payload) for elt in node.elts)

    if isinstance(node, ast.Dict):
        keys = [_eval_protocol_ast(key, names, payload) for key in node.keys]
        values = [_eval_protocol_ast(value, names, payload) for value in node.values]
        return dict(zip(keys, values, strict=True))

    if isinstance(node, ast.Constant):
        if isinstance(node.value, int):
            return node.value
        if isinstance(node.value, float):
            rational = sp.nsimplify(node.value)
            if isinstance(rational, (sp.Integer, sp.Rational)):
                return rational
            return node.value
        if isinstance(node.value, str):
            return node.value
        raise WolframBackendError(f"Unsupported constant in Wolfram payload: {payload}")

    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_eval_protocol_ast(node.operand, names, payload)

    if isinstance(node, ast.BinOp):
        left = _eval_protocol_ast(node.left, names, payload)
        right = _eval_protocol_ast(node.right, names, payload)

        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Div):
            if isinstance(left, int) and isinstance(right, int):
                return sp.Rational(left, right)
            return left / right
        if isinstance(node.op, ast.Mult):
            if _is_operator_like(left) and _is_operator_like(right):
                raise WolframBackendError(
                    f"Illegal operator multiplication in Wolfram payload: {payload}"
                )
            product = left * right
            if isinstance(product, sp.Expr):
                assert_no_illegal_operator_mul(
                    product, context="wolfram_backend.decode"
                )
            return product

    raise WolframBackendError(f"Unsupported Wolfram payload syntax: {payload}")


def _is_operator_like(value: Any) -> bool:
    return isinstance(value, Operator) or (
        isinstance(value, sp.Expr) and value.has(Operator)
    )
