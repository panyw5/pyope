from __future__ import annotations

from typing import Any

import sympy as sp

from .api import OPE, bracket
from .constants import One, Zero
from .local_operator import extract_scalar_operator
from .operators import Operator, d
from .simplify import simplify
from .utils import _ascending_pochhammer


def _infer_conformal_weight(expr: Any) -> sp.Expr | None:
    if expr is Zero:
        return None
    if expr is One:
        return sp.Integer(0)

    weight = getattr(expr, "conformal_weight", None)
    if weight is not None:
        return sp.nsimplify(weight)

    if isinstance(expr, sp.Add):
        weights = [_infer_conformal_weight(term) for term in expr.args]
        if any(weight is None for weight in weights):
            return None
        first = weights[0]
        if any(sp.simplify(first - weight) != 0 for weight in weights[1:]):
            raise ValueError("Expression contains inconsistent conformal weights")
        return first

    if isinstance(expr, sp.Mul):
        _, operator = extract_scalar_operator(expr)
        return _infer_conformal_weight(operator)

    if isinstance(expr, Operator):
        weight = getattr(expr, "conformal_weight", None)
        return None if weight is None else sp.nsimplify(weight)

    return None


def quasiprimary_product(
    left: Any,
    right: Any,
    n: int = 0,
    *,
    left_weight: Any | None = None,
    right_weight: Any | None = None,
) -> Any:
    """Compute the quasiprimary completion ``(AB)_n``.

    This implements eqs. (A.7) and (A.8) of Bonetti-Meneghelli-Rastelli,
    i.e. the quasiprimary completion of the OPE coefficient ``{AB}_n``.
    In particular, ``quasiprimary_product(A, B, 0)`` is the paper's
    ``(AB)_0``, which is generally different from the plain normal-ordered
    product ``NO(A, B)``.
    """

    h1 = (
        _infer_conformal_weight(left)
        if left_weight is None
        else sp.nsimplify(left_weight)
    )
    h2 = (
        _infer_conformal_weight(right)
        if right_weight is None
        else sp.nsimplify(right_weight)
    )

    if h1 is None or h2 is None:
        raise ValueError(
            "Could not infer conformal weights; pass left_weight/right_weight"
        )

    ope = OPE(left, right)
    max_shift = max(0, ope.max_pole - n)

    result: Any = Zero
    for p in range(max_shift + 1):
        numerator = (-1) ** p * _ascending_pochhammer(2 * h1 - n - p, p)
        denominator = sp.factorial(p) * _ascending_pochhammer(
            2 * h1 + 2 * h2 - 2 * n - p - 1,
            p,
        )
        coeff = sp.simplify(numerator / denominator)
        term = bracket(left, right, n + p)
        if p > 0:
            term = d(term, p)
        result = result + coeff * term

    return simplify(result)


def qp(
    left: Any,
    right: Any,
    n: int = 0,
    *,
    left_weight: Any | None = None,
    right_weight: Any | None = None,
) -> Any:
    """Short alias for ``quasiprimary_product``."""

    return quasiprimary_product(
        left,
        right,
        n,
        left_weight=left_weight,
        right_weight=right_weight,
    )
