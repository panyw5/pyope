from __future__ import annotations

from typing import Any, Iterable

import sympy as sp

from .constants import Zero
from .operators import d
from .simplify import simplify
from .utils import _ascending_pochhammer


def compact_family_poles(
    left_weight: Any,
    right_weight: Any,
    families: Iterable[tuple[Any, Any, Any]],
) -> dict[int, Any]:
    """Recover explicit OPE poles from compact sl(2)-covariant notation.

    Each family is a tuple ``(coefficient, quasiprimary, family_weight)``
    representing a term in the compact notation of appendix A, i.e.

    O1(z1) O2(z2) ~ coeff * z12^(-(h1+h2-h)) D_{h1,h2;h}(z12, d/dz2) quasiprimary(z2)

    with the differential operator ``D`` from eq. (A.6).

    Returns a dictionary mapping pole order ``q >= 1`` to the explicit local
    operator multiplying ``z12^{-q}``.
    """

    h1 = sp.nsimplify(left_weight)
    h2 = sp.nsimplify(right_weight)
    poles: dict[int, Any] = {}

    for coefficient, quasiprimary, family_weight in families:
        coeff = sp.nsimplify(coefficient)
        h = sp.nsimplify(family_weight)
        highest_pole = sp.nsimplify(h1 + h2 - h)

        if highest_pole < 1:
            continue

        if h == 0:
            if int(highest_pole) != highest_pole:
                raise ValueError("Identity family highest pole must be integral")
            pole_order = int(highest_pole)
            poles[pole_order] = simplify(
                poles.get(pole_order, Zero) + coeff * quasiprimary
            )
            continue

        if int(highest_pole) != highest_pole:
            raise ValueError(
                "Highest pole must be integral in compact family expansion"
            )

        for k in range(int(highest_pole)):
            pole_order = int(highest_pole) - k
            descendant_coeff = sp.simplify(
                _ascending_pochhammer(h + h1 - h2, k)
                / (sp.factorial(k) * _ascending_pochhammer(2 * h, k))
            )
            descendant = quasiprimary if k == 0 else d(quasiprimary, k)
            poles[pole_order] = simplify(
                poles.get(pole_order, Zero) + coeff * descendant_coeff * descendant
            )

    return {pole: simplify(expr) for pole, expr in poles.items()}
