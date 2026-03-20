from typing import Any

import sympy as sp

from .operators import Operator

class ConstantOperator(Operator):
    _name: str
    @property
    def is_bosonic(self) -> bool: ...
    @property
    def is_fermionic(self) -> bool: ...
    @property
    def parity(self) -> int: ...

One: ConstantOperator
Zero: ConstantOperator

class DeltaFunction(sp.Function):
    @classmethod
    def eval(cls, i: Any, j: Any) -> Any: ...

def Delta(i: Any, j: Any) -> Any: ...
