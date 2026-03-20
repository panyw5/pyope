from typing import Any, Callable, List, Optional

def check_jacobi_identity(
    A: Any, B: Any, C: Any, simplify_func: Optional[Callable[..., Any]] = None
) -> List[List[Any]]: ...
def verify_jacobi_identity(
    A: Any, B: Any, C: Any, simplify_func: Optional[Callable[..., Any]] = None
) -> bool: ...
