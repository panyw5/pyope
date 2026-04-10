import warnings

from pyope import BasicOperator, Fermionic, MakeOPE, One, OPE, clear_registry
from pyope.registry import ope_registry


def setup_function():
    clear_registry()


def teardown_function():
    clear_registry()


def test_clear_registry_removes_registered_operators_and_opes():
    b = BasicOperator("b", fermionic=True)
    c = BasicOperator("c", fermionic=True)

    Fermionic(b, c)
    OPE[b, c] = MakeOPE([One])

    assert ope_registry.is_registered(b)
    assert ope_registry.is_registered(c)
    assert ope_registry.has_ope(b, c)

    clear_registry()

    assert not ope_registry.is_registered(b)
    assert not ope_registry.is_registered(c)
    assert not ope_registry.has_ope(b, c)


def test_clear_registry_allows_reregister_without_warning():
    b = BasicOperator("b", fermionic=True)

    Fermionic(b)
    clear_registry()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        Fermionic(b)

    assert caught == []
