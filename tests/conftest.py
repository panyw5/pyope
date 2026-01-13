"""
Pytest configuration for pyope tests

This file contains shared fixtures and configuration for all tests.
"""

import pytest
from pyope.cache import get_ope_cache
from pyope.registry import ope_registry


@pytest.fixture(autouse=True, scope="function")
def disable_cache_for_tests():
    """
    Disable OPE cache for all tests to ensure test isolation.

    This fixture automatically runs before each test function and:
    1. Disables the global OPE cache
    2. Clears the registry
    3. Re-enables the cache after the test

    This ensures that tests don't interfere with each other through
    cached results.
    """
    cache = get_ope_cache()

    # Disable cache and clear registry before test
    cache.disable()
    cache.clear()
    ope_registry.clear()

    yield

    # Re-enable cache after test (for performance tests if needed)
    cache.enable()
    cache.clear()
    ope_registry.clear()
