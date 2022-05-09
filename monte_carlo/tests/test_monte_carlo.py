"""
Unit and regression test for the monte_carlo package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import monte_carlo


def test_monte_carlo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "monte_carlo" in sys.modules
