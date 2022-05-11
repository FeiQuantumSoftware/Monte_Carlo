"""
Unit and regression test for the monte_carlo package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import monte_carlo

import numpy as np


def test_spinconfig():
    test_spin = monte_carlo.SpinConfig(8)

    """basic properties"""
    expected_N_length = 8
    calculated_N_length = test_spin.N_length

    expected_iMax = 256
    calculated_iMax = test_spin.iMax

    # initialization
    expected_init_input_decimal = [0, 0, 0, 0, 1, 0, 1, 0]
    calculated_init_input_decimal = test_spin.init_input_decimal(10)

    # property
    expected_magnetization = -4
    calculated_magnetization = test_spin.magnetization()

    expected_hamiltonian = -4.4
    calculated_hamiltonian = test_spin.hamiltonian()

    # observables
    expected_observable_theory = (-3.6772068591549063, -0.5894627003462397, 0.32593415709340545, 0.5351140013397603)
    test_spin_observable_theory = test_spin.observable_theory()

    # translation ability
    expected_input_str = [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1]
    calculated_input_str = test_spin.input_str("++-+---+--+")


#     """assert"""
    assert expected_N_length == calculated_N_length
    assert expected_iMax == calculated_iMax
    assert expected_init_input_decimal == calculated_init_input_decimal

    assert expected_magnetization == calculated_magnetization
    assert expected_hamiltonian == calculated_hamiltonian

    assert expected_observable_theory == test_spin_observable_theory

    assert expected_input_str == calculated_input_str
