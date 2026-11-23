"""Simulation checks. Mass-balance closure is mandatory."""

from aquasimu.checks.massbalance import (
    assert_unit_balance,
    component_balance,
    train_balance,
)
from aquasimu.core.errors import MassBalanceError

__all__ = [
    "MassBalanceError",
    "assert_unit_balance",
    "component_balance",
    "train_balance",
]
