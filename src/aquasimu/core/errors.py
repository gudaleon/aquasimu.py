"""Errors raised by the simulator."""


class MassBalanceError(RuntimeError):
    """Raised when a water, nitrogen or phosphorus balance does not close."""
