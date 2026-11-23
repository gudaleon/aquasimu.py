"""Reaction and sorption rates used by the reactor kernel.

The numba kernel in `aquasimu.compartments.kernel` evaluates the same
expressions. These functions are the readable form of those rates (per m3 of
pore water, except `langmuir_dqdt`, which is g P per kg media per day).

Nitrogen removal
    first order:  r = k20 * theta^(T - T_ref) * C
    Michaelis-Menten:  r = vmax * theta^(T - T_ref) * C / (K + C)

Phosphorus
    first order:  r = k20 * theta^(T - T_ref) * C   (permanent removal)
    kinetic Langmuir (finite capacity):
        dq/dt = theta^(T - T_ref) * (k_ads * C * (q_max - q) - k_des * q)
        dissolved sink = (rho_bulk / theta_w) * dq/dt
    relaxation to a Langmuir isotherm:
        dq/dt = theta^(T - T_ref) * k_rel * (q_max * C / (K + C) - q)
    relaxation to a linear isotherm:
        dq/dt = theta^(T - T_ref) * k_rel * (K_d * C - q)

Rates use max(C, 0). T_ref is 20 degC in the bioreactor case.
"""
from __future__ import annotations


def temperature_factor(theta: float, temperature: float, t_ref: float = 20.0) -> float:
    """Arrhenius factor theta^(T - T_ref)."""
    return theta ** (temperature - t_ref)


def nitrate_rate(
    concentration: float,
    temperature: float,
    *,
    mode: str,
    k20: float = 0.0,
    vmax: float = 0.0,
    half_saturation: float = 1.0,
    theta: float = 1.0,
    t_ref: float = 20.0,
) -> float:
    """NO3-N removal rate (mass per pore volume per time)."""
    factor = temperature_factor(theta, temperature, t_ref)
    c = concentration if concentration > 0.0 else 0.0
    if mode == "off":
        return 0.0
    if mode == "first_order":
        return k20 * factor * c
    if mode == "mm":
        return vmax * factor * c / (half_saturation + c)
    raise ValueError(f"unknown nitrate kinetics {mode!r}")


def langmuir_dqdt(
    concentration: float,
    sorbed: float,
    temperature: float,
    *,
    k_ads: float,
    k_des: float,
    q_max: float,
    theta: float = 1.0,
    t_ref: float = 20.0,
) -> float:
    """Kinetic Langmuir sorption rate (mass per media mass per time)."""
    factor = temperature_factor(theta, temperature, t_ref)
    c = concentration if concentration > 0.0 else 0.0
    head = q_max - sorbed
    if head < 0.0:
        head = 0.0
    return factor * (k_ads * c * head - k_des * sorbed)


def biofilm_flux_coefficient(kv: float, diffusivity: float, thickness: float, k_l: float) -> float:
    """First-order biofilm flux coefficient J/C (length per time).

    A film of thickness L_f with intrinsic first-order rate k_v and diffusivity
    D_f has Thiele modulus phi = L_f * sqrt(k_v / D_f) and an internal
    coefficient sqrt(k_v * D_f) * tanh(phi). External mass transfer k_L, when
    positive, is in series:

        J/C = 1 / (1/k_L + 1 / (sqrt(k_v D_f) tanh(phi)))

    The reactor multiplies this coefficient by the biofilm area per pore volume
    and by the bulk concentration. The term is first-order only: growth,
    detachment and a moving biofilm boundary are not in this phase.
    """
    # Local import: the kernel is the implementation used by the integrator.
    from aquasimu.compartments.kernel import _biofilm_flux_coeff

    return float(_biofilm_flux_coeff(kv, diffusivity, thickness, k_l))
