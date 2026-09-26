"""Composite laminate free thermal bending solver (pure standard library).

Model
-----
Layers are stacked bottom-up; ``z`` is the height measured from the bottom
surface. The in-plane stress of layer ``i`` at height ``z`` is

    sigma_i(z) = E_i * (eps0 + kappa * z - alpha_i * dT)

where ``eps0`` is the strain at the bottom surface (z = 0), ``kappa`` the
curvature of the laminate, ``alpha_i`` the coefficient of thermal expansion
and ``dT`` the uniform temperature change.

Free bending: ``eps0`` and ``kappa`` are chosen so that the resultant axial
force ``N`` and the bending moment ``M`` about the bottom surface both
vanish (resultants per unit width):

    N = sum_i integral(sigma_i dz)     = 0
    M = sum_i integral(sigma_i z dz)   = 0

This gives the 2x2 linear system

    [ A  B ] [ eps0  ]   [ N_T ]
    [ B  D ] [ kappa ] = [ M_T ]

with, per layer (t = thickness, zc = centroid height):

    A   = sum(E * t)
    B   = sum(E * t * zc)
    D   = sum(E * (z_top**3 - z_bot**3) / 3)
    N_T = sum(E * alpha * dT * t)
    M_T = sum(E * alpha * dT * t * zc)

For positive thicknesses and moduli the determinant A*D - B*B is strictly
positive, so the solution always exists and is unique.

Safe temperature range
----------------------
Because ``eps0`` and ``kappa`` are linear in ``dT`` and vanish at ``dT = 0``,
every surface stress is ``slope * dT``. Given per-layer tensile/compressive
allowables, ``safe_temperature_range`` intersects the analytic per-surface
bounds to find the largest ``dT`` interval around zero within all limits.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional


class SolverError(ValueError):
    """Raised when the system cannot be solved for the given inputs."""


# Relative tolerance below which a stress slope (stress per unit dT) is
# treated as zero. Rounding noise in the solved slopes is ~1e-16 of the
# per-surface stress scale, so 1e-9 separates numerical zeros (e.g. a
# homogeneous stack, which is stress-free at any dT) from genuine slopes.
_ZERO_SLOPE_RTOL = 1e-9


@dataclass(frozen=True)
class LayerInput:
    """One material layer.

    thickness: layer thickness (> 0), e.g. meters.
    modulus:   Young's modulus (> 0), e.g. Pa.
    cte:       coefficient of thermal expansion (finite), e.g. 1/K.
    """

    thickness: float
    modulus: float
    cte: float


@dataclass(frozen=True)
class LayerResult:
    """Stress at the bottom and top surfaces of one layer."""

    index: int
    z_bottom: float
    z_top: float
    stress_bottom: float
    stress_top: float


@dataclass(frozen=True)
class SolveResult:
    """Solved bottom strain, curvature and per-layer surface stresses."""

    bottom_strain: float
    curvature: float
    layers: List[LayerResult]


@dataclass(frozen=True)
class StressLimits:
    """Allowable stresses for one layer (magnitudes, must be finite and > 0).

    tensile:     maximum tensile (positive) stress allowed.
    compressive: maximum compressive stress magnitude allowed; the stress
                 must stay above ``-compressive``.
    """

    tensile: float
    compressive: float


@dataclass(frozen=True)
class RangeControl:
    """Which layer surface and limit type sets a finite interval bound.

    layer_index: index of the controlling layer (0 = bottom layer).
    surface:     "bottom" or "top" surface of that layer.
    stress_type: "tension" or "compression".
    """

    layer_index: int
    surface: str
    stress_type: str


@dataclass(frozen=True)
class SafeRangeResult:
    """Largest temperature-change interval around zero within stress limits.

    ``delta_t_min``/``delta_t_max`` are ``None`` when the interval is
    unbounded in that direction (no surface has a non-zero stress slope
    constraining it); the corresponding control is then ``None`` as well.
    """

    delta_t_min: Optional[float]
    delta_t_max: Optional[float]
    min_control: Optional[RangeControl]
    max_control: Optional[RangeControl]


def solve_free_thermal_bending(
    layers: List[LayerInput], delta_t: float
) -> SolveResult:
    """Solve bottom strain and curvature for free thermal bending.

    Args:
        layers: at least two layers, ordered bottom to top.
        delta_t: uniform temperature change applied to the whole laminate.

    Returns:
        SolveResult with bottom strain, curvature and per-layer stresses.

    Raises:
        SolverError: if inputs are inconsistent or the linear system is
            numerically singular/overflowed.
    """
    if len(layers) < 2:
        raise SolverError("at least two layers are required")

    # Interface heights: z_bounds[0] = 0 (bottom), z_bounds[-1] = total height.
    z_bounds = [0.0]
    for layer in layers:
        z_bounds.append(z_bounds[-1] + layer.thickness)

    a = 0.0   # sum(E * t)
    b = 0.0   # sum(E * t * zc)
    d = 0.0   # sum(E * (z_top^3 - z_bot^3) / 3)
    n_t = 0.0  # thermal force resultant  sum(E * alpha * dT * t)
    m_t = 0.0  # thermal moment resultant sum(E * alpha * dT * t * zc)

    for layer, z_bot, z_top in zip(layers, z_bounds[:-1], z_bounds[1:]):
        t = layer.thickness
        e = layer.modulus
        zc = 0.5 * (z_bot + z_top)
        a += e * t
        b += e * t * zc
        d += e * (z_top ** 3 - z_bot ** 3) / 3.0
        n_t += e * layer.cte * delta_t * t
        m_t += e * layer.cte * delta_t * t * zc

    det = a * d - b * b
    if not math.isfinite(det) or det <= 0.0:
        raise SolverError(
            "singular system: cannot solve for bottom strain and curvature "
            "(check that thicknesses and moduli are positive and reasonably scaled)"
        )

    eps0 = (n_t * d - m_t * b) / det
    kappa = (a * m_t - b * n_t) / det

    if not (math.isfinite(eps0) and math.isfinite(kappa)):
        raise SolverError(
            "numerical overflow while solving; rescale inputs to a consistent unit system"
        )

    results: List[LayerResult] = []
    for i, (layer, z_bot, z_top) in enumerate(zip(layers, z_bounds[:-1], z_bounds[1:])):
        free = layer.cte * delta_t
        stress_bottom = layer.modulus * (eps0 + kappa * z_bot - free)
        stress_top = layer.modulus * (eps0 + kappa * z_top - free)
        if not (math.isfinite(stress_bottom) and math.isfinite(stress_top)):
            raise SolverError(
                "numerical overflow while evaluating stresses; "
                "rescale inputs to a consistent unit system"
            )
        results.append(
            LayerResult(
                index=i,
                z_bottom=z_bot,
                z_top=z_top,
                stress_bottom=stress_bottom,
                stress_top=stress_top,
            )
        )

    return SolveResult(bottom_strain=eps0, curvature=kappa, layers=results)


def safe_temperature_range(
    layers: List[LayerInput], limits: List[StressLimits]
) -> SafeRangeResult:
    """Largest dT interval around zero that keeps all stresses within limits.

    Every surface stress is linear in the temperature change and vanishes at
    ``dT = 0``, so ``sigma(dT) = slope * dT`` with ``slope`` the stress solved
    at ``dT = 1``. Each surface with a non-zero slope therefore yields one
    analytic upper and one analytic lower bound:

        slope > 0:  dT <= tensile / slope      (heating -> tension)
                    dT >= -compressive / slope (cooling -> compression)
        slope < 0:  dT <= -compressive / slope (heating -> compression)
                    dT >= tensile / slope      (cooling -> tension)

    A zero-slope surface never limits the interval (its stress stays zero for
    any ``dT``); slopes within rounding noise of zero are treated as zero.
    The safe interval is the intersection of all bounds, i.e.
    ``[max(lowers), min(uppers)]``; it always contains zero because the
    allowables are strictly positive. A direction with no constraining
    surface is unbounded (reported as ``None``).

    Args:
        layers: at least two layers, ordered bottom to top.
        limits: per-layer allowable stresses, same order as ``layers``;
            both magnitudes must be finite and strictly positive.

    Returns:
        SafeRangeResult with the interval bounds and, for each finite bound,
        the controlling layer surface and limit type. Ties on a bound are
        resolved by layer order, bottom surface first.

    Raises:
        SolverError: if inputs are inconsistent or the underlying solve fails.
    """
    if len(layers) != len(limits):
        raise SolverError("number of stress limits must match number of layers")
    for lim in limits:
        if not (math.isfinite(lim.tensile) and lim.tensile > 0.0):
            raise SolverError("tensile allowable must be finite and > 0")
        if not (math.isfinite(lim.compressive) and lim.compressive > 0.0):
            raise SolverError("compressive allowable must be finite and > 0")

    unit = solve_free_thermal_bending(layers, 1.0)
    eps0, kappa = unit.bottom_strain, unit.curvature

    # Candidate bounds: (value, layer_index, surface_rank, control), where
    # surface_rank 0 = bottom, 1 = top (bottom wins ties on the same layer).
    upper = []
    lower = []
    for layer, r, lim in zip(layers, unit.layers, limits):
        surfaces = (
            ("bottom", r.z_bottom, r.stress_bottom),
            ("top", r.z_top, r.stress_top),
        )
        for surface_rank, (surface, z, slope) in enumerate(surfaces):
            # Zero slope (within rounding noise): the surface stress stays
            # zero for any dT, so it never limits the interval.
            noise = layer.modulus * (abs(eps0) + abs(kappa) * z + abs(layer.cte))
            if abs(slope) <= _ZERO_SLOPE_RTOL * noise:
                continue
            if slope > 0.0:
                upper.append((lim.tensile / slope, r.index, surface_rank,
                              RangeControl(r.index, surface, "tension")))
                lower.append((-lim.compressive / slope, r.index, surface_rank,
                              RangeControl(r.index, surface, "compression")))
            else:
                upper.append((-lim.compressive / slope, r.index, surface_rank,
                              RangeControl(r.index, surface, "compression")))
                lower.append((lim.tensile / slope, r.index, surface_rank,
                              RangeControl(r.index, surface, "tension")))

    # Overflow of a bound to infinity means it does not constrain the interval.
    upper = [c for c in upper if math.isfinite(c[0])]
    lower = [c for c in lower if math.isfinite(c[0])]

    if upper:
        best = min(upper, key=lambda c: (c[0], c[1], c[2]))
        delta_t_max, max_control = best[0], best[3]
    else:
        delta_t_max, max_control = None, None

    if lower:
        best = max(lower, key=lambda c: (c[0], -c[1], -c[2]))
        delta_t_min, min_control = best[0], best[3]
    else:
        delta_t_min, min_control = None, None

    return SafeRangeResult(
        delta_t_min=delta_t_min,
        delta_t_max=delta_t_max,
        min_control=min_control,
        max_control=max_control,
    )
