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
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional


class SolverError(ValueError):
    """Raised when the system cannot be solved for the given inputs."""


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


@dataclass(frozen=True)
class RatedLayerInput:
    """One material layer with stress allowables.

    thickness:              layer thickness (> 0).
    modulus:                Young's modulus (> 0).
    cte:                    coefficient of thermal expansion (finite).
    tensile_allowable:      maximum tensile stress (> 0, finite).
    compressive_allowable:  maximum compressive stress magnitude (> 0, finite).
    """

    thickness: float
    modulus: float
    cte: float
    tensile_allowable: float
    compressive_allowable: float


@dataclass(frozen=True)
class BoundControl:
    """Which layer surface and stress type sets one interval bound."""

    layer_index: int
    surface: str  # "bottom" or "top"
    stress_type: str  # "tension" or "compression"


@dataclass(frozen=True)
class SafeRangeResult:
    """Largest temperature-change interval around zero within allowables.

    delta_t_min / delta_t_max are None when unbounded (no layer surface
    constrains that side); the matching control is then None as well.
    """

    delta_t_min: Optional[float]
    delta_t_max: Optional[float]
    lower_control: Optional[BoundControl]
    upper_control: Optional[BoundControl]


def safe_temperature_range(layers: List[RatedLayerInput]) -> SafeRangeResult:
    """Largest continuous dT interval around 0 within all stress allowables.

    The laminate is stress-free at dT = 0 and every surface stress is linear
    in dT, so the stress at each layer surface is ``k * dT`` with a constant
    slope ``k``. Each surface therefore imposes

        -compressive_allowable <= k * dT <= tensile_allowable

    which, for ``k != 0``, gives exactly one lower and one upper candidate
    bound; a zero slope constrains nothing. The feasible set is the interval
    between the tightest candidates, which always contains dT = 0.

    Args:
        layers: at least two layers, ordered bottom to top, each with
            finite, strictly positive tensile/compressive allowables.

    Returns:
        SafeRangeResult with the interval bounds (None = unbounded) and,
        for each finite bound, the controlling layer, surface and stress
        type. Ties are broken by layer order, bottom surface first.

    Raises:
        SolverError: if inputs are inconsistent or the system is singular.
    """
    if len(layers) < 2:
        raise SolverError("at least two layers are required")
    for i, layer in enumerate(layers):
        for name, value in (
            ("tensile_allowable", layer.tensile_allowable),
            ("compressive_allowable", layer.compressive_allowable),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise SolverError(
                    f"layer {i}: {name} must be finite and strictly positive"
                )

    # Stresses are homogeneous linear in dT, so solving at dT = 1 yields the
    # per-surface slopes directly.
    unit = solve_free_thermal_bending(
        [LayerInput(l.thickness, l.modulus, l.cte) for l in layers], 1.0
    )

    upper = math.inf
    lower = -math.inf
    upper_control: Optional[BoundControl] = None
    lower_control: Optional[BoundControl] = None

    for i, (layer, res) in enumerate(zip(layers, unit.layers)):
        for surface, z, slope in (
            ("bottom", res.z_bottom, res.stress_bottom),
            ("top", res.z_top, res.stress_top),
        ):
            # A slope is treated as zero only when it is negligible against
            # the magnitude of its own constituents (exact cancellation, e.g.
            # identical layers); genuine slopes always clear this threshold.
            scale = layer.modulus * (
                abs(unit.bottom_strain) + abs(unit.curvature) * z + abs(layer.cte)
            )
            if abs(slope) <= 1e-9 * scale:
                continue
            if slope > 0.0:
                # Heating drives tension, cooling drives compression.
                cand_upper = layer.tensile_allowable / slope
                cand_lower = -layer.compressive_allowable / slope
                type_upper, type_lower = "tension", "compression"
            else:
                cand_upper = -layer.compressive_allowable / slope
                cand_lower = layer.tensile_allowable / slope
                type_upper, type_lower = "compression", "tension"
            # Strict comparisons: the first surface to reach a bound keeps
            # it, i.e. ties go to the lowest layer, bottom surface first.
            if cand_upper < upper:
                upper = cand_upper
                upper_control = BoundControl(i, surface, type_upper)
            if cand_lower > lower:
                lower = cand_lower
                lower_control = BoundControl(i, surface, type_lower)

    return SafeRangeResult(
        delta_t_min=None if lower_control is None else lower,
        delta_t_max=None if upper_control is None else upper,
        lower_control=lower_control,
        upper_control=upper_control,
    )
