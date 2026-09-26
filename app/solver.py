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
from typing import List


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
