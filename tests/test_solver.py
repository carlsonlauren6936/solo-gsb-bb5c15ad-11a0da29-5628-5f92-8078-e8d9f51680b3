"""Sanity checks for the thermal bending solver (pure stdlib, run with python3)."""

import math
import sys

sys.path.insert(0, ".")

from app.solver import (
    LayerInput,
    SolverError,
    StressLimits,
    safe_temperature_range,
    solve_free_thermal_bending,
)


def check(name, condition):
    print(("PASS" if condition else "FAIL"), "-", name)
    if not condition:
        sys.exit(1)


def resultants(result, layers, delta_t):
    """Recompute N and M from the returned stresses to verify equilibrium."""
    n = m = 0.0
    for r, layer in zip(result.layers, layers):
        t = r.z_top - r.z_bottom
        zc = 0.5 * (r.z_bottom + r.z_top)
        avg = 0.5 * (r.stress_bottom + r.stress_top)  # stress is linear in z
        n += avg * t
        m += avg * t * zc + (r.stress_top - r.stress_bottom) * t ** 2 / 12.0
    return n, m


# 1. Two identical layers -> free expansion, zero curvature, zero stress.
layers = [LayerInput(1.0, 200e9, 12e-6), LayerInput(1.0, 200e9, 12e-6)]
res = solve_free_thermal_bending(layers, 100.0)
check("identical layers: eps0 = alpha*dT", math.isclose(res.bottom_strain, 12e-6 * 100, rel_tol=1e-9))
check("identical layers: kappa = 0", math.isclose(res.curvature, 0.0, abs_tol=1e-20))
check(
    "identical layers: stress-free",
    all(math.isclose(r.stress_bottom, 0.0, abs_tol=1e-6) and math.isclose(r.stress_top, 0.0, abs_tol=1e-6) for r in res.layers),
)

# 2. Bimetal, equal E and t: kappa = 3*(a2-a1)*dT / (4*t) (Timoshenko).
t, e, dt = 1.0, 100e9, 50.0
a1, a2 = 10e-6, 20e-6
layers = [LayerInput(t, e, a1), LayerInput(t, e, a2)]
res = solve_free_thermal_bending(layers, dt)
expected_kappa = 3.0 * (a2 - a1) * dt / (4.0 * t)
check("bimetal: curvature matches Timoshenko", math.isclose(res.curvature, expected_kappa, rel_tol=1e-9))

# 3. Equilibrium residuals are zero for a general 4-layer stack.
layers = [
    LayerInput(0.5e-3, 210e9, 12e-6),
    LayerInput(1.2e-3, 70e9, 23e-6),
    LayerInput(0.3e-3, 110e9, 17e-6),
    LayerInput(2.0e-3, 3e9, 80e-6),
]
res = solve_free_thermal_bending(layers, -80.0)
n, m = resultants(res, layers, -80.0)
scale = sum(l.modulus * l.thickness for l in layers)
check("general stack: N = 0", abs(n) < 1e-6 * scale)
check("general stack: M = 0", abs(m) < 1e-6 * scale * 1e-3)

# 4. Stress continuity at interfaces (stress is linear within each layer).
for prev, cur in zip(res.layers, res.layers[1:]):
    check("interface heights match", math.isclose(prev.z_top, cur.z_bottom, rel_tol=1e-12))

# 5. Invalid inputs raise SolverError.
for bad in (
    [LayerInput(1.0, 1.0, 1.0)],  # single layer
):
    try:
        solve_free_thermal_bending(bad, 1.0)
        check("invalid input rejected", False)
    except SolverError:
        check("single layer rejected", True)

# 6. Safe range: identical layers -> stress-free at any dT -> unbounded.
layers = [LayerInput(1.0, 200e9, 12e-6), LayerInput(1.0, 200e9, 12e-6)]
limits = [StressLimits(1e8, 1e8), StressLimits(1e8, 1e8)]
res = safe_temperature_range(layers, limits)
check("safe range: identical layers unbounded below", res.delta_t_min is None)
check("safe range: identical layers unbounded above", res.delta_t_max is None)
check("safe range: no controls when unbounded", res.min_control is None and res.max_control is None)

# 7. Safe range: bimetal, equal E/t/limits -> analytic bounds.
#    Slopes at dT=1: l0.bottom -2.5e5, l0.top +5e5, l1.bottom -5e5, l1.top +2.5e5.
#    With T=C=1e8 the bounds are +-200 (tension/compression per surface).
t, e = 1.0, 100e9
layers = [LayerInput(t, e, 10e-6), LayerInput(t, e, 20e-6)]
limits = [StressLimits(1e8, 1e8), StressLimits(1e8, 1e8)]
res = safe_temperature_range(layers, limits)
check("safe range: bimetal upper bound", math.isclose(res.delta_t_max, 200.0, rel_tol=1e-9))
check("safe range: bimetal lower bound", math.isclose(res.delta_t_min, -200.0, rel_tol=1e-9))

# 8. Safe range: stresses at both bounds stay within limits, and the
#    controlling surface sits exactly at its allowable.
for dT_bound, control in ((res.delta_t_min, res.min_control), (res.delta_t_max, res.max_control)):
    at_bound = solve_free_thermal_bending(layers, dT_bound)
    for r in at_bound.layers:
        lim = limits[r.index]
        check(
            "safe range: bound feasible (layer %d)" % r.index,
            r.stress_bottom <= lim.tensile * (1 + 1e-9)
            and r.stress_bottom >= -lim.compressive * (1 + 1e-9)
            and r.stress_top <= lim.tensile * (1 + 1e-9)
            and r.stress_top >= -lim.compressive * (1 + 1e-9),
        )
    surf = at_bound.layers[control.layer_index]
    stress = surf.stress_bottom if control.surface == "bottom" else surf.stress_top
    lim = limits[control.layer_index]
    expected = lim.tensile if control.stress_type == "tension" else -lim.compressive
    check("safe range: controlling stress at allowable", math.isclose(stress, expected, rel_tol=1e-9))

# 9. Safe range: slightly outside a bound some limit is exceeded.
outside = solve_free_thermal_bending(layers, res.delta_t_max * 1.000001)
violated = any(
    r.stress_bottom > limits[r.index].tensile
    or r.stress_top > limits[r.index].tensile
    or r.stress_bottom < -limits[r.index].compressive
    or r.stress_top < -limits[r.index].compressive
    for r in outside.layers
)
check("safe range: outside upper bound violates", violated)

# 10. Safe range tie-break: symmetric A/B/A stack (kappa = 0, all A-surface
#     slopes equal) -> tied bounds resolve to layer 0, bottom surface.
layers = [
    LayerInput(1.0, 100.0, 2.0),
    LayerInput(1.0, 100.0, 1.0),
    LayerInput(1.0, 100.0, 2.0),
]
limits = [StressLimits(10.0, 20.0), StressLimits(100.0, 100.0), StressLimits(10.0, 20.0)]
res = safe_temperature_range(layers, limits)
# A-surface slope = 100*(5/3 - 2) = -100/3 -> upper 3*20/100 = 0.6 (compression),
# lower -3*10/100 = -0.3 (tension); loose B limits never control.
check("safe range: symmetric upper bound", math.isclose(res.delta_t_max, 0.6, rel_tol=1e-9))
check("safe range: symmetric lower bound", math.isclose(res.delta_t_min, -0.3, rel_tol=1e-9))
check(
    "safe range: tie-break layer 0 bottom (upper)",
    res.max_control.layer_index == 0
    and res.max_control.surface == "bottom"
    and res.max_control.stress_type == "compression",
)
check(
    "safe range: tie-break layer 0 bottom (lower)",
    res.min_control.layer_index == 0
    and res.min_control.surface == "bottom"
    and res.min_control.stress_type == "tension",
)

# 11. Safe range: invalid allowables raise SolverError.
for bad_limits in (
    [StressLimits(0.0, 1.0), StressLimits(1.0, 1.0)],    # tensile not > 0
    [StressLimits(1.0, -1.0), StressLimits(1.0, 1.0)],   # compressive not > 0
    [StressLimits(math.inf, 1.0), StressLimits(1.0, 1.0)],  # not finite
    [StressLimits(1.0, math.nan), StressLimits(1.0, 1.0)],  # not finite
    [StressLimits(1.0, 1.0)],                            # wrong count
):
    try:
        safe_temperature_range([LayerInput(1.0, 1.0, 1.0), LayerInput(1.0, 1.0, 2.0)], bad_limits)
        check("invalid allowables rejected", False)
    except SolverError:
        check("invalid allowables rejected", True)

print("All solver checks passed.")
