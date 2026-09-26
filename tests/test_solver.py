"""Sanity checks for the thermal bending solver (pure stdlib, run with python3)."""

import math
import sys

sys.path.insert(0, ".")

from app.solver import (
    LayerInput,
    RatedLayerInput,
    SolverError,
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

# 6. Safe temperature range: bounds match the per-surface linear relations.
rated = [
    RatedLayerInput(0.001, 2.10e11, 1.2e-5, 100e6, 50e6),
    RatedLayerInput(0.002, 7.0e10, 2.3e-5, 100e6, 50e6),
]
rng = safe_temperature_range(rated)
unit = solve_free_thermal_bending(
    [LayerInput(l.thickness, l.modulus, l.cte) for l in rated], 1.0
)
candidates_upper, candidates_lower = [], []
for layer, res in zip(rated, unit.layers):
    for slope in (res.stress_bottom, res.stress_top):
        if slope > 0:
            candidates_upper.append(layer.tensile_allowable / slope)
            candidates_lower.append(-layer.compressive_allowable / slope)
        elif slope < 0:
            candidates_upper.append(-layer.compressive_allowable / slope)
            candidates_lower.append(layer.tensile_allowable / slope)
check("range: upper bound", math.isclose(rng.delta_t_max, min(candidates_upper), rel_tol=1e-12))
check("range: lower bound", math.isclose(rng.delta_t_min, max(candidates_lower), rel_tol=1e-12))
check("range: interval contains zero", rng.delta_t_min < 0.0 < rng.delta_t_max)

# At each bound the controlling surface reaches exactly its allowable, and
# every surface stays within its allowables inside the interval.
for bound, control in ((rng.delta_t_min, rng.lower_control), (rng.delta_t_max, rng.upper_control)):
    res = solve_free_thermal_bending(
        [LayerInput(l.thickness, l.modulus, l.cte) for l in rated], bound
    )
    ctrl_layer, ctrl_res = rated[control.layer_index], res.layers[control.layer_index]
    ctrl_stress = ctrl_res.stress_bottom if control.surface == "bottom" else ctrl_res.stress_top
    expected = (
        ctrl_layer.tensile_allowable
        if control.stress_type == "tension"
        else -ctrl_layer.compressive_allowable
    )
    check(f"range: {control.stress_type} allowable met at bound", math.isclose(ctrl_stress, expected, rel_tol=1e-9))
    for layer, r in zip(rated, res.layers):
        for stress in (r.stress_bottom, r.stress_top):
            check(
                "range: stresses within allowables",
                stress <= layer.tensile_allowable * (1 + 1e-9)
                and stress >= -layer.compressive_allowable * (1 + 1e-9),
            )

# 7. Unbounded when no surface stress depends on dT (identical layers).
same = [RatedLayerInput(1.0, 200e9, 12e-6, 100e6, 50e6)] * 2
rng = safe_temperature_range(same)
check("range: unbounded below", rng.delta_t_min is None and rng.lower_control is None)
check("range: unbounded above", rng.delta_t_max is None and rng.upper_control is None)

# 8. Ties on a bound go to the lowest layer, bottom surface first.
#    Symmetric [A, B, A] stack with power-of-two data: curvature is exactly
#    zero, so both surfaces of both A layers share one exact slope. The B
#    layer gets 4x allowables so the tied A surfaces set both bounds.
sig = 100e6
rated = [
    RatedLayerInput(1.0, 2.0 ** 50, 2.0 ** -20, sig, sig),
    RatedLayerInput(1.0, 2.0 ** 52, 2.0 ** -18, 4 * sig, 4 * sig),
    RatedLayerInput(1.0, 2.0 ** 50, 2.0 ** -20, sig, sig),
]
rng = safe_temperature_range(rated)
unit = solve_free_thermal_bending(
    [LayerInput(l.thickness, l.modulus, l.cte) for l in rated], 1.0
)
check("range: symmetric stack has zero curvature", unit.curvature == 0.0)
k_a = unit.layers[0].stress_bottom  # slope shared by all A-layer surfaces
k_b = unit.layers[1].stress_bottom
check("range: A layers tied", k_a == unit.layers[0].stress_top == unit.layers[2].stress_bottom == unit.layers[2].stress_top)
check("range: B layer not controlling", abs(k_b) > abs(k_a))
check("range: tie lower control", (rng.lower_control.layer_index, rng.lower_control.surface, rng.lower_control.stress_type) == (0, "bottom", "compression"))
check("range: tie upper control", (rng.upper_control.layer_index, rng.upper_control.surface, rng.upper_control.stress_type) == (0, "bottom", "tension"))
check("range: tie lower bound", rng.delta_t_min == -sig / k_a)
check("range: tie upper bound", rng.delta_t_max == sig / k_a)

# 9. Non-positive or non-finite allowables raise SolverError.
for bad in (
    [RatedLayerInput(1.0, 1.0, 1.0, 0.0, 1.0), RatedLayerInput(1.0, 1.0, 2.0, 1.0, 1.0)],
    [RatedLayerInput(1.0, 1.0, 1.0, 1.0, -5.0), RatedLayerInput(1.0, 1.0, 2.0, 1.0, 1.0)],
    [RatedLayerInput(1.0, 1.0, 1.0, math.inf, 1.0), RatedLayerInput(1.0, 1.0, 2.0, 1.0, 1.0)],
    [RatedLayerInput(1.0, 1.0, 1.0, 1.0, math.nan), RatedLayerInput(1.0, 1.0, 2.0, 1.0, 1.0)],
    [RatedLayerInput(1.0, 1.0, 1.0, 1.0, 1.0)],  # single layer
):
    try:
        safe_temperature_range(bad)
        check("invalid allowables rejected", False)
    except SolverError:
        check("invalid allowables rejected", True)

print("All solver checks passed.")
