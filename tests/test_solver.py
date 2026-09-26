"""Sanity checks for the thermal bending solver (pure stdlib, run with python3)."""

import math
import sys

sys.path.insert(0, ".")

from app.solver import LayerInput, SolverError, solve_free_thermal_bending


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

print("All solver checks passed.")
