"""FastAPI application for composite laminate free thermal bending."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from .schemas import LayerOut, SolveRequest, SolveResponse
from .solver import LayerInput, SolverError, solve_free_thermal_bending

app = FastAPI(
    title="Composite Laminate Thermal Bending API",
    version="1.0.0",
    description=(
        "Solves free thermal bending of a composite laminate: given layers "
        "(bottom to top) with thickness, Young's modulus and CTE, plus a "
        "uniform temperature change, returns the bottom strain, curvature "
        "and the stress at the bottom/top surface of every layer."
    ),
)


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(
    request: Request, exc: RequestValidationError
) -> JSONResponse:
    """Return a clear, flat list of parameter errors (HTTP 422)."""
    errors = []
    for err in exc.errors():
        location = ".".join(str(part) for part in err["loc"] if part != "body")
        errors.append({"field": location or "body", "message": err["msg"]})
    return JSONResponse(status_code=422, content={"detail": errors})


@app.exception_handler(SolverError)
async def solver_exception_handler(request: Request, exc: SolverError) -> JSONResponse:
    """Map numerical/physical inconsistencies to a clear HTTP 422."""
    return JSONResponse(status_code=422, content={"detail": str(exc)})


@app.get("/health", summary="Health check")
def health() -> dict:
    return {"status": "ok"}


@app.post("/solve", response_model=SolveResponse, summary="Solve free thermal bending")
def solve(request: SolveRequest) -> SolveResponse:
    layers = [
        LayerInput(thickness=l.thickness, modulus=l.modulus, cte=l.cte)
        for l in request.layers
    ]
    result = solve_free_thermal_bending(layers, request.delta_t)
    return SolveResponse(
        bottom_strain=result.bottom_strain,
        curvature=result.curvature,
        layers=[
            LayerOut(
                index=r.index,
                z_bottom=r.z_bottom,
                z_top=r.z_top,
                stress_bottom=r.stress_bottom,
                stress_top=r.stress_top,
            )
            for r in result.layers
        ],
    )
