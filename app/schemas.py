"""Request/response schemas for the thermal bending API."""

from __future__ import annotations

from typing import List

from pydantic import BaseModel, Field


class LayerIn(BaseModel):
    """One material layer, ordered bottom to top in the request."""

    thickness: float = Field(
        ..., gt=0, allow_inf_nan=False,
        description="Layer thickness, must be > 0 (e.g. meters).",
    )
    modulus: float = Field(
        ..., gt=0, allow_inf_nan=False,
        description="Young's modulus, must be > 0 (e.g. Pa).",
    )
    cte: float = Field(
        ..., allow_inf_nan=False,
        description="Coefficient of thermal expansion, must be finite (e.g. 1/K).",
    )


class SolveRequest(BaseModel):
    """Free thermal bending request."""

    delta_t: float = Field(
        ..., allow_inf_nan=False,
        description="Uniform temperature change applied to the whole laminate (e.g. K).",
    )
    layers: List[LayerIn] = Field(
        ..., min_length=2,
        description="Layers ordered bottom to top; at least two required.",
    )


class LayerOut(BaseModel):
    """Stress at the bottom and top surfaces of one layer."""

    index: int = Field(..., description="Layer index, 0 = bottom layer.")
    z_bottom: float = Field(..., description="Height of the layer bottom surface.")
    z_top: float = Field(..., description="Height of the layer top surface.")
    stress_bottom: float = Field(..., description="Stress at the layer bottom surface.")
    stress_top: float = Field(..., description="Stress at the layer top surface.")


class SolveResponse(BaseModel):
    """Solved bottom strain, curvature and per-layer surface stresses."""

    bottom_strain: float = Field(..., description="Strain at the bottom surface (z = 0).")
    curvature: float = Field(..., description="Curvature of the laminate (1/length).")
    layers: List[LayerOut]
