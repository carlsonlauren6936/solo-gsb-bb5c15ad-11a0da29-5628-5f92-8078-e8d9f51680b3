"""Request/response schemas for the thermal bending API."""

from __future__ import annotations

from typing import List, Literal, Optional

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


class SafeRangeLayerIn(LayerIn):
    """One layer plus its allowable stresses, for the safe-range endpoint."""

    tensile_allowable: float = Field(
        ..., gt=0, allow_inf_nan=False,
        description="Allowable tensile stress, must be finite and > 0 (e.g. Pa).",
    )
    compressive_allowable: float = Field(
        ..., gt=0, allow_inf_nan=False,
        description=(
            "Allowable compressive stress magnitude, must be finite and > 0; "
            "the stress must stay above -compressive_allowable."
        ),
    )


class SafeRangeRequest(BaseModel):
    """Safe temperature range request."""

    layers: List[SafeRangeLayerIn] = Field(
        ..., min_length=2,
        description=(
            "Layers ordered bottom to top, each with allowable stresses; "
            "at least two required."
        ),
    )


class RangeControlOut(BaseModel):
    """Which layer surface and limit type sets a finite interval bound."""

    layer_index: int = Field(..., description="Controlling layer index, 0 = bottom layer.")
    surface: Literal["bottom", "top"] = Field(..., description="Controlling surface.")
    stress_type: Literal["tension", "compression"] = Field(
        ..., description="Limit type reached at the bound."
    )


class SafeRangeResponse(BaseModel):
    """Largest temperature-change interval around zero within stress limits."""

    delta_t_min: Optional[float] = Field(
        ..., description="Lower interval bound; null if unbounded below."
    )
    delta_t_max: Optional[float] = Field(
        ..., description="Upper interval bound; null if unbounded above."
    )
    min_control: Optional[RangeControlOut] = Field(
        ..., description="Controlling surface of the lower bound; null if unbounded."
    )
    max_control: Optional[RangeControlOut] = Field(
        ..., description="Controlling surface of the upper bound; null if unbounded."
    )
