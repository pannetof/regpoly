"""Pydantic models for request/response validation."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


# ── ParamSpec and family metadata ────────────────────────────────────────

class ParamSpecModel(BaseModel):
    name: str
    type: str
    structural: bool
    has_default: bool
    default: Any = None
    rand_type: str = ""
    rand_args: str = ""
    optimizable: bool = False


class FamilyInfo(BaseModel):
    name: str                       # C++ class name (canonical)
    aliases: list[str] = []
    params: list[ParamSpecModel]


class TransformationInfo(BaseModel):
    name: str                       # C++ type string
    params: list[ParamSpecModel]


# ── Primitive search ─────────────────────────────────────────────────────

class PrimitiveSearchCreate(BaseModel):
    family: str
    # L is no longer user-facing; the server sets L = k.  Kept for
    # backward compatibility with API clients that still send it.
    L: int = 0
    structural_params: dict[str, Any]
    fixed_params: dict[str, Any] = Field(default_factory=dict)
    max_tries: int | None = None
    max_seconds: float | None = None
    search_mode: Literal["random", "exhaustive"] = "random"
    # Required when an exhaustive run's total exceeds the huge-space
    # threshold (10**12).  Ignored for random runs.
    confirm_huge: bool = False


class PrimitiveSearchRun(BaseModel):
    id: int
    family: str
    L: int
    k: int
    structural_params: dict[str, Any]
    fixed_params: dict[str, Any]
    max_tries: int | None
    max_seconds: float | None
    status: str
    tries_done: int
    found_count: int
    elapsed_seconds: float | None
    error_message: str | None
    created_at: str
    started_at: str | None
    finished_at: str | None
    # Exhaustive-mode bookkeeping; present (with defaults) for all
    # runs after the schema migration.
    search_mode: str = "random"
    enum_index:  int = 0
    enum_total:  str | None = None     # decimal; may exceed 2^63
    enum_axes:   list[dict] | None = None


class PrimitiveGenerator(BaseModel):
    id: int
    search_run_id: int | None
    family: str
    L: int
    k: int
    structural_params: dict[str, Any]
    search_params: dict[str, Any]
    all_params: dict[str, Any]
    found_at_try: int | None
    created_at: str


class ProgressUpdate(BaseModel):
    tries_done: int
    found_count: int
    current_info: dict[str, Any] | None
    message: str | None
    updated_at: str


# ── Tempering search ─────────────────────────────────────────────────────

class TemperingComponentSpec(BaseModel):
    generator_ids: list[int]
    tempering: list[dict[str, Any]] = Field(default_factory=list)
    shared_with_component: int | None = None


class TemperingSearchCreate(BaseModel):
    components: list[TemperingComponentSpec]
    test: dict[str, Any]            # {"type": "equidistribution", ...}
    Lmax: int
    nb_tries: int = 100
    optimizer: dict[str, Any] | None = None


class TemperingSearchRun(BaseModel):
    id: int
    test_type: str
    test_config: dict[str, Any]
    Lmax: int
    nb_tries: int
    optimizer_config: dict[str, Any] | None
    status: str
    combos_total: int | None
    combos_done: int
    best_se: int | None
    elapsed_seconds: float | None
    error_message: str | None
    created_at: str
    started_at: str | None
    finished_at: str | None


# ── Tested generator ────────────────────────────────────────────────────

class TestedGeneratorComponent(BaseModel):
    component_index: int
    generator_id: int | None
    family: str
    L: int
    k: int
    all_params: dict[str, Any]
    tempering_params: list[dict[str, Any]]


class TestResult(BaseModel):
    test_type: str
    test_config: dict[str, Any]
    se: int | None
    is_me: bool | None
    secf: int | None
    is_cf: bool | None
    score: float | None
    detail: dict[str, Any]
    elapsed_seconds: float | None
    created_at: str


class TestedGenerator(BaseModel):
    id: int
    search_run_id: int | None
    Lmax: int
    k_g: int
    J: int
    components: list[TestedGeneratorComponent]
    results: list[TestResult]
    created_at: str


# ── Listing/pagination wrapper ──────────────────────────────────────────

class Page(BaseModel):
    items: list
    total: int
    page: int
    per_page: int
