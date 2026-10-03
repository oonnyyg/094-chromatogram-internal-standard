"""Input loading and validation for chromatogram batches."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


class InputError(ValueError):
    """Raised when the input JSON is missing fields or violates constraints."""


@dataclass(frozen=True)
class Window:
    start: float
    end: float


@dataclass(frozen=True)
class Sample:
    id: str
    times: tuple[float, ...]
    values: tuple[float, ...]


@dataclass(frozen=True)
class Standard(Sample):
    concentration: float = 0.0


@dataclass(frozen=True)
class BatchInput:
    analyte_window: Window
    is_window: Window
    standards: tuple[Standard, ...]
    unknowns: tuple[Sample, ...]


def _is_real(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _finite(value: Any) -> float:
    if not _is_real(value) or not math.isfinite(float(value)):
        raise InputError(f"expected a finite number, got {value!r}")
    return float(value)


def _window(raw: Any, name: str) -> Window:
    if not isinstance(raw, dict):
        raise InputError(f"window {name!r} must be an object")
    try:
        start = _finite(raw["start"])
        end = _finite(raw["end"])
    except KeyError as exc:
        raise InputError(f"window {name!r} missing key {exc.args[0]!r}") from None
    if not end > start:
        raise InputError(f"window {name!r} must have positive length")
    return Window(start, end)


def _curve(raw: Any, where: str, coverage: tuple[float, float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if not isinstance(raw, list) or len(raw) < 2:
        raise InputError(f"{where}: curve needs at least 2 sampling points")
    times: list[float] = []
    values: list[float] = []
    prev: float | None = None
    for index, point in enumerate(raw):
        if not isinstance(point, (list, tuple)) or len(point) != 2:
            raise InputError(f"{where}: point {index} must be a [time, signal] pair")
        t = _finite(point[0])
        y = _finite(point[1])
        if prev is not None and not t > prev:
            raise InputError(f"{where}: times must be strictly increasing")
        prev = t
        times.append(t)
        values.append(y)
    lo, hi = coverage
    if times[0] > lo or times[-1] < hi:
        raise InputError(f"{where}: curve must fully cover both time windows")
    return tuple(times), tuple(values)


def _sample(raw: Any, index: int, coverage: tuple[float, float], *, allow_concentration: bool) -> Sample:
    where = f"sample {index}"
    if not isinstance(raw, dict):
        raise InputError(f"{where}: must be an object")
    sid = raw.get("id")
    if not isinstance(sid, str) or not sid:
        raise InputError(f"{where}: id must be a non-empty string")
    where = f"sample {sid!r}"
    if "curve" not in raw:
        raise InputError(f"{where}: missing curve")
    times, values = _curve(raw["curve"], where, coverage)
    has_conc = "concentration" in raw and raw["concentration"] is not None
    if allow_concentration:
        if not has_conc:
            raise InputError(f"{where}: standard requires a concentration")
        concentration = _finite(raw["concentration"])
        if concentration < 0.0:
            raise InputError(f"{where}: concentration must be non-negative")
        return Standard(sid, times, values, concentration)
    if has_conc:
        raise InputError(f"{where}: unknown samples must not provide a concentration")
    return Sample(sid, times, values)


def load_input(path: str | Path) -> BatchInput:
    """Read and validate the batch JSON file."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
    except FileNotFoundError as exc:
        raise InputError(f"input file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise InputError(f"invalid JSON: {exc}") from exc

    if not isinstance(raw, dict):
        raise InputError("top-level JSON must be an object")
    windows = raw.get("windows")
    if not isinstance(windows, dict):
        raise InputError("missing 'windows' object")
    analyte = _window(windows.get("analyte"), "analyte")
    internal = _window(windows.get("internal_standard"), "internal_standard")
    # Interiors must be disjoint; touching at an endpoint is allowed.
    if not (analyte.end <= internal.start or internal.end <= analyte.start):
        raise InputError("analyte and internal-standard windows overlap")

    coverage = (min(analyte.start, internal.start), max(analyte.end, internal.end))

    raw_standards = raw.get("standards")
    if not isinstance(raw_standards, list):
        raise InputError("'standards' must be a list")
    standards = tuple(
        _sample(item, i, coverage, allow_concentration=True)
        for i, item in enumerate(raw_standards)
    )
    if len({s.concentration for s in standards}) < 3:
        raise InputError("at least three distinct standard concentrations are required")

    raw_unknowns = raw.get("unknowns", [])
    if not isinstance(raw_unknowns, list):
        raise InputError("'unknowns' must be a list")
    unknowns = tuple(
        _sample(item, i, coverage, allow_concentration=False)
        for i, item in enumerate(raw_unknowns)
    )

    ids = [s.id for s in standards] + [s.id for s in unknowns]
    if len(set(ids)) != len(ids):
        raise InputError("sample ids must be unique")

    return BatchInput(analyte, internal, standards, unknowns)
