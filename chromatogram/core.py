"""Window clipping, baseline-corrected integration and OLS calibration."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .inputspec import BatchInput, Sample, Standard, Window


class QuantitationError(ValueError):
    """Raised when the batch cannot be quantitated (e.g. non-positive IS area)."""


@dataclass(frozen=True)
class IntegrationResult:
    area: float


@dataclass(frozen=True)
class SampleResult:
    id: str
    analyte_area: float
    internal_standard_area: float
    ratio: float


@dataclass(frozen=True)
class Calibration:
    slope: float
    intercept: float
    residuals: dict[str, float]


def _value_at(times: tuple[float, ...], values: tuple[float, ...], t: float) -> float:
    if t <= times[0]:
        return values[0]
    if t >= times[-1]:
        return values[-1]
    hi = 1
    while times[hi] < t:
        hi += 1
    lo = hi - 1
    t0, t1 = times[lo], times[hi]
    if t == t0:
        return values[lo]
    if t == t1:
        return values[hi]
    frac = (t - t0) / (t1 - t0)
    return values[lo] + frac * (values[hi] - values[lo])


def _clip_to_window(times: tuple[float, ...], values: tuple[float, ...], window: Window
                    ) -> tuple[list[float], list[float]]:
    """Return nodes strictly inside the window plus interpolated boundaries."""
    pts_t: list[float] = [window.start]
    pts_y: list[float] = [_value_at(times, values, window.start)]
    for t, y in zip(times, values):
        if window.start < t < window.end:
            pts_t.append(t)
            pts_y.append(y)
    pts_t.append(window.end)
    pts_y.append(_value_at(times, values, window.end))
    return pts_t, pts_y


def integrate_window(times: tuple[float, ...], values: tuple[float, ...], window: Window) -> float:
    """Integrate the positive part of the baseline-corrected signal.

    Baseline is the straight chord joining the signals at the window edges.
    Segments crossing zero are split at the crossing so a piecewise-negative
    region never contributes; this is not the same as clamping node values.
    """
    pts_t, pts_y = _clip_to_window(times, values, window)
    t0, t1 = pts_t[0], pts_t[-1]
    y0, y1 = pts_y[0], pts_y[-1]
    span = t1 - t0
    corrected = [y - (y0 + (y1 - y0) * (t - t0) / span) for t, y in zip(pts_t, pts_y)]

    total_parts: list[float] = []
    for i in range(len(pts_t) - 1):
        a, b = corrected[i], corrected[i + 1]
        dt = pts_t[i + 1] - pts_t[i]
        if a >= 0.0 and b >= 0.0:
            total_parts.append(0.5 * (a + b) * dt)
        elif a <= 0.0 and b <= 0.0:
            continue
        else:
            # Strict sign change: zero at fraction -a/(b-a) along the segment.
            frac = -a / (b - a)
            if a > 0.0:
                total_parts.append(0.5 * a * (frac * dt))
            else:
                total_parts.append(0.5 * b * ((1.0 - frac) * dt))
    return math.fsum(total_parts)


def _sample_result(sample: Sample, batch: BatchInput) -> SampleResult:
    analyte_area = integrate_window(sample.times, sample.values, batch.analyte_window)
    is_area = integrate_window(sample.times, sample.values, batch.is_window)
    if not is_area > 0.0:
        raise QuantitationError(
            f"sample {sample.id!r}: internal-standard area must be positive, got {is_area}"
        )
    return SampleResult(sample.id, analyte_area, is_area, analyte_area / is_area)


def calibrate(standards: tuple[Standard, ...], ratios: dict[str, float]) -> Calibration:
    """Ordinary least squares of ratio (y) on concentration (x), with intercept."""
    xs = [s.concentration for s in standards]
    ys = [ratios[s.id] for s in standards]
    n = len(xs)
    mx = math.fsum(xs) / n
    my = math.fsum(ys) / n
    sxx = math.fsum((x - mx) ** 2 for x in xs)
    sxy = math.fsum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    if not slope > 0.0:
        raise QuantitationError(f"calibration slope must be positive, got {slope}")
    residuals = {s.id: ratios[s.id] - (slope * s.concentration + intercept) for s in standards}
    return Calibration(slope, intercept, residuals)


@dataclass(frozen=True)
class Quantitation:
    standards: tuple[SampleResult, ...]
    unknowns: tuple[SampleResult, ...]
    calibration: Calibration
    concentrations: dict[str, float]
    extrapolated: dict[str, bool]
    concentration_range: tuple[float, float]


def process_batch(batch: BatchInput) -> Quantitation:
    """Integrate every sample, fit standards and back-calculate unknowns."""
    standard_results = tuple(_sample_result(s, batch) for s in batch.standards)
    ratios = {r.id: r.ratio for r in standard_results}
    calibration = calibrate(batch.standards, ratios)

    unknown_results = tuple(_sample_result(s, batch) for s in batch.unknowns)
    cmin = min(s.concentration for s in batch.standards)
    cmax = max(s.concentration for s in batch.standards)
    concentrations: dict[str, float] = {}
    extrapolated: dict[str, bool] = {}
    for result in unknown_results:
        concentration = (result.ratio - calibration.intercept) / calibration.slope
        concentrations[result.id] = concentration
        # Equality with a calibration boundary is not extrapolation.
        extrapolated[result.id] = concentration < cmin or concentration > cmax

    return Quantitation(
        standard_results,
        unknown_results,
        calibration,
        concentrations,
        extrapolated,
        (cmin, cmax),
    )
