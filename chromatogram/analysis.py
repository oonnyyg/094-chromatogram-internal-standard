"""Batch analysis: input validation, OLS calibration, and quantification."""

import math

from .curve import validate_curve, window_area
from .errors import AnalysisError, InputError


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _parse_window(data, key):
    window = data.get(key)
    if (
        not isinstance(window, list)
        or len(window) != 2
        or not all(_is_number(v) for v in window)
    ):
        raise InputError(f"'{key}' must be [start, end] numbers")
    start, end = float(window[0]), float(window[1])
    if not (math.isfinite(start) and math.isfinite(end)):
        raise InputError(f"'{key}' bounds must be finite")
    if not end > start:
        raise InputError(f"'{key}' must have positive length")
    return start, end


def _parse_samples(data, key, require_concentration):
    samples = data.get(key)
    if not isinstance(samples, list):
        raise InputError(f"'{key}' must be a list")
    parsed = []
    for sample in samples:
        if not isinstance(sample, dict):
            raise InputError(f"'{key}': each sample must be an object")
        sid = sample.get("id")
        if not isinstance(sid, str) or not sid:
            raise InputError(f"'{key}': each sample needs a non-empty string 'id'")
        if "curve" not in sample:
            raise InputError(f"'{key}' sample '{sid}': missing 'curve'")
        curve = validate_curve(sample["curve"], f"'{key}' sample '{sid}'")
        concentration = None
        if require_concentration:
            conc = sample.get("concentration")
            if not _is_number(conc) or not math.isfinite(float(conc)) or float(conc) < 0.0:
                raise InputError(
                    f"'{key}' sample '{sid}': 'concentration' must be finite and non-negative"
                )
            concentration = float(conc)
        elif "concentration" in sample:
            raise InputError(
                f"'{key}' sample '{sid}': unknowns must not provide 'concentration'"
            )
        parsed.append({"id": sid, "curve": curve, "concentration": concentration})
    return parsed


def parse_input(data):
    """Validate the full input document and return a normalized structure."""
    if not isinstance(data, dict):
        raise InputError("top-level input must be a JSON object")
    analyte_window = _parse_window(data, "analyte_window")
    is_window = _parse_window(data, "internal_standard_window")
    if analyte_window[0] < is_window[1] and is_window[0] < analyte_window[1]:
        raise InputError("analyte and internal-standard windows must not overlap")
    standards = _parse_samples(data, "standards", require_concentration=True)
    unknowns = _parse_samples(data, "unknowns", require_concentration=False)
    ids = [s["id"] for s in standards + unknowns]
    if len(ids) != len(set(ids)):
        raise InputError("sample ids must be unique across standards and unknowns")
    distinct = {s["concentration"] for s in standards}
    if len(distinct) < 3:
        raise InputError("at least three standards with distinct concentrations are required")
    return {
        "analyte_window": analyte_window,
        "internal_standard_window": is_window,
        "standards": standards,
        "unknowns": unknowns,
    }


def _area_ratios(samples, analyte_window, is_window):
    for sample in samples:
        analyte = window_area(sample["curve"], *analyte_window)
        internal = window_area(sample["curve"], *is_window)
        if not internal > 0.0:
            raise AnalysisError(
                f"sample '{sample['id']}': internal-standard area is not positive; batch aborted"
            )
        sample["analyte_area"] = analyte
        sample["internal_standard_area"] = internal
        sample["area_ratio"] = analyte / internal


def _fit_ols(points):
    """Ordinary least squares y = slope * x + intercept for (x, y) points."""
    n = len(points)
    sx = sum(p[0] for p in points)
    sy = sum(p[1] for p in points)
    sxx = sum(p[0] * p[0] for p in points)
    sxy = sum(p[0] * p[1] for p in points)
    denom = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, intercept


def run_analysis(data):
    """Run the full batch; returns (calibration report, per-unknown results)."""
    parsed = parse_input(data)
    standards = parsed["standards"]
    unknowns = parsed["unknowns"]
    analyte_window = parsed["analyte_window"]
    is_window = parsed["internal_standard_window"]
    _area_ratios(standards, analyte_window, is_window)
    _area_ratios(unknowns, analyte_window, is_window)
    points = [(s["concentration"], s["area_ratio"]) for s in standards]
    slope, intercept = _fit_ols(points)
    if not slope > 0.0:
        raise AnalysisError("calibration slope is not positive; batch aborted")
    for s in standards:
        s["fitted_ratio"] = slope * s["concentration"] + intercept
        s["residual"] = s["area_ratio"] - s["fitted_ratio"]
    conc_min = min(p[0] for p in points)
    conc_max = max(p[0] for p in points)
    for u in unknowns:
        u["concentration"] = (u["area_ratio"] - intercept) / slope
        u["extrapolated"] = bool(u["concentration"] < conc_min or u["concentration"] > conc_max)
    calibration = {
        "slope": slope,
        "intercept": intercept,
        "concentration_range": [conc_min, conc_max],
        "standards": [
            {
                "id": s["id"],
                "concentration": s["concentration"],
                "analyte_area": s["analyte_area"],
                "internal_standard_area": s["internal_standard_area"],
                "area_ratio": s["area_ratio"],
                "fitted_ratio": s["fitted_ratio"],
                "residual": s["residual"],
            }
            for s in standards
        ],
    }
    results = [
        {
            "id": u["id"],
            "analyte_area": u["analyte_area"],
            "internal_standard_area": u["internal_standard_area"],
            "area_ratio": u["area_ratio"],
            "concentration": u["concentration"],
            "extrapolated": u["extrapolated"],
        }
        for u in unknowns
    ]
    return calibration, results
