"""CSV / JSON report rendering."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .core import Quantitation, SampleResult
from .inputspec import BatchInput


def _sample_row(result: SampleResult) -> dict[str, object]:
    return {
        "sample_id": result.id,
        "analyte_area": repr(result.analyte_area),
        "internal_standard_area": repr(result.internal_standard_area),
        "area_ratio": repr(result.ratio),
    }


def write_reports(batch: BatchInput, quant: Quantitation, out_dir: str | Path) -> dict[str, Path]:
    """Write per-role CSV files and the calibration JSON. Returns their paths."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    standards_path = out / "standards.csv"
    with open(standards_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "concentration",
                "analyte_area",
                "internal_standard_area",
                "area_ratio",
                "residual",
            ],
        )
        writer.writeheader()
        for sample, result in zip(batch.standards, quant.standards):
            row = _sample_row(result)
            row["concentration"] = repr(sample.concentration)
            row["residual"] = repr(quant.calibration.residuals[sample.id])
            writer.writerow(row)

    unknowns_path = out / "unknowns.csv"
    with open(unknowns_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "analyte_area",
                "internal_standard_area",
                "area_ratio",
                "concentration",
                "outside_calibration_range",
            ],
        )
        writer.writeheader()
        for result in quant.unknowns:
            row = _sample_row(result)
            row["concentration"] = repr(quant.concentrations[result.id])
            row["outside_calibration_range"] = "true" if quant.extrapolated[result.id] else "false"
            writer.writerow(row)

    calibration_path = out / "calibration.json"
    payload = {
        "model": "ratio = slope * concentration + intercept",
        "slope": quant.calibration.slope,
        "intercept": quant.calibration.intercept,
        "concentration_range": {
            "min": quant.concentration_range[0],
            "max": quant.concentration_range[1],
        },
        "standards": [
            {
                "sample_id": sample.id,
                "concentration": sample.concentration,
                "area_ratio": result.ratio,
                "fitted_ratio": quant.calibration.slope * sample.concentration + quant.calibration.intercept,
                "residual": quant.calibration.residuals[sample.id],
            }
            for sample, result in zip(batch.standards, quant.standards)
        ],
    }
    with open(calibration_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")

    return {
        "standards": standards_path,
        "unknowns": unknowns_path,
        "calibration": calibration_path,
    }
