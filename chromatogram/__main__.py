"""Command-line interface: python -m chromatogram input.json --out OUTPUT_DIR"""

import argparse
import csv
import json
import sys
from pathlib import Path

from .analysis import run_analysis
from .errors import AnalysisError, InputError

CSV_FIELDS = [
    "id",
    "analyte_area",
    "internal_standard_area",
    "area_ratio",
    "concentration",
    "extrapolated",
]


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="chromatogram",
        description="Internal-standard quantification for chromatography batches.",
    )
    parser.add_argument("input", help="path to the input JSON file")
    parser.add_argument("--out", required=True, help="output directory for reports")
    args = parser.parse_args(argv)

    try:
        with open(args.input, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"error: cannot read input JSON: {exc}", file=sys.stderr)
        return 1

    try:
        calibration, results = run_analysis(data)
    except (InputError, AnalysisError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "results.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    json_path = out_dir / "calibration.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(calibration, fh, indent=2, ensure_ascii=False)
        fh.write("\n")

    print(f"wrote {csv_path} and {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
