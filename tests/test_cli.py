from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _curve(analyte_height: float, is_height: float = 8.0) -> list[list[float]]:
        # Flat baseline, flat-topped analyte peak [12,18], IS peak [32,38].
        # Windowed trapezoidal area is 8 * height (windows span [10,20]/[30,40]).
    return [
        [0.0, 0.0],
        [10.0, 0.0],
        [12.0, analyte_height],
        [18.0, analyte_height],
        [20.0, 0.0],
        [30.0, 0.0],
        [32.0, is_height],
        [38.0, is_height],
        [40.0, 0.0],
    ]


def _payload() -> dict:
    # ratio = 3*h_a/24 = h_a/8 ; standards satisfy h_a/8 = 0.2 c + 0.1
    heights = {0.0: 0.8, 1.0: 2.4, 2.0: 4.0}
    standards = [
        {"id": f"std{c}", "concentration": c, "curve": _curve(h)}
        for c, h in heights.items()
    ]
    return {
        "windows": {"analyte": {"start": 10.0, "end": 20.0},
                    "internal_standard": {"start": 30.0, "end": 40.0}},
        "standards": standards,
        "unknowns": [
            {"id": "unk-in", "curve": _curve(2.4)},       # ratio 0.3 -> c = 1.0
            {"id": "unk-high", "curve": _curve(6.0)},     # ratio 0.75 -> c = 3.25
            {"id": "unk-neg", "curve": _curve(0.1)},      # ratio 0.0125 -> c = -0.4375
        ],
    }


def _run(input_path: Path, out_dir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "chromatogram", str(input_path), "--out", str(out_dir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


class CliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(self.id().replace(".", "_"))

    def tearDown(self) -> None:
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_end_to_end_reports(self) -> None:
        work = self.tmp / "ok"
        work.mkdir(parents=True)
        payload = _payload()
        input_path = work / "input.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        out_dir = work / "reports"

        proc = _run(input_path, out_dir)
        self.assertEqual(proc.returncode, 0, proc.stderr)

        calib = json.loads((out_dir / "calibration.json").read_text(encoding="utf-8"))
        self.assertAlmostEqual(calib["slope"], 0.2, places=12)
        self.assertAlmostEqual(calib["intercept"], 0.1, places=12)
        for entry, c in zip(calib["standards"], (0.0, 1.0, 2.0)):
            self.assertAlmostEqual(entry["area_ratio"], 0.2 * c + 0.1, places=12)
            self.assertAlmostEqual(entry["residual"], 0.0, places=12)

        with open(out_dir / "standards.csv", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(float(rows[1]["analyte_area"]), 19.2, places=10)
        self.assertAlmostEqual(float(rows[0]["internal_standard_area"]), 64.0, places=10)

        with open(out_dir / "unknowns.csv", encoding="utf-8", newline="") as handle:
            unknowns = {row["sample_id"]: row for row in csv.DictReader(handle)}
        self.assertAlmostEqual(float(unknowns["unk-in"]["concentration"]), 1.0, places=12)
        self.assertEqual(unknowns["unk-in"]["outside_calibration_range"], "false")
        self.assertAlmostEqual(float(unknowns["unk-high"]["concentration"]), 3.25, places=12)
        self.assertEqual(unknowns["unk-high"]["outside_calibration_range"], "true")
        neg_conc = float(unknowns["unk-neg"]["concentration"])
        self.assertAlmostEqual(neg_conc, -0.4375, places=12)
        self.assertEqual(unknowns["unk-neg"]["outside_calibration_range"], "true")

    def _assert_fails_without_reports(self, mutate) -> None:
        work = self.tmp / "bad"
        work.mkdir(parents=True)
        payload = _payload()
        mutate(payload)
        input_path = work / "input.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        out_dir = work / "reports"

        proc = _run(input_path, out_dir)
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse(out_dir.exists(), "no reports may be written for invalid input")

    def test_overlapping_windows_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            payload["windows"]["internal_standard"]["start"] = 15.0
        self._assert_fails_without_reports(mutate)

    def test_non_positive_is_area_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            for std in payload["standards"]:
                for point in std["curve"]:
                    if point[0] >= 30.0:
                        point[1] = 0.0
        self._assert_fails_without_reports(mutate)

    def test_too_few_concentrations_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            payload["standards"][2]["concentration"] = 1.0
        self._assert_fails_without_reports(mutate)

    def test_unknown_with_concentration_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            payload["unknowns"][0]["concentration"] = 1.0
        self._assert_fails_without_reports(mutate)

    def test_non_increasing_times_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            payload["standards"][0]["curve"][3][0] = 11.0
        self._assert_fails_without_reports(mutate)

    def test_duplicate_ids_rejected(self) -> None:
        def mutate(payload: dict) -> None:
            payload["unknowns"][0]["id"] = payload["standards"][0]["id"]
        self._assert_fails_without_reports(mutate)

    def test_numeric_tolerance_on_synthetic_linear_data(self) -> None:
        work = self.tmp / "tol"
        work.mkdir(parents=True)
        input_path = work / "input.json"
        input_path.write_text(json.dumps(_payload()), encoding="utf-8")
        out_dir = work / "reports"
        proc = _run(input_path, out_dir)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        calib = json.loads((out_dir / "calibration.json").read_text(encoding="utf-8"))
        tolerance = 1e-8 * max(1.0, abs(0.2))
        self.assertTrue(math.isfinite(calib["slope"]))
        self.assertLessEqual(abs(calib["slope"] - 0.2), tolerance)


if __name__ == "__main__":
    unittest.main()
