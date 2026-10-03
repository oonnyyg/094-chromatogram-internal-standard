import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.test_analysis import make_input


class CliTests(unittest.TestCase):
    def run_cli(self, data, out_dir):
        with tempfile.NamedTemporaryFile(
            "w", suffix=".json", delete=False, encoding="utf-8"
        ) as fh:
            json.dump(data, fh)
            input_path = fh.name
        return subprocess.run(
            [sys.executable, "-m", "chromatogram", input_path, "--out", str(out_dir)],
            capture_output=True,
            text=True,
        )

    def test_success_writes_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            proc = self.run_cli(make_input(), out_dir)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            csv_path = out_dir / "results.csv"
            json_path = out_dir / "calibration.json"
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            with open(csv_path, newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["id"], "unk0")
            self.assertAlmostEqual(float(rows[0]["concentration"]), 1.5, places=12)
            self.assertEqual(rows[0]["extrapolated"], "False")
            calibration = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertAlmostEqual(calibration["slope"], 2.0, places=12)
            self.assertAlmostEqual(calibration["intercept"], 0.5, places=12)
            self.assertEqual(len(calibration["standards"]), 4)

    def test_invalid_input_fails_without_reports(self):
        data = make_input()
        data["standards"] = data["standards"][:2]
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            proc = self.run_cli(data, out_dir)
            self.assertNotEqual(proc.returncode, 0)
            self.assertFalse((out_dir / "results.csv").exists())
            self.assertFalse((out_dir / "calibration.json").exists())

    def test_missing_input_file_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            proc = subprocess.run(
                [sys.executable, "-m", "chromatogram", "no-such.json", "--out", str(out_dir)],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertFalse((out_dir / "results.csv").exists())


if __name__ == "__main__":
    unittest.main()
