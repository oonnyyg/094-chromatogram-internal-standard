import unittest

from chromatogram.analysis import run_analysis
from chromatogram.errors import AnalysisError, InputError


def make_curve(analyte_height, is_height):
    """Triangle peaks: analyte in [0,1] peaking at 0.5, IS in [2,3] at 2.5."""
    return [
        [-0.5, 0.0],
        [0.0, 0.0],
        [0.5, analyte_height],
        [1.0, 0.0],
        [1.5, 0.0],
        [2.0, 0.0],
        [2.5, is_height],
        [3.0, 0.0],
        [3.5, 0.0],
    ]


def make_sample(sid, ratio, concentration=None, is_height=2.0):
    sample = {"id": sid, "curve": make_curve(ratio * is_height, is_height)}
    if concentration is not None:
        sample["concentration"] = concentration
    return sample


SLOPE = 2.0
INTERCEPT = 0.5


def make_input(standard_concs=(0.0, 1.0, 2.0, 3.0), unknown_ratios=(3.5,)):
    standards = [
        make_sample(f"std{i}", SLOPE * c + INTERCEPT, concentration=c)
        for i, c in enumerate(standard_concs)
    ]
    unknowns = [make_sample(f"unk{i}", r) for i, r in enumerate(unknown_ratios)]
    return {
        "analyte_window": [0.0, 1.0],
        "internal_standard_window": [2.0, 3.0],
        "standards": standards,
        "unknowns": unknowns,
    }


class RunAnalysisTests(unittest.TestCase):
    def test_recovers_calibration_and_concentration(self):
        calibration, results = run_analysis(make_input())
        self.assertAlmostEqual(calibration["slope"], SLOPE, places=12)
        self.assertAlmostEqual(calibration["intercept"], INTERCEPT, places=12)
        self.assertEqual(len(calibration["standards"]), 4)
        for entry in calibration["standards"]:
            self.assertAlmostEqual(entry["residual"], 0.0, places=12)
        self.assertEqual(len(results), 1)
        self.assertAlmostEqual(results[0]["concentration"], 1.5, places=12)
        self.assertFalse(results[0]["extrapolated"])
        self.assertAlmostEqual(results[0]["analyte_area"], 0.5 * 1.0 * 7.0, places=12)
        self.assertAlmostEqual(results[0]["internal_standard_area"], 0.5 * 1.0 * 2.0, places=12)
        self.assertAlmostEqual(results[0]["area_ratio"], 3.5, places=12)

    def test_negative_concentration_kept(self):
        _, results = run_analysis(make_input(unknown_ratios=(0.1,)))
        self.assertLess(results[0]["concentration"], 0.0)
        self.assertTrue(results[0]["extrapolated"])

    def test_boundary_not_extrapolated(self):
        ratios = (INTERCEPT, SLOPE * 3.0 + INTERCEPT)
        _, results = run_analysis(make_input(unknown_ratios=ratios))
        self.assertAlmostEqual(results[0]["concentration"], 0.0, places=12)
        self.assertAlmostEqual(results[1]["concentration"], 3.0, places=12)
        self.assertFalse(results[0]["extrapolated"])
        self.assertFalse(results[1]["extrapolated"])

    def test_residuals_measured_minus_fitted(self):
        data = make_input()
        data["standards"][1]["curve"] = make_curve((SLOPE * 1.0 + INTERCEPT + 0.4) * 2.0, 2.0)
        calibration, _ = run_analysis(data)
        for entry in calibration["standards"]:
            expected = entry["area_ratio"] - (
                calibration["slope"] * entry["concentration"] + calibration["intercept"]
            )
            self.assertAlmostEqual(entry["residual"], expected, places=12)

    def test_non_positive_is_area_aborts(self):
        data = make_input()
        data["unknowns"][0]["curve"] = make_curve(7.0, 0.0)
        with self.assertRaises(AnalysisError):
            run_analysis(data)

    def test_non_positive_slope_aborts(self):
        data = make_input(standard_concs=(0.0, 1.0, 2.0))
        for s in data["standards"]:
            s["curve"] = make_curve(1.0, 2.0)
        with self.assertRaises(AnalysisError):
            run_analysis(data)


class ValidationTests(unittest.TestCase):
    def test_overlapping_windows_rejected(self):
        data = make_input()
        data["internal_standard_window"] = [0.5, 2.5]
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_zero_length_window_rejected(self):
        data = make_input()
        data["analyte_window"] = [0.5, 0.5]
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_duplicate_ids_rejected(self):
        data = make_input()
        data["unknowns"][0]["id"] = data["standards"][0]["id"]
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_too_few_distinct_standard_concentrations(self):
        data = make_input(standard_concs=(1.0, 1.0, 2.0))
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_negative_standard_concentration_rejected(self):
        data = make_input(standard_concs=(-1.0, 1.0, 2.0))
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_unknown_with_concentration_rejected(self):
        data = make_input()
        data["unknowns"][0]["concentration"] = 1.0
        with self.assertRaises(InputError):
            run_analysis(data)

    def test_window_not_covered_rejected(self):
        data = make_input()
        data["analyte_window"] = [-1.0, 0.5]
        with self.assertRaises(InputError):
            run_analysis(data)


if __name__ == "__main__":
    unittest.main()
