import unittest

from chromatogram.curve import interpolate, validate_curve, window_area
from chromatogram.errors import InputError


class ValidateCurveTests(unittest.TestCase):
    def test_rejects_single_point(self):
        with self.assertRaises(InputError):
            validate_curve([[0.0, 1.0]], "ctx")

    def test_rejects_non_increasing_time(self):
        with self.assertRaises(InputError):
            validate_curve([[0.0, 1.0], [0.0, 2.0]], "ctx")
        with self.assertRaises(InputError):
            validate_curve([[1.0, 1.0], [0.5, 2.0]], "ctx")

    def test_rejects_non_finite(self):
        with self.assertRaises(InputError):
            validate_curve([[0.0, 1.0], [1.0, float("nan")]], "ctx")


class InterpolateTests(unittest.TestCase):
    def test_linear_between_nodes(self):
        curve = validate_curve([[0.0, 0.0], [2.0, 4.0]], "ctx")
        self.assertAlmostEqual(interpolate(curve, 0.5), 1.0)
        self.assertAlmostEqual(interpolate(curve, 2.0), 4.0)


class WindowAreaTests(unittest.TestCase):
    def test_flat_curve_zero_area(self):
        curve = validate_curve([[0.0, 5.0], [10.0, 5.0]], "ctx")
        self.assertAlmostEqual(window_area(curve, 2.0, 8.0), 0.0)

    def test_triangle_above_linear_baseline(self):
        # Peak of height 4 at t=5 on top of a constant offset of 1.
        curve = validate_curve(
            [[0.0, 1.0], [5.0, 5.0], [10.0, 1.0]], "ctx"
        )
        # Baseline is flat at 1; residual is a triangle base 10 height 4.
        self.assertAlmostEqual(window_area(curve, 0.0, 10.0), 0.5 * 10.0 * 4.0)

    def test_interpolated_window_boundaries(self):
        # Signal 2t; window [1, 3] -> baseline equals signal, area 0.
        curve = validate_curve([[0.0, 0.0], [4.0, 8.0]], "ctx")
        self.assertAlmostEqual(window_area(curve, 1.0, 3.0), 0.0)

    def test_zero_crossing_split_not_clamped(self):
        # Residual crosses zero inside the window: signal 1 at t=0, -1 at t=1,
        # 1 at t=2 with a flat-zero baseline curve segment around it.
        curve = validate_curve(
            [[-1.0, 0.0], [0.0, 1.0], [1.0, -1.0], [2.0, 1.0], [3.0, 0.0]], "ctx"
        )
        # Window [0, 2]: baseline joins (0,1) and (2,1) -> flat at 1.
        # Residual: 0, -2, 0 -> positive part is zero.
        self.assertAlmostEqual(window_area(curve, 0.0, 2.0), 0.0)
        # Window [-1, 0]: baseline joins (-1,0) and (0,1); residual is the
        # signal minus baseline; signal 0->1 linearly equals baseline, so 0.
        self.assertAlmostEqual(window_area(curve, -1.0, 0.0), 0.0)

    def test_positive_part_with_crossing(self):
        # Baseline flat at 0 via window endpoints; residual +2 then -2.
        curve = validate_curve(
            [[0.0, 0.0], [1.0, 2.0], [2.0, -2.0], [3.0, 0.0]], "ctx"
        )
        # Window [0, 3]: baseline joins (0,0) and (3,0) -> 0.
        # Positive part: triangle [0,1.5] height 2 -> 1.5.
        self.assertAlmostEqual(window_area(curve, 0.0, 3.0), 1.5)

    def test_window_outside_curve_rejected(self):
        curve = validate_curve([[0.0, 0.0], [1.0, 1.0]], "ctx")
        with self.assertRaises(InputError):
            window_area(curve, 0.5, 2.0)


if __name__ == "__main__":
    unittest.main()
