from __future__ import annotations

import math
import unittest

from chromatogram.core import Window, integrate_window


class IntegrateWindowTests(unittest.TestCase):
    def test_triangle_with_interpolated_boundaries(self) -> None:
        # Window edges fall between nodes; baseline chord is flat at zero.
        times = (0.0, 4.0, 6.0, 10.0)
        values = (0.0, 4.0, 4.0, 0.0)
        area = integrate_window(times, values, Window(2.0, 8.0))
        # Window-edge chord is flat at 2; corrected heights 0, 2, 2, 0.
        self.assertAlmostEqual(area, 8.0, places=12)

    def test_slanted_baseline_is_subtracted(self) -> None:
        # Constant 1.0 signal over a window whose chord is also 1.0 -> area 0.
        times = (0.0, 5.0, 10.0)
        values = (1.0, 1.0, 1.0)
        self.assertEqual(integrate_window(times, values, Window(0.0, 10.0)), 0.0)

        # Linear ramp equals its own chord -> area 0 even with slope.
        values2 = (0.0, 5.0, 10.0)
        self.assertEqual(integrate_window(times, values2, Window(0.0, 10.0)), 0.0)

    def test_negative_parts_excluded_and_crossing_split(self) -> None:
        # Endpoints sit on the baseline chord (zero here); internal values
        # +4 then -4 make the signal cross zero inside a segment.
        times = (0.0, 5.0, 10.0, 15.0)
        values = (0.0, 4.0, -4.0, 0.0)
        area = integrate_window(times, values, Window(0.0, 15.0))
        # Triangle [0,5]: 10; positive half of [5,10]: 4*2.5/2 = 5.
        self.assertAlmostEqual(area, 15.0, places=12)

    def test_diagonal_crossing_between_internal_nodes(self) -> None:
        # Corrected: 0, +1, -1, 0 across 2-wide segments.
        times = (0.0, 2.0, 4.0, 6.0)
        values = (0.0, 1.0, -1.0, 0.0)
        area = integrate_window(times, values, Window(0.0, 6.0))
        # [0,2] trapezoid 1; positive half of [2,4]: 0.5 -> 1.5.
        # Clamping node values first would instead yield 2.0.
        self.assertAlmostEqual(area, 1.5, places=12)

    def test_all_negative_gives_zero(self) -> None:
        times = (0.0, 1.0, 2.0)
        values = (-3.0, -5.0, -4.0)
        area = integrate_window(times, values, Window(0.0, 2.0))
        self.assertEqual(area, 0.0)
        self.assertTrue(math.isfinite(area))


if __name__ == "__main__":
    unittest.main()
