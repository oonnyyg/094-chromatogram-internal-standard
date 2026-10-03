"""Curve handling: validation, interpolation, and baseline-corrected integration."""

import math

from .errors import InputError


def validate_curve(points, context):
    """Validate a curve given as a list of [time, signal] pairs.

    Returns a list of (time, signal) tuples with strictly increasing times
    and at least two finite samples.
    """
    if not isinstance(points, list) or len(points) < 2:
        raise InputError(f"{context}: curve must contain at least two sample points")
    curve = []
    for point in points:
        if (
            not isinstance(point, list)
            or len(point) != 2
            or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in point)
        ):
            raise InputError(f"{context}: each sample point must be [time, signal] numbers")
        t, s = float(point[0]), float(point[1])
        if not (math.isfinite(t) and math.isfinite(s)):
            raise InputError(f"{context}: sample points must be finite")
        curve.append((t, s))
    for i in range(1, len(curve)):
        if not curve[i][0] > curve[i - 1][0]:
            raise InputError(f"{context}: sample times must be strictly increasing")
    return curve


def interpolate(curve, t):
    """Linear interpolation of the curve at time t (t must be within range)."""
    if t < curve[0][0] or t > curve[-1][0]:
        raise InputError("interpolation time outside curve range")
    for i in range(1, len(curve)):
        t0, s0 = curve[i - 1]
        t1, s1 = curve[i]
        if t <= t1:
            if t1 == t0:
                return s1
            frac = (t - t0) / (t1 - t0)
            return s0 + frac * (s1 - s0)
    return curve[-1][1]


def _positive_segment_area(dt, y0, y1):
    """Area of the positive part of a linear segment from y0 to y1 over dt.

    Zero crossings are handled analytically by splitting the segment at the
    crossing point, so negative values are never simply clamped to zero.
    """
    if y0 <= 0.0 and y1 <= 0.0:
        return 0.0
    if y0 >= 0.0 and y1 >= 0.0:
        return 0.5 * (y0 + y1) * dt
    # One endpoint positive, one negative: split at the zero crossing.
    frac = y0 / (y0 - y1)  # fraction of dt where the signal crosses zero
    if y0 > 0.0:
        return 0.5 * y0 * (dt * frac)
    return 0.5 * y1 * (dt * (1.0 - frac))


def window_area(curve, start, end):
    """Baseline-corrected positive area of the curve over [start, end].

    The baseline is the straight line joining the (interpolated) signal at
    the window boundaries. After subtracting it, only the positive part of
    the residual is integrated, splitting segments at zero crossings.
    """
    if start < curve[0][0] or end > curve[-1][0]:
        raise InputError(
            f"window [{start}, {end}] is not fully covered by curve "
            f"range [{curve[0][0]}, {curve[-1][0]}]"
        )
    b0 = interpolate(curve, start)
    b1 = interpolate(curve, end)
    span = end - start

    # Collect sub-segment boundaries: window ends plus interior nodes.
    times = [start]
    times.extend(t for t, _ in curve if start < t < end)
    times.append(end)

    area = 0.0
    for i in range(1, len(times)):
        ta, tb = times[i - 1], times[i]
        ya = interpolate(curve, ta) - (b0 + (b1 - b0) * (ta - start) / span)
        yb = interpolate(curve, tb) - (b0 + (b1 - b0) * (tb - start) / span)
        area += _positive_segment_area(tb - ta, ya, yb)
    return area
