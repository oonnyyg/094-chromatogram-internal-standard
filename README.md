# chromatogram-internal-standard

Internal-standard quantitation for chromatographic experiments, using only the
Python standard library (Python 3.12+).

## Usage

```
python -m chromatogram INPUT.json --out OUTPUT_DIR
```

The tool reads one JSON file describing shared time windows, calibration
standards and unknown samples, and writes three reports:

- `standards.csv` — per-standard analyte/IS areas, area ratio, concentration
  and calibration residual (measured ratio minus fitted ratio).
- `unknowns.csv` — per-unknown areas, ratio, back-calculated concentration
  (negative values are kept, never clamped to zero) and an
  `outside_calibration_range` flag.
- `calibration.json` — the OLS model `ratio = slope * concentration +
  intercept`, the calibration concentration range, and each standard's
  measured ratio, fitted ratio and residual.

Run the example with `python -m chromatogram examples/batch.json --out out`.

## Input format

```json
{
  "windows": {
    "analyte": {"start": 10.0, "end": 20.0},
    "internal_standard": {"start": 30.0, "end": 40.0}
  },
  "standards": [
    {"id": "std-0", "concentration": 0.0, "curve": [[t, signal], ...]}
  ],
  "unknowns": [
    {"id": "unk-a", "curve": [[t, signal], ...]}
  ]
}
```

Rules enforced on input (any violation aborts the whole batch with a non-zero
exit code and writes no reports):

- Sample ids are unique and non-empty.
- Every curve has at least two finite points with strictly increasing times and
  fully covers both windows.
- Both windows have positive length and do not overlap (touching endpoints is
  allowed).
- Standard concentrations are finite and non-negative; at least three distinct
  concentrations are required.
- Unknown samples must not carry a `concentration`.

## Integration and baseline convention

Curves vary linearly between nodes. Window boundaries that do not land on a
node are evaluated by linear interpolation.

For each window the baseline is the straight chord connecting the (possibly
interpolated) signals at the two window edges. The baseline is subtracted
before integration, and only the **positive part** of the corrected signal is
integrated by the trapezoidal rule. When a segment crosses zero it is split at
the exact crossing point, so negative excursions contribute nothing — this is
not equivalent to clamping negative node values to zero first (which would
still count the negative region up to a node as area).

Every sample is assumed to receive the same internal-standard concentration.
The response used for calibration is `analyte_area / internal_standard_area`.
If any internal-standard area is non-positive, or the fitted slope is not
positive, the whole batch is rejected with a non-zero exit code and no reports
are produced.

## Calibration and extrapolation convention

Standards are fitted by ordinary least squares with a free intercept, with
concentration on the x-axis and area ratio on the y-axis. Residuals are
defined as measured minus fitted.

Unknown concentrations are computed as
`(ratio - intercept) / slope`; negative results are reported as-is. An unknown
is flagged `outside_calibration_range` only when its result is strictly below
the minimum or strictly above the maximum standard concentration — landing
exactly on a calibration boundary is not considered extrapolation.

## Tests

```
python -m unittest discover -s tests
```
