"""Command line entry point: python -m chromatogram INPUT.json --out DIR."""

from __future__ import annotations

import argparse
import sys

from .core import QuantitationError, process_batch
from .inputspec import InputError, load_input
from .report import write_reports


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m chromatogram",
        description="Internal-standard quantitation from chromatogram curves.",
    )
    parser.add_argument("input", help="Path to the batch JSON file")
    parser.add_argument("--out", required=True, help="Directory for CSV/JSON reports")
    args = parser.parse_args(argv)

    try:
        batch = load_input(args.input)
        quant = process_batch(batch)
        write_reports(batch, quant, args.out)
    except (InputError, QuantitationError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
