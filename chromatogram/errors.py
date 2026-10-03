class InputError(Exception):
    """Raised when the input JSON fails validation."""


class AnalysisError(Exception):
    """Raised when the batch cannot be quantified (e.g. non-positive IS area)."""
