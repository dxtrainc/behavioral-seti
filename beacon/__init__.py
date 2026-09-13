"""Shared analysis library. Search scripts import from here and define no surrogate,
window, injection or control of their own -- see docs/beacon-methodology-review.md."""
from . import surrogates, windows, inject, controls, seeds   # noqa: F401
