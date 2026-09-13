"""
DETERMINISTIC SEEDS. Python's hash() is salted per process.

Five sites derived a surrogate seed from hash((i, j, form)):

    search/sweepT.py:327 and :349     the pair and triple sweeps
    search/sweep13.py:242
    search/sweep30.py:266
    search/confirm.py:321

PYTHONHASHSEED is randomised at interpreter start unless it is set, so hash() on a tuple
containing a string returns a DIFFERENT value in every process. Three consecutive runs
gave 1689220225, 56146563 and 1311111642 for the same key. Every p-value produced through
those seeds is therefore unreproducible: the surrogates cannot be regenerated, so the
numbers cannot be checked by anyone, including us.

This is not a correctness problem -- a seed drawn arbitrarily is still a valid seed, and
the p-values were computed from genuine surrogate draws. It is a REPRODUCIBILITY problem,
and for a paper whose central claim is that a null result is worth the fraction of a space
it excludes, an unreproducible null is worth less than it looks.

stable_seed() uses BLAKE2b over a canonical repr, which is fixed across processes,
interpreter versions and platforms.

NOTE ON THE COMMITTED NUMBERS. Switching to a stable seed changes which surrogates are
drawn, so a re-run produces different p-values from those in the paper. Nothing is lost by
that: the committed values were already unreproducible, and could not have been
regenerated under the old scheme either. What changes is that everything from here is
checkable.
"""
import hashlib
import numpy as np


def stable_seed(*key):
    """A 32-bit seed that is identical in every process, forever.

    >>> stable_seed(0, 1, "ratio") == stable_seed(0, 1, "ratio")
    True
    """
    canon = "|".join(repr(k) for k in key).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(canon, digest_size=4).digest(), "big")


def rng_for(*key):
    """The generator a search should use for one test, reproducibly."""
    return np.random.default_rng(stable_seed(*key))
