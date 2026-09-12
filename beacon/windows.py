"""
EVERY WINDOW IS A FILTER, AND ITS TRANSFER FUNCTION IS COMPUTED, NOT JUDGED.

Rule W' in its general form. Four of seven control failures in one session were the
same error -- a smoothing or masking window chosen without checking it against the
scale of the thing it operates on -- and three were made after the first was found:

    GOLF continuum      501 bins  =   0.61 uHz  vs a mode linewidth      -> ACF r = 0.004
    VIRGO envelope     2001 bins  =   2.3  uHz  vs a ~1000 uHz envelope  -> gate failed
    centroid window      +/-8 uHz             across a 9 uHz separation  -> scored 0.87x
    row 5 detrend       1441 min  =   1 day    vs a 6 h transit          -> ate the signal

"501 bins" reads as wide until it is 0.61 uHz against a linewidth. The rule is not a
number: it is that the caller states the frequency it means to KEEP and the frequency
it means to REMOVE, and the transfer function is measured at both.
"""
import numpy as np


def transfer(op, n, fs, freqs):
    """|H(f)| of an arbitrary 1-D filter `op`, measured by probing with unit tones."""
    t = np.arange(n)/float(fs)
    out = []
    for f in np.atleast_1d(freqs):
        x = np.cos(2*np.pi*f*t)
        y = np.asarray(op(x), float)
        r = x - y if _looks_like_detrend(op) else y
        out.append(2.0*np.abs(np.fft.rfft(r*np.hanning(n))).max()/(n/2.0))
    return np.array(out)


def _looks_like_detrend(op):
    return bool(getattr(op, "is_detrend", False))


def checked(op, n, fs, f_signal, f_nuisance, keep=0.9, kill=0.1, label="", budget=None):
    """Run `op` only if it passes the signal it claims to keep and removes the nuisance
    it claims to remove. Otherwise raise -- or, if a `budget` list is supplied, log the
    loss to the sensitivity budget and continue with the loss recorded.

    f_signal    the frequency the analysis is FOR   -- require |H| >= keep
    f_nuisance  the frequency it is removing        -- require |H| <= kill
    """
    hs = float(transfer(op, n, fs, f_signal)[0])
    hn = float(transfer(op, n, fs, f_nuisance)[0])
    ok = (hs >= keep) and (hn <= kill)
    rec = dict(label=label, f_signal=f_signal, H_signal=hs,
               f_nuisance=f_nuisance, H_nuisance=hn, passed=ok)
    if budget is not None:
        budget.append(rec)
        return rec
    if not ok:
        raise ValueError(
            "window %s is the wrong filter: |H| = %.3f at the signal frequency %.4g Hz "
            "(needs >= %.2f) and %.3f at the nuisance frequency %.4g Hz (needs <= %.2f)"
            % (label or repr(op), hs, f_signal, keep, hn, f_nuisance, kill))
    return rec


def running_median(size, mode="nearest"):
    """A detrend as a filter object, so `checked` knows to measure x - median(x)."""
    from scipy.ndimage import median_filter
    def op(x):
        return median_filter(np.asarray(x, float), size=size, mode=mode)
    op.is_detrend = True
    op.size = size
    return op
