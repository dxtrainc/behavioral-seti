"""
THE ONLY PLACE SURROGATE CONSTRUCTION IS ALLOWED TO LIVE.

Before this module there were four copies of the same phase screen -- two on the
compressed index axis (wl_sweep.common_phase, wl_twin.screen) and two on the full grid
(wl_cal3.full_phase, wl_cal4.screen) -- and two copies of the masked roll. The
compressed-axis form was measured at 9.39x over-firing on a gapped daily series; the
fix was written in wl_cal4 and never reached wl_twin. A fix that lands in one copy is
not a fix. Every search script imports from here; none defines its own.
"""
import numpy as np


def phase_screen(Y, rng, mask=None):
    """Common phase screen across channels: preserves every pairwise cross-spectrum
    exactly and destroys higher-order alignment.

    Y is (channels, N) on the FULL time grid. If `mask` is given (bool, length N) the
    result is subsampled through it, so data and surrogate share a calendar and not
    merely a spectrum.

    Refuses a 1-D input: the caller must be explicit about the channel axis.
    """
    Y = np.asarray(Y, float)
    if Y.ndim != 2:
        raise ValueError("phase_screen expects (channels, N) on the full grid; "
                         "got shape %r" % (Y.shape,))
    n = Y.shape[1]
    F = np.fft.rfft(Y, axis=1)
    ph = rng.uniform(0, 2*np.pi, F.shape[1])
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    S = np.fft.irfft(F*np.exp(1j*ph)[None, :], n=n, axis=1)
    return S if mask is None else S[:, mask]


def roll_masked(y, ok, rng, margin=50):
    """Circular shift on the full time grid with the real mask re-applied.

    The mask must be re-applied AFTER the roll: rolling moves the gaps, and a gap that
    has moved is a stretch of real data in the wrong place.
    """
    y = np.asarray(y, float)
    s = np.roll(y, int(rng.integers(margin, len(y)-margin))).copy()
    s[~ok] = 0.0
    return s


def fill_to_grid(values, finite_mask):
    """Interpolate a gapped series onto its full grid FOR SURROGATE GENERATION ONLY.

    The interpolated values must never be used as data -- that is what the mask is for.
    viewpoint_search.py interpolated MAVEN's missing 5.7% and then used the result in
    the spectrum, the continuum, the threshold and the injections; interpolated
    stretches are smoother than data, which lowers apparent noise.
    """
    idx = np.arange(len(values), dtype=float)
    return np.interp(idx, idx[finite_mask], np.asarray(values, float)[finite_mask])


def gap_report(mask, shortest_period_samples):
    """Rule S': index-axis surrogates are permitted only when the gaps cannot matter.

    Returns index_axis_ok=True only if the gap fraction is under 2% AND the longest gap
    is shorter than the shortest period searched. On the daily sweep (10.9% gaps) this
    is False, which is the condition under which the compressed-axis screen over-fired
    by 9.39x; on 1-minute GOES it is typically True.
    """
    mask = np.asarray(mask, bool)
    idx = np.where(mask)[0]
    gaps = np.diff(idx) - 1 if len(idx) > 1 else np.array([0])
    longest = int(gaps.max()) if len(gaps) else 0
    frac = 1.0 - float(mask.mean())
    return dict(gap_fraction=frac,
                longest_gap=longest,
                n_breaks=int((gaps > 0).sum()),
                index_axis_ok=bool(frac < 0.02 and longest < shortest_period_samples))
