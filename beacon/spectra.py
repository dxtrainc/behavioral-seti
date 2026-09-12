"""Spectrum, continuum and threshold. Three copies of these existed across the search
scripts; the threshold in particular is the quantity that sets every fast-band limit,
and section 4.2 already measured its exponential tail model as optimistic by ~4e3 on
GOES. It is therefore stated here with its assumption attached."""
import numpy as np
from scipy.ndimage import median_filter


def spec(r, dt=1.0, window=np.hanning):
    n = len(r)//2*2
    F = np.fft.rfft(r[:n]*window(n))
    return np.fft.rfftfreq(n, d=dt)[1:], (np.abs(F)**2)[1:], float(window(n).sum())


def cont(P, w=801):
    """Local median of log power: the red-noise floor, robust to the lines on it.
    `w` is a WINDOW and therefore a filter -- it must be wide compared with a line and
    narrow compared with the continuum's own curvature. Pass it through
    windows.checked if either scale is in doubt."""
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))


def thr_of(R, alpha=0.05):
    """Family-wise threshold on R = P/continuum, assuming R is exponential.

    THE ASSUMPTION IS THE EXPOSURE. Periodogram power over a smooth continuum is
    exponential only if the continuum is right and the noise is Gaussian; section 4.2
    measured the real tail on GOES as exceeding this model by about 4e3. Any limit set
    by this threshold should carry a measured band-wide false-alarm rate beside it --
    controls.negative_arm over ALL judged bins, not one.
    """
    mu = np.median(R)/np.log(2.0)
    return mu*np.log(len(R)/alpha), mu
