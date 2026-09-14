"""
ROW 7, THE ESTIMATOR THE CENTROID VERSION DECLINED TO BUILD.

golf_peakbag.py chose power-weighted centroids over Lorentzian profile fits, and
said why: "the fits ... fail unpredictably on low-SNR modes and the failures are
hard to detect automatically; the centroid is stable, and if it already closes
most of the gap the extra complexity is not earned."

The centroid did NOT close most of the gap. Row 7 stands at 1.46e-5 fractional
against the 1e-6 a designer would set -- fifteen times short -- so the extra
complexity is now earned, and the objection has to be answered rather than
avoided. This module answers it by making the failures detectable:

  1. Every fit is bounded. The centre may not move more than half the fitting
     window, the linewidth is held between 0.3 and 12 uHz (p-modes at nu_max are
     ~1 uHz; anything outside that range is not a mode), and the amplitude must
     be positive.
  2. Every fit must beat its own null. The Lorentzian must reduce chi-square
     against a flat-background fit of the same window by a stated factor, or the
     mode is dropped for that segment.
  3. Every fit is compared with the centroid on the same window. If the two
     disagree by more than the fitted linewidth the fit is treated as failed --
     a profile fit that wanders away from where the power actually is has found
     a neighbouring peak or a noise spike.
  4. A segment needs a quorum of surviving modes or it returns NaN, exactly as
     the centroid version does.

The failure RATE is reported, not hidden. If most fits fail the honest result is
that the method does not work on this data, and that is a finding rather than a
disappointment.

Interface matches golf_bench.shift_centroid so the existing injection and
activity-regression machinery can drive it unchanged.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "search"))

import golf_row7 as G          # noqa: E402
import golf_peakbag as PB      # noqa: E402

# p-mode linewidths at nu_max are of order 1 uHz. These bounds are deliberately
# loose enough to admit real modes and tight enough to reject fits that have run
# away onto the continuum.
GAMMA_MIN, GAMMA_MAX = 0.3e-6, 12.0e-6
CHI2_GAIN = 1.15               # the Lorentzian must beat flat by this factor
MIN_MODES = 8                  # quorum per segment, as in shift_centroid


def _lorentz(nu, nu0, gamma, amp, base):
    return base + amp / (1.0 + (2.0 * (nu - nu0) / gamma) ** 2)


def fit_mode(fb, P, j, half):
    """Fit one Lorentzian in a window. Returns (centre, weight) or (nan, 0).

    Weight is the fitted amplitude, so strong modes dominate the combination in
    the same way power weighting does for the centroid.
    """
    lo, hi = max(0, j - half), min(len(P), j + half + 1)
    if hi - lo < 7:
        return np.nan, 0.0
    x, y = fb[lo:hi], P[lo:hi]
    if not np.all(np.isfinite(y)) or np.ptp(y) <= 0:
        return np.nan, 0.0

    base0 = float(np.median(y))
    amp0 = float(np.max(y) - base0)
    if amp0 <= 0:
        return np.nan, 0.0
    nu0_0 = float(x[int(np.argmax(y))])
    p0 = [nu0_0, 1.0e-6, amp0, base0]

    span = float(x[-1] - x[0])
    bounds = ([nu0_0 - span / 2, GAMMA_MIN, 0.0, 0.0],
              [nu0_0 + span / 2, GAMMA_MAX, amp0 * 20 + 1e-300, base0 * 10 + 1e-300])

    try:
        from scipy.optimize import curve_fit
        popt, _ = curve_fit(_lorentz, x, y, p0=p0, bounds=bounds, maxfev=4000)
    except Exception:
        return np.nan, 0.0

    nu0, gamma, amp, base = popt
    if not np.isfinite(nu0) or amp <= 0:
        return np.nan, 0.0

    # gate 2: must beat a flat fit of the same window
    chi_fit = float(np.sum((y - _lorentz(x, *popt)) ** 2))
    chi_flat = float(np.sum((y - np.mean(y)) ** 2))
    if chi_fit <= 0 or chi_flat / chi_fit < CHI2_GAIN:
        return np.nan, 0.0

    # gate 3: must agree with the centroid to within a linewidth
    c, _ = PB.centroid(fb, P, j, half)
    if np.isfinite(c) and abs(nu0 - c) > gamma:
        return np.nan, 0.0

    return float(nu0), float(amp)


def shift_lorentzian(X, seg_d, half_uHz, report=None):
    """Per-segment frequency shift from Lorentzian profile fits.

    Signature matches golf_bench.shift_centroid. `report` optionally receives
    the fit statistics, because a method whose failure rate is unknown has not
    been validated.
    """
    nss = int(seg_d * 86400 / G.DT)
    fb, S, duty = PB.seg_spectra(X, nss)
    good = [s for s in S if s is not None]
    if len(good) < 8:
        return None, None
    mean_P = np.mean(good, axis=0)
    modes = PB.find_modes(fb, mean_P)
    df = fb[1] - fb[0]
    half = max(3, int(round(half_uHz * 1e-6 / df)))

    # reference centres from the mean spectrum, fitted the same way
    ref = {}
    for j in modes:
        c, w = fit_mode(fb, mean_P, j, half)
        if np.isfinite(c):
            ref[j] = c

    tried = failed = 0
    sh = np.full(len(S), np.nan)
    for k, P in enumerate(S):
        if P is None:
            continue
        d, w = [], []
        for j, c0 in ref.items():
            tried += 1
            c, ww = fit_mode(fb, P, j, half)
            if np.isfinite(c) and ww > 0:
                d.append(c - c0)
                w.append(ww)
            else:
                failed += 1
        if len(d) < MIN_MODES:
            continue
        d, w = np.array(d), np.array(w)
        m = np.abs(d) < 5e-6          # reject modes that jumped to a neighbour
        if m.sum() < MIN_MODES:
            continue
        sh[k] = float(np.sum(w[m] * d[m]) / np.sum(w[m]))

    if report is not None:
        report.update({
            "modes_in_reference": len(ref),
            "modes_located": int(len(modes)),
            "fits_attempted": tried,
            "fits_failed": failed,
            "fit_failure_rate": (failed / tried) if tried else 1.0,
            "segments": len(S),
            "segments_usable": int(np.isfinite(sh).sum()),
        })
    return sh, seg_d


if __name__ == "__main__":
    import json
    import time

    t0 = time.time()
    SEG_D, HW = 180, 4.0
    print("GOLF row 7 -- Lorentzian profile fits vs power-weighted centroids")
    print("  segments %d d, window +/-%.1f uHz\n" % (SEG_D, HW), flush=True)

    import golf_bench as B

    out = {}
    for tag in ("PM1", "PM2"):
        X = G.load(tag)
        rep = {}
        sl, _ = shift_lorentzian(X, SEG_D, HW, report=rep)
        sc, _ = B.shift_centroid(X, SEG_D, HW)
        out[tag] = {"lorentz": sl, "centroid": sc, "report": rep}
        print("  %s: %d/%d modes in reference, fit failure rate %.1f%%, %d/%d segments usable"
              % (tag, rep["modes_in_reference"], rep["modes_located"],
                 100 * rep["fit_failure_rate"], rep["segments_usable"], rep["segments"]),
              flush=True)

    def comb(key):
        a, b = out["PM1"][key], out["PM2"][key]
        return np.nanmean(np.vstack([a, b]), axis=0)

    res = {"seg_d": SEG_D, "half_uHz": HW,
           "reports": {t: out[t]["report"] for t in out}}
    for key in ("centroid", "lorentz"):
        sh = comb(key)
        ok = np.isfinite(sh)
        res[key] = {"usable": int(ok.sum()), "rms_uHz": float(np.nanstd(sh) * 1e6)}
        print("\n  %-9s usable %d  raw rms %.4f uHz" % (key, ok.sum(), np.nanstd(sh) * 1e6))

    np.savez(os.path.expanduser("~/golf_lorentz.npz"),
             lorentz=comb("lorentz"), centroid=comb("centroid"))
    json.dump(res, open(os.path.expanduser("~/golf_lorentz.json"), "w"), indent=1)
    print("\n  saved ~/golf_lorentz.npz and .json   (%.0f s)" % (time.time() - t0))
