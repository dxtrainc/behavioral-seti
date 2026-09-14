"""
THE SIDEREAL FOLD, INJECTED AGAINST ITS OWN GATE.

This is the one row in the paper reported as an analytic sensitivity rather than
an injection-verified limit, and section 4.2 says why: the detector fires on
unmodified data. The seasonal leakage of section 4.4 is real -- sidereal sits at
4.9 sigma and ANTI-sidereal at 6.6 sigma, at a frequency where nothing physical
lives -- so a random-phase injection into the fold cancels as often as it adds,
and the recovery curve measures the leakage rather than the injection.

WHAT WAS WRONG WITH THE OLD INJECTION. inj2more.py injects at the sidereal
frequency and thresholds on amp(sidereal). That statistic is already above
threshold with nothing injected. A test whose null state is "detected" cannot
measure a detection.

THE FIX IS TO INJECT AGAINST THE GATE, NOT THE FOLD. The control that killed the
original candidate is the anti-sidereal line, and the reason it works is that the
two lines respond differently:

    a real sidereal signal  ->  raises sidereal ONLY
    seasonal leakage        ->  raises sidereal AND anti-sidereal together

so the contrast D = amp(sidereal) - amp(anti-sidereal) is blind to the leakage
and sensitive to the signal. Thresholding on D rather than on amp(sidereal)
gives a statistic whose null state is "not detected", which is the precondition
for measuring anything.

Frequencies, in cycles per year:
    solar         365.2422
    sidereal      366.2422        = solar + 1
    anti-sidereal 364.2422        = solar - 1, the mirror of sidereal
The anti-sidereal line is the standard control in cosmic-ray anisotropy work for
exactly this reason: nothing physical produces it, so whatever appears there is
instrumental or seasonal, and it appears at sidereal in equal measure.

The first amplitude tested is ZERO. The previous row-7 injection reported 100%
recovery at its smallest amplitude because it had no zero arm and was measuring
its own residual; that mistake is not repeated here.
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import median_filter

YR = 365.2422 * 86400.0
F_SOLAR = 365.2422 / YR
F_SID = 366.2422 / YR
F_ANTI = 364.2422 / YR
NOISE_BAND = np.linspace(358.0, 363.0, 26)     # cyc/yr, away from all three lines
N_TRIAL = 200
STEP = 5                                        # decimation, as in inj2more.py


def _p(n):
    for c in (n, os.path.join(os.path.expanduser("~"), os.path.basename(n))):
        if os.path.exists(c):
            return c
    return n


def load_series():
    z = np.load(_p("nm_OULU_60s.npz"))
    w = z["v"].astype(float)
    m = np.isfinite(w)
    w = np.where(m, w, np.nanmedian(w[m]))
    w = w / np.median(w) - 1.0
    w = w - median_filter(w, size=1441, mode="nearest")
    t = np.arange(0, len(w), STEP) * 60.0
    return t, w[::STEP]


def amp(x, t, f):
    ph = 2 * np.pi * f * t
    return 2.0 * np.hypot(np.dot(x, np.cos(ph)), np.dot(x, np.sin(ph))) / len(x)


def contrast(x, t):
    """The gated statistic: sidereal power above the anti-sidereal mirror."""
    return amp(x, t, F_SID) - amp(x, t, F_ANTI)


def main():
    t, x = load_series()
    print("sidereal fold, injected against the anti-sidereal gate")
    print("  %d samples at %d s\n" % (len(x), STEP * 60))

    a_sid = amp(x, t, F_SID)
    a_anti = amp(x, t, F_ANTI)
    a_sol = amp(x, t, F_SOLAR)
    noise = np.array([amp(x, t, c / YR) for c in NOISE_BAND])
    sd = float(np.std(noise))

    print("  unmodified data:")
    print("    solar         %.3e   (%.1f sd)" % (a_sol, a_sol / sd))
    print("    sidereal      %.3e   (%.1f sd)" % (a_sid, a_sid / sd))
    print("    anti-sidereal %.3e   (%.1f sd)" % (a_anti, a_anti / sd))
    print("    -> the fold statistic amp(sidereal) is ABOVE a 5 sd threshold "
          "with nothing injected" if a_sid > 5 * sd else "")

    # null distribution of the contrast, from pairs in the noise band
    pairs = []
    for i in range(0, len(NOISE_BAND) - 1, 2):
        f1, f2 = NOISE_BAND[i] / YR, NOISE_BAND[i + 1] / YR
        pairs.append(amp(x, t, f1) - amp(x, t, f2))
    sd_d = float(np.std(pairs))
    thr = 5.0 * sd_d
    d0 = contrast(x, t)
    print("\n  gated statistic D = amp(sidereal) - amp(anti-sidereal):")
    print("    null sd from %d noise-band pairs: %.3e" % (len(pairs), sd_d))
    print("    threshold 5 sd = %+.3e" % thr)
    print("    observed D on unmodified data = %+.3e  (%.1f sd)  -> %s"
          % (d0, d0 / sd_d, "FIRES" if d0 > thr else "does not fire, as required"))

    rng = np.random.default_rng(11)
    amps = np.array([0.0, 5e-5, 1e-4, 1.8e-4, 3e-4, 5e-4, 1e-3, 2e-3])
    rec, curve = [], []
    print("\n  injection scan, first amplitude ZERO:")
    for a in amps:
        hit = 0
        for _ in range(N_TRIAL):
            ph0 = rng.uniform(0, 2 * np.pi)
            xi = x + a * np.cos(2 * np.pi * F_SID * t + ph0)
            if contrast(xi, t) > thr:
                hit += 1
        r = hit / N_TRIAL
        rec.append(r)
        curve.append({"amp": float(a), "rec": r})
        print("    amp %8.1e -> %5.1f%%" % (a, 100 * r), flush=True)

    mono = all(rec[i] >= rec[i - 1] - 0.05 for i in range(1, len(rec)))

    def at(p):
        for i in range(1, len(rec)):
            if rec[i] >= p and rec[i - 1] < p:
                u = (p - rec[i - 1]) / max(rec[i] - rec[i - 1], 1e-9)
                return float(np.exp(np.log(max(amps[i - 1], 1e-12))
                                    + u * (np.log(amps[i]) - np.log(max(amps[i - 1], 1e-12)))))
        return None

    a95 = at(0.95)
    print("\n  zero-amplitude recovery: %.1f%%  %s"
          % (100 * rec[0], "(correct: the gate does not fire on nothing)"
             if rec[0] <= 0.10 else "*** STILL FIRES ON NOTHING ***"))
    print("  monotonic: %s" % mono)
    print("  95%% recovery at %s   (paper quotes 1.8e-4 as an analytic sensitivity)"
          % ("%.2e" % a95 if a95 else "above the grid"))

    out = {
        "statistic": "amp(sidereal) - amp(anti-sidereal)",
        "unmodified": {"solar": a_sol, "sidereal": a_sid, "anti_sidereal": a_anti,
                       "noise_sd": sd, "sidereal_sd": a_sid / sd,
                       "anti_sidereal_sd": a_anti / sd},
        "gate": {"null_sd": sd_d, "threshold": thr, "observed": d0,
                 "fires_on_unmodified": bool(d0 > thr)},
        "curve": curve, "zero_arm": rec[0], "monotonic": bool(mono),
        "a95": a95, "quoted_analytic": 1.8e-4,
    }
    json.dump(out, open(os.path.expanduser("~/sidereal_gate.json"), "w"), indent=1)
    print("\n  saved ~/sidereal_gate.json")


if __name__ == "__main__":
    sys.exit(main())
