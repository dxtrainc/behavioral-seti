"""
Is the second-half null a real absence, or has the test simply lost power there?

The chip sequence is the SIGN OF THE DAY-TO-DAY CHANGE in F10.7. At solar
maximum that sign tracks active-region evolution and carries information; near
minimum F10.7 is nearly flat and the sign is dominated by measurement noise, so
the code degenerates toward a coin flip. If that is what happens, the second
half is not evidence of absence -- it is a half with no sensitivity, and the
split-half test cannot falsify anything there.

Measured directly: the injection power curve, run separately in each half.
"""
import json, datetime as dt
import numpy as np
from selfkey_replicate import (detrend, chips, corr, load_f107, build,
                               nulls, to_series, ace_daily)
import os

DATA = os.environ.get("BEACON_DATA") or os.path.expanduser("~")   # data root; see README

rng = np.random.default_rng(11)

def power(x, c, depths=(0.0, 0.05, 0.10, 0.25, 0.5), N=400):
    thr = np.percentile(nulls(x, c, "arbitrary"), 95)
    sd = x.std(); n = len(x); out = {}
    for d in depths:
        hits = 0
        for _ in range(N):
            sur = np.roll(x, int(rng.integers(30, n-30)))
            if abs(corr(sur + d*sd*c, c)) > thr: hits += 1
        out[d] = 100.*hits/N
    return out

def keyqual(fk):
    """how much information the chip rule can carry: the size of the daily
    change relative to the flux's own measurement scatter"""
    d = np.diff(fk)
    return float(np.median(np.abs(d))), float(np.std(d))

def main():
    jd_f, flux = load_f107()
    res = {}
    for name, dd in (("ACE", ace_daily()),
                     ("WIND", json.load(open(os.path.join(DATA, "wind_daily.json"))))):
        jd, b = to_series(dd)
        x, c = build(jd, b, jd_f, flux)
        fk = np.interp(jd, jd_f, flux)
        h = len(x)//2; hf = len(fk)//2
        print("\n=== %s ===" % name, flush=True)
        for tag, xs, cs, fs in (("full", x, c, fk),
                                ("first half", x[:h], c[:h], fk[:hf]),
                                ("second half", x[h:], c[h:], fk[hf:])):
            pw = power(xs, cs)
            md, sd = keyqual(fs)
            print("  %-12s n=%4d  F10.7 |dF| median %5.2f sfu   power:  %s"
                  % (tag, len(xs), md,
                     "  ".join("%.2f->%.0f%%" % (k, v) for k, v in pw.items())),
                  flush=True)
            res["%s %s" % (name, tag)] = {"power": pw, "dF_median": md,
                                          "dF_std": sd, "n": len(xs)}
    json.dump(res, open(os.path.join(DATA, "selfkey_power.json"), "w"), indent=1)
    print("\nsaved ~/selfkey_power.json", flush=True)

main()
