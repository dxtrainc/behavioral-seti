"""
ROW 6 INJECTION -- what the VIRGO null is null TO, and where the archive goes blind.

§4.15 reports row 6 as a null but NOT as an injection-verified limit, for a reason
that lives in the archive rather than the analysis: the VIRGO L2 product already has
a seven-degree polynomial fit and a TWO MONTH HIGHPASS applied, plus correction for
orbit, degradation, outliers and "attractors". A limit cannot be quoted until someone
measures what survives that.

This cannot inject before the VDC's processing -- L2 is what exists. What it CAN do,
and what the paper needs, is two things:

  1. the 95% recovery amplitude at periods well inside the searchable band, which is
     the limit proper; and
  2. the recovery as a function of PERIOD across the highpass corner, which turns
     "blind by construction to modulation slower than about two months" from a quoted
     header note into a measured curve.

The second is the point. A caveat stated from documentation is worth less than a
caveat measured, and the blind band is the single largest restriction on this row.

Periods avoid the 180 s and 360 s instrument lines, the 2.5-4.0 mHz p-mode band, the
12 h spacecraft thermal line and exact divisors of a day, after the lesson that a
600 s injection sat on the 144th harmonic of a daily fold and was annihilated by it.
"""
import os, sys, json, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import virgo_fast as V

CH = "BLUE"
import os as _o
_slow = _o.environ.get("VIRGO_DETREND_MIN","1440") != "1440"
PERIODS = ([86413.0, 259207.0, 864013.0, 2592011.0, 5184013.0, 10368017.0] if _slow
           else [307.0, 1009.0, 3607.0, 36011.0, 864013.0, 2592011.0, 5184013.0, 10368017.0])
AMPS    = ([0.10, 0.20, 0.40, 0.80, 1.60, 3.20, 6.40] if _slow
           else [0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.60, 3.20])     # ppm
NREP    = 12

_t, _v = V.load(CH)
T = _t - _t[0]

def once(arg):
    per, amp, seed = arg
    rng = np.random.default_rng(seed)
    v = _v + amp*np.sin(2*np.pi*T/per + rng.uniform(0, 2*np.pi))
    r, ok = V.prep(v)
    f, P = V.spectrum(r)
    keep = V.line_mask(f)
    C = V.continuum(P)
    R = np.where(keep, P/C, 0.0)
    N = int(keep.sum())
    mu = np.median((P/C)[keep])/np.log(2.0)
    thr = mu*np.log(N/0.05)
    j = int(np.argmin(np.abs(f - 1.0/per)))
    lo, hi = max(0, j-3), min(len(R), j+4)
    return per, amp, float(R[lo:hi].max()), float(thr)

if __name__ == "__main__":
    jobs = [(p, a, 900+i) for p in PERIODS for a in AMPS for i in range(NREP)]
    print("VIRGO %s, %d samples, %.2f yr" % (CH, len(_v), (_t[-1]-_t[0])/3.156e7))
    print("L2 highpass corner ~2 months = %.3e Hz" % (1.0/(60*86400)))
    print("%d jobs (%d periods x %d amplitudes x %d reps)\n" % (len(jobs), len(PERIODS), len(AMPS), NREP))
    t0 = time.time()
    with Pool(20) as pool: out = pool.map(once, jobs, chunksize=1)
    print("%.0f s\n" % (time.time()-t0))

    print("recovery fraction (R above the search's own threshold)")
    hdr = "  %-14s" % "period"
    for a in AMPS: hdr += "%7.2f" % a
    print(hdr + "   ppm")
    lim = {}
    for per in PERIODS:
        rows = [o for o in out if o[0] == per]
        line = "  %-14s" % (("%.0f s" % per) if per < 86400 else ("%.1f d" % (per/86400)))
        for a in AMPS:
            rr = [o for o in rows if o[1] == a]
            rec = np.mean([o[2] >= o[3] for o in rr])
            line += "%6.0f%%" % (100*rec)
            if rec >= 0.95 and per not in lim: lim[per] = a
        print(line)

    print("\n95%% recovery amplitude by period:")
    for per in PERIODS:
        tag = ("%.0f s" % per) if per < 86400 else ("%.1f d" % (per/86400))
        print("  %-12s %s" % (tag, ("%.2f ppm" % lim[per]) if per in lim
                              else "NOT REACHED in %.2f ppm" % AMPS[-1]))
    inband = [lim[p] for p in PERIODS if p in lim]
    if inband:
        print("\n  best recovered amplitude across the scanned band: %.2f ppm = %.2e fractional"
              % (min(inband), min(inband)*1e-6))
    json.dump({"channel": CH, "periods": PERIODS, "amps": AMPS,
               "limit_ppm": {str(k): v for k, v in lim.items()}},
              open(os.path.expanduser("~/virgo_row6_injection.json"), "w"), indent=1)
