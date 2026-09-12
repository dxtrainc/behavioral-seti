"""
ROW 5, GATE 1 REBUILT AS A MATCHED FILTER.

TWO FAILURES PRECEDED THIS, AND BOTH ARE INSTRUCTIVE.

  v1  compared the deepest hour on the transit day against 0.4x the predicted
      depth. Residual rms is 62.4 ppm/min, so noise cleared that threshold and all
      six transits "passed" -- Mercury reading 50-70 ppm against a predicted 12.3.
      A control that cannot fail is not a control.

  v2  added a control-day null and all six FAILED, including Venus at 76 ppm. That
      was the honest answer for the statistic used, but the statistic was wrong:
      "deepest hour anywhere in a day" searches ~1380 positions and rides the red
      noise the daily-median detrend leaves behind.

THE ARITHMETIC SAYS VENUS SHOULD BE EASY. 76 ppm over ~6 h against 62.4 ppm/min white
noise is SNR = 76/(62.4/sqrt(360)) = 23. So a detector that cannot see it is losing
roughly an order of magnitude somewhere, and the two places it can be lost are the
detrend and the filter.

  DETREND. A running median of one day is only 4x the transit duration, so it
  partially absorbs the very dip it should preserve. The window here is 6 days, long
  compared with the transit and short compared with active-region evolution.

  FILTER. A box matched to the KNOWN duration, evaluated at the KNOWN time. There is
  no search, so there is no look-elsewhere penalty -- this is a gate, not a hunt.

THE NULL IS THE SAME FILTER AT THE SAME DURATION ON TRANSIT-FREE DAYS, which measures
whatever red noise survives at six-hour timescales instead of assuming it is white.
If Venus still fails against that, the archive cannot support this row at 10 ppm and
the paper should say so.
"""
import os, sys, json, datetime as dt
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), ".."))
from beacon.controls import denominator_safe

def _check_denominator(tr, name):
    """RULE D. A relative residual x/trend - 1 is only dimensionless while the trend
    stays away from zero. EVE ESP CH_36 crosses zero -- median 2.6e-4, minimum -5.2e-4 --
    so its residual rms came out at 765,522 ppm, 76%, and every ratio formed against it
    was a division by something near zero rather than a carrier. The check is cheap and
    the failure is silent without it."""
    res = denominator_safe(tr, name)
    if not res:
        print("  *** DENOMINATOR UNSAFE: %s -- %s" % (name, res["note"]), flush=True)
    return bool(res)


VD = os.path.expanduser("~/virgo")
TAI0 = dt.datetime(1958, 1, 1)
RSUN = 696000.0

# mid-transit (UT) and duration in hours, from published contact times
TRANSITS = [("Venus",   "2004-06-08 08:20", 6.22, (6052.0/RSUN)**2),
            ("Venus",   "2012-06-06 01:29", 6.67, (6052.0/RSUN)**2),
            ("Mercury", "2003-05-07 07:52", 5.32, (2440.0/RSUN)**2),
            ("Mercury", "2006-11-08 21:41", 4.97, (2440.0/RSUN)**2),
            ("Mercury", "2016-05-09 14:57", 7.50, (2440.0/RSUN)**2),
            ("Mercury", "2019-11-11 15:20", 5.48, (2440.0/RSUN)**2)]

def load():
    d = np.asarray(fits.open(os.path.join(VD, "VIRGO_TSI_Minute.fits"))[0].data, float)
    return d[:, 0], d[:, 1]

def matched(t, v, centre_s, dur_h, detrend_days=6):
    """box matched filter of width dur_h at centre_s; returns depth in ppm.
    Detrend window is long compared with the transit so it cannot absorb it."""
    half = detrend_days*43200.0
    m = (t >= centre_s-half) & (t <= centre_s+half)
    tt, vv = t[m], v[m]
    ok = np.isfinite(vv)
    if ok.sum() < 2000: return np.nan
    w = int(detrend_days*1440) | 1
    base = median_filter(np.where(ok, vv, np.nanmedian(vv[ok])), size=w, mode="nearest")
    _check_denominator(base, "6-day median TSI")
    r = np.where(ok, vv/base - 1.0, np.nan)
    k = int(round(dur_h*60))
    inb = (tt >= centre_s-dur_h*1800) & (tt <= centre_s+dur_h*1800) & ok
    if inb.sum() < k*0.5: return np.nan
    return -np.nanmean(r[inb])*1e6

if __name__ == "__main__":
    t, v = load()
    print("VIRGO TSI minute, matched filter at the known duration and time\n")
    print("  %-9s %-18s %5s %9s %9s %9s %8s %6s"
          % ("body", "mid-transit UT", "dur", "expected", "measured", "null sd", "SNR", "gate"))
    out = []
    for body, when, dur, depth in TRANSITS:
        c = (dt.datetime.strptime(when, "%Y-%m-%d %H:%M") - TAI0).total_seconds()
        meas = matched(t, v, c, dur)
        # null: same filter, same duration, offsets of +/-8..40 days (no transit)
        null = []
        for off in list(range(-40, -7)) + list(range(8, 41)):
            x = matched(t, v, c + off*86400.0, dur)
            if np.isfinite(x): null.append(x)
        if not np.isfinite(meas) or len(null) < 15:
            print("  %-9s %-18s %5.2f %9.1f %9s" % (body, when, dur, depth*1e6, "no data")); continue
        null = np.array(null)
        sd = null.std()
        snr = (meas - null.mean())/sd if sd > 0 else np.nan
        ok = snr > 5
        out.append((body, when, depth*1e6, float(meas), float(sd), float(snr), bool(ok)))
        print("  %-9s %-18s %5.2f %9.1f %9.1f %9.1f %8.1f %6s"
              % (body, when, dur, depth*1e6, meas, sd, snr, "PASS" if ok else "FAIL"))
    n = sum(1 for o in out if o[-1])
    print("\n  %d of %d known transits recovered at SNR > 5" % (n, len(out)))
    if out:
        vv = [o for o in out if o[0] == "Venus"]
        if vv:
            print("  Venus: predicted %.1f ppm, measured %s ppm"
                  % (vv[0][2], " / ".join("%.1f" % o[3] for o in vv)))
    json.dump([{"body": b, "when": w, "expected_ppm": e, "measured_ppm": m,
                "null_sd_ppm": s, "snr": r, "pass": p} for b, w, e, m, s, r, p in out],
              open(os.path.expanduser("~/virgo_row5_matched.json"), "w"), indent=1)
