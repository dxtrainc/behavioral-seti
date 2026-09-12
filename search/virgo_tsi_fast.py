"""
ROW 6, ON THE RIGHT PRODUCT -- VIRGO TSI at one minute, 27 years, NO highpass.

§4.15 searched VIRGO SPM and reported a null with a restriction: the L2 SPM product
carries a seven-degree polynomial fit and a TWO MONTH HIGHPASS applied upstream, so it
is blind by construction to anything slower. But §3.6 row 6 is "broadband irradiance,
coherent" -- and TSI *is* the broadband product. SPM is the three-colour spectral one.

VIRGO_TSI_Minute carries absolute irradiance in W/m^2 on the WRR scale, 14,199,837
samples at 60 s from 1996 to 2023, 92.8% present, and **no highpass**. It therefore
reaches the band SPM cannot, and it is the more appropriate instrument for this row.

GATE 1 is the same five-minute oscillation, which must come out of a blind search.
GATE 2 has no twin here -- TSI is a single radiometric series -- so a candidate is
checked against the three SPM photometers instead: a real irradiance modulation
appears in the colour channels too, a radiometer artefact does not. That is a
cross-product gate rather than a twin-detector one, and weaker; it is labelled as such.

The instrument lines identified in §4.15 at exactly 180 s and 360 s are excluded here
as well, a priori, being exact multiples of the sample interval.
"""
import os, sys, json, time
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import virgo_fast as V
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

def load_tsi():
    d = np.asarray(fits.open(os.path.join(VD, "VIRGO_TSI_Minute.fits"))[0].data, float)
    return d[:, 0], d[:, 1]

TSI_DETREND = int(os.environ.get("VIRGO_DETREND_MIN", "1441"))

def prep_tsi(v, wmin=None):
    if wmin is None: wmin = TSI_DETREND
    ok = np.isfinite(v) & (v > 1000)
    y = np.where(ok, v, np.nan)
    med = np.nanmedian(y)
    s = 1.4826*np.nanmedian(np.abs(y-med))
    y = np.clip(np.nan_to_num(y, nan=med), med-6*s, med+6*s)
    tr = median_filter(y, size=wmin, mode="nearest")
    _check_denominator(tr, "TSI trend")
    return np.where(ok, y/tr - 1.0, 0.0), ok

if __name__ == "__main__":
    t0 = time.time()
    t, v = load_tsi()
    r, ok = prep_tsi(v)
    print("VIRGO TSI minute: %d samples, %.1f%% present, %.2f yr"
          % (len(v), 100*ok.mean(), (t[-1]-t[0])/3.156e7))
    print("  NO highpass in this product -- reaches the band SPM cannot")
    print("  residual rms %.1f ppm per minute\n" % (r[ok].std()*1e6), flush=True)

    f, P = V.spectrum(r)
    keep = V.line_mask(f)
    C = V.continuum(P)
    Rr = np.where(keep, P/C, 0.0)
    N = int(keep.sum())
    mu = np.median((P/C)[keep])/np.log(2.0)
    thr = mu*np.log(N/0.05)

    print("GATE 1 -- the five-minute oscillation must be recovered")
    b = (f > 2.0e-3) & (f < 4.5e-3) & keep
    fb_, Pb_ = V.rebin(f[b], P[b], 0.5e-6)
    Cb = median_filter(Pb_, size=61, mode="nearest")
    fe, Re = V.rebin(fb_, Pb_/Cb, 50e-6)
    j = int(np.argmax(Re))
    print("    envelope peak %.4f mHz (R=%.2f)   %s\n"
          % (fe[j]*1e3, Re[j], "PASS" if 2.7e-3 < fe[j] < 3.6e-3 else "*** FAIL ***"), flush=True)

    print("BLIND NARROWBAND SEARCH -- %d bins (%.3f%% of band vetoed), threshold R > %.1f"
          % (N, 100*V.masked_fraction(f), thr))
    idx = np.where(Rr > thr)[0]
    cand = [(float(f[i]), float(Rr[i])) for i in idx[np.argsort(-Rr[idx])][:20]]
    if not cand: print("    no bins above threshold")
    for fr, Rv in cand:
        per = 1.0/fr
        tag = ("%.2f min" % (per/60)) if per < 86400 else ("%.2f d" % (per/86400))
        print("    f = %.6e Hz   period %14s   R = %8.1f" % (fr, tag, Rv))

    print("\nGATE 2 -- cross-product: does a candidate appear in the SPM photometers?")
    if cand:
        SPM = {}
        for c in V.CH:
            _, vv = V.load(c)
            rr, _ = V.prep(vv)
            ff, PP = V.spectrum(rr)
            SPM[c] = (ff, PP/V.continuum(PP))
        for fr, Rv in cand[:8]:
            row = []
            for c in V.CH:
                ff, RR = SPM[c]
                row.append(RR[max(0, int(np.argmin(np.abs(ff-fr)))-2):int(np.argmin(np.abs(ff-fr)))+3].max())
            allc = all(x > 10 for x in row)
            print("    %.6e Hz  TSI R=%8.1f   B %7.1f  G %7.1f  R %7.1f  %s"
                  % (fr, Rv, row[0], row[1], row[2], "ALL THREE" if allc else ""))
    else:
        print("    no candidates to check")
    json.dump({"n_bins": N, "threshold": thr, "cand": cand},
              open(os.path.expanduser("~/virgo_tsi_row6.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
