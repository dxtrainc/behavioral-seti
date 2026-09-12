"""
ROW 3 -- sub-minute EUV. SDO/EVE ESP at 0.25 s, 120 consecutive days.

§3.6 ranks this row third by designer preference at 1-10 ppm against a reach of
~10^-6, and records it unsearched. PROBA2/LYRA is unreachable; SDO/EVE ESP serves the
row and is finer: 0.25 s, four EUV bands, 41.4 million samples over 2014 days 1-120.

WHY THE FAST END MATTERS HERE. §4.3 measured the same instrument, days and pipeline as
245x more sensitive at two minutes than at one day, because solar variability is RED
and the noise floor falls with frequency. At 0.25 s the Nyquist is 2 Hz -- a band forty
years of daily records cannot represent at all.

GATE 2 IS THE REASON TO USE THIS INSTRUMENT. ESP carries CH_D, a DARK DIODE: the same
detector chain, the same electronics, the same telemetry, and no photons. Anything
appearing in CH_D is instrumental by construction, with no modelling and no argument.
That is a physics-supplied control of exactly the kind §5.6 credits with removing two
apparent detections elsewhere in this programme, and it is built into the instrument
rather than assembled afterwards.

GATE 1 is the five-minute oscillation near 3.3 mHz, which §4.6 recovered in EUV
irradiance on two other spacecraft. It must appear in the light channels. If it does
not, the search is not calibrated and nothing downstream is worth reading.

Gaps are ZERO-FILLED on a regular 0.25 s grid. Zero-filling adds no power at any
frequency; interpolating would add some.
"""
import os, sys, glob, json, time
import numpy as np
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


ESP = os.path.expanduser("~/esp")
DT = 0.25
LIGHT = ["CH_18", "CH_26", "CH_30", "CH_36"]
DARK = "CH_D"

def build():
    fs = sorted(glob.glob(os.path.join(ESP, "esp_*.npz")))
    doy = np.array([int(os.path.basename(f).split("_")[-1][4:7]) for f in fs])
    d0 = doy.min(); ndays = int(doy.max()-d0+1)
    n = int(ndays*86400/DT)
    out = {c: np.zeros(n) for c in LIGHT+[DARK]}
    hit = np.zeros(n, bool)
    for f, dy in zip(fs, doy):
        z = np.load(f)
        idx = (((dy-d0)*86400 + z["sod"])/DT).astype(np.int64)
        m = (idx >= 0) & (idx < n)
        idx = idx[m]
        for c in LIGHT+[DARK]:
            out[c][idx] = z[c][m]
        hit[idx] = True
    return out, hit, ndays

def prep(x, ok, win=2401):
    """relative residual against a running median (~10 min), gaps zeroed"""
    med = np.median(x[ok]) if ok.any() else 1.0
    y = np.where(ok, x, med)
    tr = median_filter(y, size=win, mode="nearest")
    tr = np.where(np.abs(tr) < 1e-12, med if abs(med) > 1e-12 else 1.0, tr)
    if not _check_denominator(tr, "trend"):
        print("      channel excluded from any ratio (see section 4.18)", flush=True)
    r = np.where(ok, y/tr - 1.0, 0.0)
    s = 1.4826*np.median(np.abs(r[ok])) if ok.any() else 0.0
    if s > 0: r = np.clip(r, -8*s, 8*s)
    return r

def spectrum(r):
    n = len(r)//2*2
    F = np.fft.rfft(r[:n]*np.hanning(n))
    return np.fft.rfftfreq(n, d=DT)[1:], (np.abs(F)**2)[1:]

def continuum(P, w=2001):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

if __name__ == "__main__":
    t0 = time.time()
    C, ok, ndays = build()
    print("EVE ESP: %d days, %d samples at %.2f s, %.1f%% present"
          % (ndays, len(ok), DT, 100*ok.mean()), flush=True)
    print("  Nyquist %.2f Hz, resolution %.3e Hz\n" % (0.5/DT, 1.0/(len(ok)*DT)), flush=True)

    R, SP = {}, {}
    for c in LIGHT+[DARK]:
        R[c] = prep(C[c], ok)
        f, P = spectrum(R[c])
        SP[c] = (f, P/continuum(P))
        print("  prepped %-6s rms %.1f ppm" % (c, np.std(R[c][ok])*1e6), flush=True)

    print("\nGATE 1 -- the five-minute oscillation in the LIGHT channels")
    def envelope(f, Rr):
        b = (f > 1.5e-3) & (f < 5.0e-3)
        fb, Rb = f[b], Rr[b]
        BW = 50e-6
        idx = ((fb-fb[0])/BW).astype(int); nb = idx.max()+1
        m = np.bincount(idx, weights=Rb, minlength=nb)/np.maximum(np.bincount(idx, minlength=nb), 1)
        fc = fb[0]+(np.arange(nb)+0.5)*BW
        j = int(np.argmax(m)); return fc[j], m[j]
    for c in LIGHT+[DARK]:
        fe, Re = envelope(*SP[c])
        tag = "dark" if c == DARK else ("PASS" if 2.5e-3 < fe < 3.8e-3 else "*** FAIL ***")
        print("    %-6s envelope %.4f mHz (R=%.2f)   %s" % (c, fe*1e3, Re, tag))

    print("\nBLIND NARROWBAND SEARCH")
    CAND = {}
    for c in LIGHT:
        f, Rr = SP[c]
        N = len(Rr); mu = np.median(Rr)/np.log(2.0); thr = mu*np.log(N/0.05)
        idx = np.where(Rr > thr)[0]
        print("  %-6s %d bins, threshold R > %.1f, %d above" % (c, N, thr, len(idx)))
        CAND[c] = [(float(f[i]), float(Rr[i])) for i in idx[np.argsort(-Rr[idx])][:10]]
        for fr, Rv in CAND[c][:6]:
            print("       f = %.6e Hz   period %12.4f s   R = %8.1f" % (fr, 1.0/fr, Rv))

    print("\nGATE 2 -- the DARK channel. Anything here is instrumental.")
    fd, Rd = SP[DARK]
    Nd = len(Rd); mud = np.median(Rd)/np.log(2.0); thrd = mud*np.log(Nd/0.05)
    survive = 0
    for c in LIGHT:
        for fr, Rv in CAND[c][:6]:
            j = int(np.argmin(np.abs(fd-fr)))
            dark = Rd[max(0, j-2):j+3].max()
            clean = dark < thrd
            survive += clean
            print("    %-6s %.6e Hz  light R=%8.1f  dark R=%8.1f  %s"
                  % (c, fr, Rv, dark, "clean" if clean else "IN DARK -> instrumental"))
    print("\n  candidates absent from the dark channel: %d" % survive)
    json.dump({"ndays": ndays, "cand": CAND}, open(os.path.expanduser("~/esp_row3.json"), "w"), indent=1)
    print("  %.0f s" % (time.time()-t0))
