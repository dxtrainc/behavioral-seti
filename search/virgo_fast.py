"""
ROW 6 -- SOHO/VIRGO SPM, broadband irradiance, coherent, at 0.1 ppm.

Section 3.6 orders the unsearched rows 5, 6, 3, 7 and notes that rows 5 and 6 are
the same archive twice. This is row 6: a blind coherent narrowband search of the
three VIRGO sun-photometer channels over 27.3 years at 60 s cadence.

THE ARCHIVE WAS NOT ACTUALLY BLOCKED. The programme recorded rows 5-6 as unreachable
because ftp.pmodwrc.ch publishes only an IPv6 AAAA record and no host here has IPv6
egress. That is still true of PMOD's FTP -- three variants were retried and all fail
-- but NASA mirrors the entire mission over HTTPS at soho.nascom.nasa.gov. The
blocker was real; the conclusion drawn from it was not.

GATE 1, AND WHY IT IS UNUSUALLY STRONG HERE. A blind search that cannot recover a
known real line is not calibrated. For VIRGO the control is the best available
anywhere: VIRGO SPM is the instrument solar p-modes are classically measured with,
so the five-minute oscillation -- an envelope peaking near 3.09 mHz with modes spaced
by 135 uHz -- MUST come out of this search. If it does not, nothing else the search
says is worth reading.

GATE 2. BLUE, GREEN and RED are three independent photometers behind three
independent filters. A solar signal appears in all three; a detector artefact does
not. This is the twin-instrument logic of fastsearch.py with three witnesses.

DIMENSIONLESS CARRIERS. The data is already a relative signal in ppm, so
log(I/I0) = ppm/1e6 to first order and the DIFFERENCE of two channels in ppm IS
their log colour ratio. BLUE-GREEN, BLUE-RED and GREEN-RED are therefore
dimensionless by construction -- the class a sender is restricted to -- and cost
nothing extra to form.

THE BLIND BAND, WHICH MUST BE QUOTED WITH ANY LIMIT. The L2 header states a
seven-degree polynomial fit and a TWO MONTH HIGHPASS FILTER, plus correction for
orbit, degradation, outliers and "attractors". This product is blind by construction
to modulation slower than ~2 months, and the outlier correction may remove exactly
the impulsive or quantised structure a search like this looks for. The first is a
hard limit. The second is testable by injection and is NOT assumed here.
"""
import os, sys, time, json
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

VDIR = os.path.expanduser("~/virgo")
CH = {"BLUE": "VIRGO-SPM-BLUE-L2-MISSIONLONG.fits",
      "GREEN": "VIRGO-SPM-GREEN-L2-MISSIONLONG.fits",
      "RED": "VIRGO-SPM-RED-L2-MISSIONLONG.fits"}

def load(name):
    d = fits.open(os.path.join(VDIR, CH[name]))[0].data
    return d[:, 0].astype(float), d[:, 1].astype(float)      # seconds TAI, ppm

def prep(v, wmin=1440):
    """gaps zeroed (adds no power at any frequency; interpolating would add some),
    outliers clipped at the robust level, residual slow trend removed."""
    ok = np.isfinite(v)
    y = np.where(ok, v, np.nan)
    med = np.nanmedian(y); s = 1.4826*np.nanmedian(np.abs(y-med))
    y = np.clip(y, med-6*s, med+6*s)
    m = ok.astype(np.float64); z = np.nan_to_num(y)
    cv = np.concatenate([[0], np.cumsum(z)]); cm = np.concatenate([[0], np.cumsum(m)])
    h = wmin//2; i = np.arange(len(y))
    lo = np.maximum(0, i-h); hi = np.minimum(len(y), i+h)
    den = cm[hi]-cm[lo]
    trend = np.where(den > 10, (cv[hi]-cv[lo])/np.maximum(den, 1), med)
    return np.where(ok, y-trend, 0.0), ok

def spectrum(r, dt=60.0):
    n = len(r)//2*2
    F = np.fft.rfft(r[:n]*np.hanning(n))
    return np.fft.rfftfreq(n, d=dt)[1:], (np.abs(F)**2)[1:]


# ---------------------------------------------------------------- INSTRUMENT LINES
# The first run returned 18 "gate 2 candidates" and every one sat at 5.555556e-3 Hz
# (180.0000 s) or 2.777778e-3 Hz (360.0000 s) -- EXACTLY three and six minutes, i.e.
# exactly 2/3 and 1/3 of the 60 s Nyquist. A solar signal has no reason to land on a
# round number in OUR units; these are a resampling or telemetry artefact. They reach
# R = 19,486, thousands of times the p-mode envelope, and they dragged the gate-1
# envelope estimator to 3.86 and 4.50 mHz on GREEN and RED, failing the gate.
#
# They are identified A PRIORI -- exact multiples of the sample interval -- and not
# chosen after seeing which bins came out large, and they are EXCLUDED rather than
# subtracted: the search does not run there and the blind band is reported.
LINES = [1.0/180.0, 1.0/360.0]
def line_mask(f, halfwidth_rel=3e-4):
    """True where the search is allowed to look"""
    m = np.ones(len(f), bool)
    for L in LINES:
        for h in (1, 2, 3):
            m &= np.abs(f - L*h) > L*h*halfwidth_rel
    return m



def rebin(f, P, width):
    """rebin a power spectrum to `width` Hz per bin, by mean"""
    idx = ((f - f[0])/width).astype(np.int64)
    nb = int(idx.max())+1
    cnt = np.bincount(idx, minlength=nb).astype(float)
    val = np.bincount(idx, weights=P, minlength=nb)
    ok = cnt > 0
    return f[0] + (np.arange(nb)[ok]+0.5)*width, val[ok]/cnt[ok]


def continuum(P, w=801):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

def search(P, f, label, alpha=0.05, top=15):
    keep = line_mask(f)
    C = continuum(P); R = np.where(keep, P/C, 0.0); N = int(keep.sum())
    mu = np.median((P/C)[keep])/np.log(2.0)
    thr = mu*np.log(N/alpha)
    idx = np.where(R > thr)[0]
    print("  %-14s %d bins, threshold R > %.1f, %d above" % (label, N, thr, len(idx)))
    out = []
    for i in idx[np.argsort(-R[idx])][:top]:
        out.append((float(f[i]), float(R[i])))
        print("       f = %.6e Hz   period %11.4f min   R = %8.1f"
              % (f[i], 1.0/f[i]/60, R[i]))
    if not len(idx): print("       no bins above threshold")
    return out, R, f, mu, thr

if __name__ == "__main__":
    t0 = time.time()
    T, V = {}, {}
    for c in CH:
        t, v = load(c); T[c] = t; V[c] = v
        print("%-6s %d samples, %.1f%% finite, %.2f yr"
              % (c, len(v), 100*np.isfinite(v).mean(), (t[-1]-t[0])/3.156e7))
    assert all(np.array_equal(T["BLUE"], T[c]) for c in CH), "time axes differ"

    R_, OK = {}, {}
    for c in CH:
        R_[c], OK[c] = prep(V[c])
    print("\nprepared in %.1f s\n" % (time.time()-t0))

    print("GATE 1 -- the five-minute oscillation MUST be recovered")
    print("  expect an envelope near 3.09 mHz (period 5.4 min)")
    for c in CH:
        f, P = spectrum(R_[c])
        # THE FIRST TWO RUNS FAILED THIS GATE ON A WINDOW SIZE. VIRGO runs 14,342,400
        # samples at 60 s, so df = 1.16 nHz and a median_filter of 2001 bins is
        # 2.3 uHz -- against a p-mode envelope roughly 1000 uHz wide. It was
        # measuring noise, and it railed at the band edge on RED. The envelope is a
        # BROAD feature: flatten against a continuum wide compared with a mode, then
        # rebin to 50 uHz and take the peak of that.
        band = (f > 2.0e-3) & (f < 4.5e-3) & line_mask(f)
        fb_, Pb_ = rebin(f[band], P[band], 0.5e-6)
        Cb = median_filter(Pb_, size=61, mode="nearest")
        fe, Re = rebin(fb_, Pb_/Cb, 50e-6)
        j = int(np.argmax(Re))
        i = int(np.argmax(Pb_/Cb))
        print("    %-6s band peak %.4f mHz (R=%.1f)   envelope peak %.4f mHz (R=%.2f)  %s"
              % (c, fb_[i]*1e3, (Pb_/Cb)[i], fe[j]*1e3, Re[j],
                 "PASS" if 2.7e-3 < fe[j] < 3.6e-3 else "*** FAIL ***"))

    print("\nBLIND NARROWBAND SEARCH -- single channels (row 6)")
    CAND = {}
    SPEC = {}
    for c in CH:
        f, P = spectrum(R_[c])
        CAND[c], Rv, fv, mu, thr = search(P, f, c)
        SPEC[c] = (fv, Rv)

    print("\nDIMENSIONLESS COLOUR RATIOS -- ppm differences are log ratios")
    for a, b in (("BLUE", "GREEN"), ("BLUE", "RED"), ("GREEN", "RED")):
        d = R_[a] - R_[b]
        f, P = spectrum(d)
        CAND[a+"-"+b], Rv, fv, _, _ = search(P, f, a+"-"+b)
        SPEC[a+"-"+b] = (fv, Rv)

    print("\nGATE 2 -- does any single-channel candidate appear in all three?")
    hits = 0
    for c in CH:
        for fr, Rv in CAND[c][:8]:
            row = []
            for c2 in CH:
                fv, Rr = SPEC[c2]
                row.append(Rr[np.argmin(np.abs(fv-fr))])
            allthree = all(x > 10 for x in row)
            hits += allthree
            print("    %-6s %.6e Hz (%9.3f min)  B %8.1f  G %8.1f  R %8.1f  %s"
                  % (c, fr, 1/fr/60, row[0], row[1], row[2],
                     "ALL THREE" if allthree else ""))
    print("\n  candidates passing gate 2: %d" % hits)
    json.dump({k: v for k, v in CAND.items()},
              open(os.path.expanduser("~/virgo_row6.json"), "w"), indent=1)
    print("\ntotal %.1f s" % (time.time()-t0))
