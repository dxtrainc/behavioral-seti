"""
ROW 7, SHARPENED -- per-mode frequency tracking instead of a band cross-correlation.

§4.16 put row 7's limit at 1.94e-5 fractional, twenty times short of the 1e-6 a
designer would set, and said the shortfall was in the METHOD rather than the archive.
This tests that claim instead of asserting it.

WHAT THE FIRST VERSION DID. Cross-correlate the whole 2.5-4.0 mHz band of each 90-day
segment against a mean reference and read off the lag. That treats ~40 independent
oscillation modes as a single rigid template, so every mode contributes equally
regardless of its signal, and the estimate is a bulk shift of the whole comb.

WHAT PEAK-BAGGING DOES. Locate the individual modes once, in the mean spectrum. Then,
in each segment, measure each mode's own centroid in a narrow window and combine the
per-mode shifts weighted by mode power. Three gains, all real:

  - a centroid is limited by linewidth/SNR, not by the 0.129 uHz bin width;
  - strong modes dominate instead of being diluted by weak ones;
  - the scatter ACROSS modes gives a per-segment uncertainty the band method
    could not produce at all, so each point carries an error bar.

Power-weighted centroids are used rather than nonlinear Lorentzian fits. The fits are
the textbook method and would do slightly better, but they fail unpredictably on
low-SNR modes and the failures are hard to detect automatically; the centroid is
stable, and if it already closes most of the gap the extra complexity is not earned.

THE GATE IS UNCHANGED AND MUST STILL PASS: the comb at 135 uHz, the known solar-cycle
shift, and PM1 against PM2. A sharper estimator that loses the controls is not sharper.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G

SEG_D = 90
NU_LO, NU_HI = 2.3e-3, 4.2e-3

def seg_spectra(x, nss):
    """power spectrum of each segment on a common frequency grid"""
    out, duty = [], []
    n = len(x)//nss
    for k in range(n):
        s = x[k*nss:(k+1)*nss]
        duty.append(float(np.mean(s != 0.0)))
        if duty[-1] < 0.5: out.append(None); continue
        f, P = G.spec(s)
        b = (f > NU_LO) & (f < NU_HI)
        out.append(P[b])
    f, _ = G.spec(x[:nss])
    return f[(f > NU_LO) & (f < NU_HI)], out, np.array(duty)

def find_modes(fb, mean_P, dnu=135e-6, min_sep=40e-6):
    """locate modes in the mean spectrum: local maxima above the local continuum"""
    from scipy.ndimage import median_filter
    from scipy.signal import find_peaks
    cont = median_filter(mean_P, size=401, mode="nearest")
    R = mean_P/np.maximum(cont, 1e-300)
    df = fb[1]-fb[0]
    pk, props = find_peaks(R, height=2.0, distance=int(min_sep/df/4))
    order = np.argsort(-props["peak_heights"])
    return pk[order][:60]

def centroid(fb, P, j, half):
    lo, hi = max(0, j-half), min(len(P), j+half+1)
    w = P[lo:hi] - np.median(P[lo:hi])
    w = np.where(w > 0, w, 0.0)
    if w.sum() <= 0: return np.nan, 0.0
    return float((fb[lo:hi]*w).sum()/w.sum()), float(w.sum())

if __name__ == "__main__":
    t0 = time.time()
    nss = int(SEG_D*86400/G.DT)
    X = G.load("MEAN")
    fb, S, duty = seg_spectra(X, nss)
    good = [s for s in S if s is not None]
    mean_P = np.mean(good, axis=0)
    modes = find_modes(fb, mean_P)
    df = fb[1]-fb[0]
    half = max(3, int(round(8e-6/df)))          # +/-8 uHz window per mode
    print("GOLF peak-bagging: %d segments of %d d, %d usable" % (len(S), SEG_D, len(good)))
    print("  bin width %.4f uHz; %d modes located; window +/-%.1f uHz\n"
          % (df*1e6, len(modes), half*df*1e6), flush=True)

    ref = {}
    for j in modes:
        c, w = centroid(fb, mean_P, j, half)
        if np.isfinite(c): ref[j] = c

    shift, err = np.full(len(S), np.nan), np.full(len(S), np.nan)
    for k, P in enumerate(S):
        if P is None: continue
        d, w = [], []
        for j, c0 in ref.items():
            c, ww = centroid(fb, P, j, half)
            if np.isfinite(c) and ww > 0:
                d.append(c-c0); w.append(ww)
        if len(d) < 10: continue
        d, w = np.array(d), np.array(w)
        m = np.abs(d) < 5e-6                      # reject modes that jumped a neighbour
        if m.sum() < 10: continue
        d, w = d[m], w[m]
        shift[k] = float(np.sum(w*d)/np.sum(w))
        err[k] = float(np.std(d)/np.sqrt(len(d)))

    ok = np.isfinite(shift)
    print("  usable segments %d; per-segment uncertainty median %.4f uHz"
          % (ok.sum(), np.nanmedian(err)*1e6))
    print("  shift p-p %.4f uHz (band method gave 0.584)" % ((np.nanmax(shift)-np.nanmin(shift))*1e6))

    import importlib
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    SS = importlib.import_module("golf_search")
    act = SS.activity_on_grid(len(shift), SEG_D)
    m = ok & np.isfinite(act)
    A = np.vstack([np.ones(m.sum()), act[m]]).T
    coef, *_ = np.linalg.lstsq(A, shift[m], rcond=None)
    resid = np.full(len(shift), np.nan); resid[m] = shift[m] - A@coef
    rms = np.nanstd(resid)
    print("\n  activity regression R^2 = %.3f; slope %.3f nHz/sfu"
          % (1.0-np.var(resid[m])/np.var(shift[m]), coef[1]*1e9))
    print("  RESIDUAL RMS %.4f uHz   (band method: 0.0760 uHz)" % (rms*1e6))
    print("  improvement factor %.2fx" % (0.0760/(rms*1e6)))
    print("\n  implied limit ~ %.2e fractional at nu_max (band method 1.94e-5)"
          % (1.94e-5*(rms*1e6)/0.0760))
    np.savez(os.path.expanduser("~/golf_peakbag.npz"), shift=shift, err=err, resid=resid)
    json.dump({"n_modes": len(ref), "resid_rms_uHz": float(rms*1e6),
               "band_rms_uHz": 0.0760, "improvement": float(0.0760/(rms*1e6))},
              open(os.path.expanduser("~/golf_peakbag.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
