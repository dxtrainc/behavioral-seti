"""
ROW 7, THE SEARCH -- now that both gates pass.

golf_row7.py established the pipeline is calibrated: the comb comes out at 135.00 uHz
against 134.9 predicted (r ~ 0.90 on all three channels), and the known solar-cycle
frequency shift comes out at 0.584 uHz peak-to-peak, fractional 1.89e-4 against a
literature ~1.3e-4. A search that can see a known 1e-4 effect has earned the right to
ask about 1e-6.

THE QUESTION. Mode frequencies track solar activity. Is there anything in the shift
series BEYOND that term? The activity dependence is removed by regression -- and the
amount removed is reported, so the reader can see how much of the series was spent.

THE NULL. The shift series is 104 points on a 90-day grid, so the trials count is
small and honest: the searchable band runs from 1/(25.9 yr) to the 180-day Nyquist,
about 52 independent frequencies. Significance is by permutation of the residual,
which is legitimate here because the regression has already removed the one term with
a known time structure.

GATE 2. PM1 and PM2 are two photomultipliers inside one instrument. A real frequency
shift appears in both; a detector artefact does not. Any candidate must survive in
both, and the two series are also compared directly as a consistency check.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G

def activity_on_grid(nseg, seg_days, t0_jd=2450186.5):
    """F10.7 averaged onto the same 90-day segments as the shift series"""
    import datetime as dt
    p = G.os.path.expanduser("~/f107.csv")
    d = np.loadtxt(p, delimiter=",", skiprows=1)
    jd, v = d[:, 0], d[:, 1]
    ok = np.isfinite(v) & (v > 0)
    jd, v = jd[ok], v[ok]
    out = np.full(nseg, np.nan)
    for k in range(nseg):
        a = t0_jd + k*seg_days; b = a + seg_days
        m = (jd >= a) & (jd < b)
        if m.sum() > seg_days*0.5: out[k] = v[m].mean()
    return out

def periodogram(y, ok):
    """unevenly-sampled: Lomb-Scargle over the resolvable band"""
    from scipy.signal import lombscargle
    t = np.where(ok)[0].astype(float)*G.SEG_D          # days
    z = y[ok] - y[ok].mean()
    T = t[-1]-t[0]
    f = np.linspace(1.0/T, 1.0/(2*G.SEG_D), 400)       # cycles/day
    w = 2*np.pi*f
    P = lombscargle(t, z, w, normalize=True)
    return f, P

if __name__ == "__main__":
    z = np.load(os.path.expanduser("~/golf_shift.npz"))
    sh = z["shift"]; seg = int(z["seg_days"])
    print("shift series: %d segments of %d d, %d finite" % (len(sh), seg, np.isfinite(sh).sum()))

    SH = {}
    for tag in ("MEAN", "PM1", "PM2"):
        x = G.load(tag)
        s, duty, _ = G.shift_series(x, int(seg*86400/G.DT))
        SH[tag] = s
        print("  %-5s p-p %.3f uHz" % (tag, (np.nanmax(s)-np.nanmin(s))*1e6), flush=True)

    print("\nGATE 2 -- do the two photomultipliers agree?")
    ok = np.isfinite(SH["PM1"]) & np.isfinite(SH["PM2"])
    r = np.corrcoef(SH["PM1"][ok], SH["PM2"][ok])[0, 1]
    print("    PM1 vs PM2 correlation over %d segments: r = %.3f  %s"
          % (ok.sum(), r, "PASS" if r > 0.5 else "*** FAIL -- not measuring the same thing ***"))

    s = SH["MEAN"]
    act = activity_on_grid(len(s), seg)
    ok = np.isfinite(s) & np.isfinite(act)
    print("\nACTIVITY REGRESSION -- how much of the series is the known term?")
    A = np.vstack([np.ones(ok.sum()), act[ok]]).T
    coef, *_ = np.linalg.lstsq(A, s[ok], rcond=None)
    fit = A @ coef
    resid = np.full(len(s), np.nan); resid[ok] = s[ok] - fit
    ss = 1.0 - np.var(s[ok]-fit)/np.var(s[ok])
    print("    shift = %.3e + %.3e x F10.7      R^2 = %.3f" % (coef[0], coef[1], ss))
    print("    slope: %.3f nHz per sfu; residual rms %.4f uHz"
          % (coef[1]*1e9, np.nanstd(resid)*1e6))

    print("\nSEARCH -- residual, after the activity term")
    rok = np.isfinite(resid)
    f, P = periodogram(resid, rok)
    rng = np.random.default_rng(7)
    NPERM = 20000
    mx = np.empty(NPERM)
    rv = resid[rok].copy()
    for i in range(NPERM):
        pr = np.full(len(resid), np.nan); pr[rok] = rng.permutation(rv)
        mx[i] = periodogram(pr, rok)[1].max()
    thr = np.quantile(mx, 0.95)
    j = int(np.argmax(P))
    p = (1.0 + (mx >= P[j]).sum())/(1.0+NPERM)
    print("    %d independent frequencies, %d permutations" % (len(f)//2, NPERM))
    print("    strongest peak: period %.1f d, power %.4f, threshold %.4f, p = %.4f"
          % (1.0/f[j], P[j], thr, p))
    print("    -> %s" % ("CANDIDATE" if p < 0.05 else "NULL -- no structure beyond the activity term"))
    json.dump({"r_pm1_pm2": float(r), "R2_activity": float(ss),
               "best_period_d": float(1.0/f[j]), "power": float(P[j]),
               "p": float(p), "resid_rms_uHz": float(np.nanstd(resid)*1e6)},
              open(os.path.expanduser("~/golf_row7_search.json"), "w"), indent=1)
