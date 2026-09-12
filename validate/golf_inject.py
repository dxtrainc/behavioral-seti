"""
ROW 7 INJECTION -- turning a null into a limit.

The paper's standard is the amplitude recovered 95% of the time under injection,
measured against the threshold the search itself uses. Row 7 is null; this measures
what it is null TO.

The signal is injected into the FREQUENCY SHIFT SERIES, which is where a modulation
of mode frequencies would appear -- i.e. before the activity regression, so the
injected signal has to survive the same detrending a real one would. Injecting after
the regression would exempt it from a filter it must actually pass.

Period 1,100 d: incommensurate with the 90-day segment grid, with the year, and well
clear of the 11-year activity term the regression removes.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G, golf_search as S

NU_MAX = 3.09e-3
PER_D = 1100.0

if __name__ == "__main__":
    z = np.load(os.path.expanduser("~/golf_shift.npz"))
    sh = z["shift"]; seg = int(z["seg_days"])
    act = S.activity_on_grid(len(sh), seg)
    base_ok = np.isfinite(sh) & np.isfinite(act)
    t = np.arange(len(sh))*seg

    rng = np.random.default_rng(31)
    # threshold from the permutation null of the REAL residual, once
    A = np.vstack([np.ones(base_ok.sum()), act[base_ok]]).T
    coef, *_ = np.linalg.lstsq(A, sh[base_ok], rcond=None)
    resid0 = np.full(len(sh), np.nan); resid0[base_ok] = sh[base_ok] - A@coef
    rok = np.isfinite(resid0)
    rv = resid0[rok].copy()
    NPERM = 5000
    mx = np.empty(NPERM)
    for i in range(NPERM):
        pr = np.full(len(resid0), np.nan); pr[rok] = rng.permutation(rv)
        mx[i] = S.periodogram(pr, rok)[1].max()
    THR = float(np.quantile(mx, 0.95))
    print("permutation threshold (alpha=0.05): power %.4f" % THR)
    print("residual rms %.4f uHz\n" % (np.nanstd(resid0)*1e6))

    print("  %10s %12s %10s" % ("amp (uHz)", "fractional", "recovery"))
    NREP = 200
    out = []
    for amp_u in (0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.24, 0.32):
        amp = amp_u*1e-6
        rec = 0
        for k in range(NREP):
            ph = rng.uniform(0, 2*np.pi)
            inj = sh + amp*np.sin(2*np.pi*t/PER_D + ph)     # injected BEFORE regression
            ok = np.isfinite(inj) & np.isfinite(act)
            Ai = np.vstack([np.ones(ok.sum()), act[ok]]).T
            ci, *_ = np.linalg.lstsq(Ai, inj[ok], rcond=None)
            r = np.full(len(sh), np.nan); r[ok] = inj[ok] - Ai@ci
            rk = np.isfinite(r)
            if S.periodogram(r, rk)[1].max() >= THR: rec += 1
        frac = amp/NU_MAX
        out.append((amp_u, frac, rec/NREP))
        print("  %10.2f %12.2e %9.1f%%" % (amp_u, frac, 100*rec/NREP))

    hit = [o for o in out if o[2] >= 0.95]
    if hit:
        a, f, r = hit[0]
        print("\n  95%% recovery at %.2f uHz = %.2e fractional at nu_max" % (a, f))
        print("  designer level for row 7 (section 3.6): 1e-6   ->  margin %.2f" % (1e-6/f))
    else:
        print("\n  95%% recovery not reached in the scanned range")
    json.dump({"threshold": THR, "period_d": PER_D,
               "curve": [{"amp_uHz": a, "frac": f, "recovery": r} for a, f, r in out]},
              open(os.path.expanduser("~/golf_row7_injection.json"), "w"), indent=1)
