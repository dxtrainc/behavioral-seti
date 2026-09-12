"""
ROW 7, STEP 3 -- inject through the SELECTED estimator, so the limit is measured.

The bench chose PM1+PM2 mean, 180-day segments, +/-4 uHz centroid, on the KNOWN
solar-cycle term alone: R^2 0.874 against the baseline's 0.517, and a recovered
amplitude of 1.08x the literature 0.4 uHz p-p, which is the check that it measures
frequencies rather than damping variance.

Its implied limit of 1.05e-5 is an EXTRAPOLATION from the baseline's injected figure
and is not a limit. Selection and sensitivity must not share the same evidence: the
estimator was chosen on the activity term, and the limit is now measured by injection
through that estimator, with the choice already locked.

Injection goes into the shift series BEFORE the activity regression, so the signal
meets the same detrending a real one would. Period 2,300 d: incommensurate with the
180-day segment grid and with the year, and clear of the ~11 yr activity term the
regression removes.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G, golf_peakbag as PB, golf_search as SS, golf_bench as B

SEG_D, HW, PER_D, NU_MAX = 180, 4.0, 2300.0, 3.09e-3

if __name__ == "__main__":
    t0 = time.time()
    s1, _ = B.shift_centroid(G.load("PM1"), SEG_D, HW)
    s2, _ = B.shift_centroid(G.load("PM2"), SEG_D, HW)
    sh = np.nanmean(np.vstack([s1, s2]), axis=0)
    act = SS.activity_on_grid(len(sh), SEG_D)
    t = np.arange(len(sh))*SEG_D
    base = np.isfinite(sh) & np.isfinite(act)
    print("selected estimator: PM1+PM2 mean, %d d segments, +/-%.0f uHz" % (SEG_D, HW))
    print("  %d usable segments of %d\n" % (base.sum(), len(sh)), flush=True)

    def resid_of(series):
        m = np.isfinite(series) & np.isfinite(act)
        A = np.vstack([np.ones(m.sum()), act[m]]).T
        co, *_ = np.linalg.lstsq(A, series[m], rcond=None)
        r = np.full(len(series), np.nan); r[m] = series[m] - A@co
        return r

    r0 = resid_of(sh)
    rok = np.isfinite(r0); rv = r0[rok].copy()
    rng = np.random.default_rng(4242)
    NPERM = 5000
    mx = np.empty(NPERM)
    for i in range(NPERM):
        pr = np.full(len(r0), np.nan); pr[rok] = rng.permutation(rv)
        mx[i] = SS.periodogram(pr, rok)[1].max()
    THR = float(np.quantile(mx, 0.95))
    print("  permutation threshold %.4f; residual rms %.4f uHz\n" % (THR, np.nanstd(r0)*1e6), flush=True)

    print("  %10s %13s %10s" % ("amp (uHz)", "fractional", "recovery"))
    out = []
    for amp_u in (0.005, 0.010, 0.015, 0.020, 0.030, 0.045, 0.070, 0.100):
        amp = amp_u*1e-6; rec = 0; NREP = 200
        for k in range(NREP):
            inj = sh + amp*np.sin(2*np.pi*t/PER_D + rng.uniform(0, 2*np.pi))
            r = resid_of(inj); rk = np.isfinite(r)
            if SS.periodogram(r, rk)[1].max() >= THR: rec += 1
        frac = amp/NU_MAX
        out.append((amp_u, frac, rec/NREP))
        print("  %10.3f %13.2e %9.1f%%" % (amp_u, frac, 100*rec/NREP), flush=True)

    hit = [o for o in out if o[2] >= 0.95]
    print()
    if hit:
        a, f, _ = hit[0]
        print("  95%% recovery at %.3f uHz = %.2e fractional at nu_max" % (a, f))
        print("  previous (band estimator, injected):     1.94e-05")
        print("  improvement: %.2fx" % (1.94e-5/f))
        print("  designer level 1e-6  ->  margin %.3f (was 0.05)" % (1e-6/f))
    else:
        print("  95%% recovery not reached within %.3f uHz" % out[-1][0])
    json.dump({"estimator": "PM1+PM2 %dd +/-%.0fuHz" % (SEG_D, HW), "period_d": PER_D,
               "threshold": THR,
               "curve": [{"amp_uHz": a, "frac": f, "rec": r} for a, f, r in out]},
              open(os.path.expanduser("~/golf_inject2.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
