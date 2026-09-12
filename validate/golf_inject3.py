"""
ROW 7, STEP 3 DONE PROPERLY -- with the zero-amplitude control that was missing.

TWO THINGS WENT WRONG IN THE PREVIOUS ATTEMPT AND BOTH ARE RECORDED HERE.

  1. NO ZERO-AMPLITUDE ARM. The scan began at 0.005 uHz and reported 100% recovery
     there -- a 12x improvement and margin 0.618. It was spurious: with no injection
     at all the test still "recovered" 100% of the time, because the residual's own
     peak already sat above the permutation threshold. The curve was measuring the
     residual. A test with no arm that can fail is not a test, and this is the second
     time today that exact omission produced a headline number.

  2. THE RESIDUAL PEAK WAS UNMODELLED SOLAR PHYSICS. It sat at 1811 d on 180-day
     segments and 3694 d on 90-day ones -- half the solar cycle and the solar cycle.
     The activity regression was LINEAR in F10.7, but the frequency/activity relation
     has curvature and hysteresis, so a linear fit leaves power at exactly those
     periods. Adding a quadratic term drops the peak from the 99.8th percentile to the
     50.6th; adding dA/dt as well takes it to the 30.4th, and the residual rms falls
     from 0.0411 to 0.0353 uHz. The apparent candidate dissolved AND the estimator
     improved, which is the order those two things should happen in.

The activity model here is therefore quadratic plus hysteresis, and the first
amplitude tested is ZERO.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G, golf_search as SS, golf_bench as B

SEG_D, HW, PER_D, NU_MAX = 180, 4.0, 2300.0, 3.09e-3

def build():
    s1, _ = B.shift_centroid(G.load("PM1"), SEG_D, HW)
    s2, _ = B.shift_centroid(G.load("PM2"), SEG_D, HW)
    sh = np.nanmean(np.vstack([s1, s2]), axis=0)
    act = SS.activity_on_grid(len(sh), SEG_D)
    dact = np.gradient(np.where(np.isfinite(act), act, np.nan))
    return sh, act, dact

def resid_of(x, act, dact):
    m = np.isfinite(x) & np.isfinite(act) & np.isfinite(dact)
    A = np.vstack([np.ones(m.sum()), act[m], act[m]**2, dact[m]]).T
    co, *_ = np.linalg.lstsq(A, x[m], rcond=None)
    r = np.full(len(x), np.nan); r[m] = x[m] - A@co
    return r

if __name__ == "__main__":
    t0 = time.time()
    sh, act, dact = build()
    t = np.arange(len(sh))*SEG_D
    r0 = resid_of(sh, act, dact)
    rok = np.isfinite(r0); rv = r0[rok].copy()
    rng = np.random.default_rng(808)
    NPERM = 6000
    mx = np.empty(NPERM)
    for i in range(NPERM):
        pr = np.full(len(r0), np.nan); pr[rok] = rng.permutation(rv)
        mx[i] = SS.periodogram(pr, rok)[1].max()
    THR = float(np.quantile(mx, 0.95))
    obs = float(SS.periodogram(r0, rok)[1].max())
    print("estimator: PM1+PM2, %d d, +/-%.0f uHz; activity model quadratic + dA/dt" % (SEG_D, HW))
    print("  %d usable segments; residual rms %.4f uHz" % (rok.sum(), np.nanstd(r0)*1e6))
    print("  threshold %.4f; observed residual max %.4f (%.1f pctile)  %s\n"
          % (THR, obs, 100*(mx < obs).mean(), "OK" if obs < THR else "*** STILL ABOVE ***"), flush=True)

    print("  %10s %13s %10s" % ("amp (uHz)", "fractional", "recovery"))
    out = []
    for amp_u in (0.000, 0.010, 0.020, 0.030, 0.045, 0.065, 0.090, 0.130):
        amp = amp_u*1e-6; rec = 0; NREP = 300
        for k in range(NREP):
            inj = sh + amp*np.sin(2*np.pi*t/PER_D + rng.uniform(0, 2*np.pi))
            r = resid_of(inj, act, dact); rk = np.isfinite(r)
            if SS.periodogram(r, rk)[1].max() >= THR: rec += 1
        frac = amp/NU_MAX
        flag = ""
        if amp_u == 0.0:
            flag = "  <-- CONTROL, must be ~5%" + ("" if rec/NREP < 0.15 else "   *** FAILS ***")
        out.append((amp_u, frac, rec/NREP))
        print("  %10.3f %13.2e %9.1f%%%s" % (amp_u, frac, 100*rec/NREP, flag), flush=True)

    mono = all(out[i][2] <= out[i+1][2] + 0.08 for i in range(len(out)-1))
    print("\n  monotonic in amplitude: %s" % ("yes" if mono else "NO -- curve is not trustworthy"))
    hit = [o for o in out if o[2] >= 0.95 and o[0] > 0]
    if out[0][2] < 0.15 and mono and hit:
        a, f, _ = hit[0]
        print("  95%% recovery at %.3f uHz = %.2e fractional at nu_max" % (a, f))
        print("  band estimator, injected: 1.94e-05   ->  improvement %.2fx" % (1.94e-5/f))
        print("  designer level 1e-6  ->  margin %.3f (was 0.05)" % (1e-6/f))
    else:
        print("  NO LIMIT CLAIMED -- control or monotonicity failed")
    json.dump({"threshold": THR, "obs": obs, "monotonic": bool(mono),
               "curve": [{"amp_uHz": a, "frac": f, "rec": r} for a, f, r in out]},
              open(os.path.expanduser("~/golf_inject3.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
