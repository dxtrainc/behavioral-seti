"""
SECTION 4.12 FAST BAND -- REPLACE THE ANALYTIC THRESHOLD WITH A MEASURED ONE.

The band-wide false-alarm control added on review returned 197.28 exceedances per
surrogate over 1,336,602 judged bins where 0.05 was expected -- the analytic threshold
mu*ln(N/alpha) is about 4,000x too permissive on this data. That factor matches, from a
different dataset and a different method, section 4.2's independently measured finding
that the exponential tail model is optimistic by ~4e3 on GOES.

So the 4.2e-3 limit is not quoted at family-wise 0.05. It is quoted against a threshold
that fires two hundred times per realisation of pure noise.

THE FIX IS THE ONE ALREADY USED FOR THE QUADRUPLE SWEEP: take the threshold from the
distribution of the MAXIMUM of R over the judged band across surrogates, and read the
95th percentile off it. That is Westfall-Young, it is exact under whatever correlation
and tail shape the data actually has, and it assumes nothing about exponentiality.

A SECOND FINDING FROM THE SAME CONTROL. The cached continuum gives 375 false alarms
against 197 for a continuum recomputed on each surrogate -- the shortcut nearly doubles
the rate. It was justified in the source by appeal to the depth-0 control, which tested
one bin against a threshold set for half a million and could not have detected it. The
continuum is recomputed here.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import viewpoint_fast as V

NSUR = int(os.environ.get("NSUR", "300"))

if __name__ == "__main__":
    t0 = time.time()
    grid, E, M = V.load_grid()
    okE = np.isfinite(E) & (E > 0); okM = np.isfinite(M) & (M > 0)
    rE = V.prep(np.nan_to_num(E, nan=1.0), okE)
    rM = V.prep(np.nan_to_num(M, nan=1.0), okM)
    rE = V.fold_out(rE, okE, 1440)
    rM, _segP = V.fold_drifting(rM, okM)
    rM = V.fold_out(rM, okM, 1440)
    fE, PE, W = V.spec(rE); fM, PM, _ = V.spec(rM)
    inband = (1.0/fE <= 6*3600.0) & (1.0/fE >= 125.0)
    print("judged band: %d bins" % inband.sum(), flush=True)

    rng = np.random.default_rng(5150)
    out = {}
    for tag, r, ok, P0 in (("Earth", rE, okE, PE), ("Mars", rM, okM, PM)):
        n = len(r)
        analytic, _ = V.thr_of(P0/V.cont(P0))
        mx = np.empty(NSUR)
        for i in range(NSUR):
            sur = np.roll(r, int(rng.integers(1000, n-1000))).copy(); sur[~ok] = 0.0
            _, P2, _ = V.spec(sur)
            mx[i] = float((P2/V.cont(P2))[inband].max())
        emp = float(np.quantile(mx, 0.95))
        print("\n  %s" % tag, flush=True)
        print("    analytic threshold  mu*ln(N/alpha)     %8.1f" % analytic)
        print("    EMPIRICAL threshold 95th pct of max-R  %8.1f" % emp)
        print("    ratio %.1fx -- the analytic value is that much too permissive" % (emp/analytic))
        # amplitude scales as sqrt(threshold) for a coherent tone against a fixed continuum
        print("    limits scale as sqrt(threshold): %.2fx looser" % np.sqrt(emp/analytic))
        out[tag] = dict(analytic=float(analytic), empirical=emp,
                        ratio=float(emp/analytic), amp_factor=float(np.sqrt(emp/analytic)))
    print("\n  SECTION 4.12 FAST-BAND LIMITS, REQUOTED")
    for tag, old in (("Earth", 7.2e-6), ("Mars", 4.2e-3)):
        f = out[tag]["amp_factor"]
        print("    %-6s %.1e  ->  %.1e   (x%.2f)" % (tag, old, old*f, f))
    print("\n  The binding number is Mars: a viewpoint claim needs both lines to see it.")
    json.dump(out, open(os.path.expanduser("~/vp_empirical.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
