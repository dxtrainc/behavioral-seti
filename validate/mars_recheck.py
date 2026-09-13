"""
RE-MEASURE THE MAVEN MARS THRESHOLD WITH A CORRECTED SURROGATE.

Commit 265a064 changed section 4.12's fast-band limit from 4.2e-3 to 1.6e-2 on the
strength of an empirical threshold of 377.3 against an analytic 24.8. That measurement
used np.roll followed by re-masking, which on a gapped series inflates the surrogate's
zero fraction -- on EUVS at 73.1% coverage it took 26.9% zeros to 40.4%, and the measured
threshold from 484.5 to 279,228, a factor of 576 entirely of my own making.

So the committed requote is very likely overstated and this re-derives it the correct
way: roll only the OBSERVED samples and leave the gap pattern exactly where it is.

Whatever comes back, the direction is not in doubt -- the analytic threshold is too
permissive on every channel tested. What is in doubt is by how much, and a published
number was changed before that was established.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.expanduser("~/beacon-repo/search"))
import viewpoint_fast as V

NSUR = int(os.environ.get("NSUR", "120"))

if __name__ == "__main__":
    t0 = time.time()
    grid, E, M = V.load_grid()
    okE = np.isfinite(E) & (E > 0); okM = np.isfinite(M) & (M > 0)
    rE = V.fold_out(V.prep(np.nan_to_num(E, nan=1.0), okE), okE, 1440)
    rM = V.prep(np.nan_to_num(M, nan=1.0), okM)
    rM, _ = V.fold_drifting(rM, okM)
    rM = V.fold_out(rM, okM, 1440)
    fE, PE, W = V.spec(rE)
    inband = (1.0/fE <= 6*3600.0) & (1.0/fE >= 125.0)
    print("judged band: %d bins\n" % inband.sum(), flush=True)

    out = {}
    for tag, r, ok in (("Earth", rE, okE), ("Mars", rM, okM)):
        f, P, _ = V.spec(r)
        analytic, _ = V.thr_of(P/V.cont(P))
        idx = np.where(ok)[0]; vals = r[idx]
        rng = np.random.default_rng(8675309)
        mx, zf = [], []
        for _ in range(NSUR):
            s = np.zeros_like(r)
            s[idx] = np.roll(vals, int(rng.integers(100, len(vals)-100)))
            _, P2, _ = V.spec(s)
            mx.append(float((P2/V.cont(P2))[inband].max()))
            zf.append(float((s == 0).mean()))
        emp = float(np.quantile(mx, 0.95))
        fac = float(np.sqrt(max(emp/analytic, 1.0)))
        print("  %-6s zeros: data %.1f%%  surrogate %.1f%%" %
              (tag, 100*(r == 0).mean(), 100*np.mean(zf)))
        print("         analytic %.1f   empirical %.1f   ratio %.1fx   limits %.2fx looser"
              % (analytic, emp, emp/analytic, fac), flush=True)
        out[tag] = dict(analytic=float(analytic), empirical=emp, factor=fac)

    print("\n  SECTION 4.12 FAST BAND, RE-DERIVED")
    for tag, printed in (("Earth", 7.2e-6), ("Mars", 4.2e-3)):
        f = out[tag]["factor"]
        print("    %-6s as printed %.1e -> %.1e   (x%.2f)" % (tag, printed, printed*f, f))
    print("\n    commit 265a064 claimed Mars 4.2e-3 -> 1.6e-2 (x3.90, broken surrogate)")
    print("    corrected                     -> %.1e (x%.2f)"
          % (4.2e-3*out["Mars"]["factor"], out["Mars"]["factor"]))
    json.dump(out, open(os.path.expanduser("~/mars_recheck.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
