"""
Fast-band viewpoint test: turning the null into a limit.

The search found no Earth-line peak at all in the judged band (125 s to 6 h),
so there is no candidate to test for viewpoint-specificity. What the search then
bounds is the amplitude at which a modulation present in one line and absent in
the other WOULD have been seen. Measured by injection, at 95% recovery, into
both series separately.
"""
import os, json
import numpy as np
import viewpoint_fast as V

D = os.path.expanduser("~")
rng = np.random.default_rng(41)

def main():
    grid, E, M = V.load_grid()
    okE = np.isfinite(E) & (E > 0); okM = np.isfinite(M) & (M > 0)
    rE = V.prep(np.nan_to_num(E, nan=1.0), okE)
    rM = V.prep(np.nan_to_num(M, nan=1.0), okM)
    fM0, PM0, _ = V.spec(rM); RM0 = PM0/V.cont(PM0)
    b = (fM0 > 1/(8*3600.)) & (fM0 < 1/(3*3600.))
    Porb = 1.0/float(fM0[b][np.argmax(RM0[b])])/60.0
    rE = V.fold_out(rE, okE, 1440)
    rM = V.fold_drifting(rM, okM); rM = V.fold_out(rM, okM, 1440)
    n = len(rE); t = np.arange(n, dtype=float)*V.STEP

    out = {}
    for tag, r, ok in (("Earth (GOES EUVS)", rE, okE), ("Mars (MAVEN EUVM)", rM, okM)):
        f, P, W = V.spec(r); C = V.cont(P); R = P/C
        thr, _ = V.thr_of(R)
        sd = r[ok].std()
        print("\n%s   threshold %.1f   sigma %.3e" % (tag, thr, sd), flush=True)
        rows = {}
        for Pd in (300.0, 900.0, 3600.0, 10800.0):
            f0 = 1.0/Pd
            k = int(np.argmin(np.abs(f-f0))); sl = slice(max(0, k-2), k+3)
            lim = None; curve = []
            for a in (0.0, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3):
                hits = 0; N = 60
                for _ in range(N):
                    sur = np.roll(r, int(rng.integers(1000, n-1000)))
                    y = sur + a*sd*np.cos(2*np.pi*f0*t + rng.uniform(0, 2*np.pi))
                    y[~ok] = 0.0
                    _, P2, _ = V.spec(y)
                    if (P2[sl]/C[sl]).max() > thr: hits += 1
                pw = 100.0*hits/N
                curve.append((a, pw))
                if lim is None and a > 0 and pw >= 95: lim = a
            print("  P=%6.0f s  %s  95%% at %s" %
                  (Pd, " ".join("%.0e:%3.0f%%" % c for c in curve),
                   ("%.1e sigma" % lim) if lim else ">1e-3"), flush=True)
            rows["%.0f" % Pd] = {"curve": curve, "lim95": lim}
        out[tag] = rows
    json.dump(out, open(os.path.join(D, "viewpoint_fast_limit.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_fast_limit.json", flush=True)

main()
