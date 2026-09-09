"""
Fast-band viewpoint test: turning the null into a limit.

The search found no Earth-line peak at all in the judged band (125 s to 6 h),
so there is no candidate to test for viewpoint-specificity. What the search then
bounds is the amplitude at which a modulation present in one line and absent in
the other WOULD have been seen. Measured by injection, at 95% recovery.

WHERE THE INJECTION GOES, which the first version of this file got wrong. It
injected into the series AFTER the daily-profile fold, so the injected tone was
never subject to the fold while a real signal would have been. That overstates
sensitivity, and at a daily harmonic it overstates it completely. The injection
now enters BEFORE the fold and passes through the identical pipeline.

AND WHICH PERIODS. The fold removes every harmonic of 1/86400 s. The first run
used 300 s and 3600 s, which are the 288th and 24th harmonics -- so the honest
version of that run reports no sensitivity at all there, and it is a real blind
spot of the pipeline rather than a fact about the Sun. Quotable limits are
therefore measured at periods NOT commensurate with the day.
"""
import os, json
import numpy as np
import viewpoint_fast as V

D = os.path.expanduser("~")
rng = np.random.default_rng(41)

def main():
    grid, E, M = V.load_grid()
    okE = np.isfinite(E) & (E > 0); okM = np.isfinite(M) & (M > 0)
    baseE = V.prep(np.nan_to_num(E, nan=1.0), okE)
    baseM = V.prep(np.nan_to_num(M, nan=1.0), okM)
    fM0, PM0, _ = V.spec(V.fold_out(baseM, okM, 1440)); RM0 = PM0/V.cont(PM0)
    b = (fM0 > 1/(8*3600.)) & (fM0 < 1/(3*3600.))
    Porb = 1.0/float(fM0[b][np.argmax(RM0[b])])/60.0
    n = len(baseE); t = np.arange(n, dtype=float)*V.STEP

    _, SEGP = V.fold_drifting(baseM, okM, seg_days=30)   # measured once

    def pipeline(r, ok, mars):
        """identical to the search: fold, then spectrum"""
        if mars:
            r = V.fold_drifting(r, ok, seg_days=30, periods_in=SEGP)
            r = V.fold_out(r, ok, 1440)
        else:
            r = V.fold_out(r, ok, 1440)
        f, P, W = V.spec(r)
        return f, P

    out = {}
    for tag, base, ok, mars in (("Earth (GOES EUVS)", baseE, okE, False),
                                ("Mars (MAVEN EUVM)", baseM, okM, True)):
        f, P = pipeline(base, ok, mars)
        C = V.cont(P); thr, _ = V.thr_of(P/C)
        sd = base[ok].std()
        print("\n%s   threshold %.1f   sigma %.3e" % (tag, thr, sd), flush=True)
        rows = {}
        for Pd, label in ((307.0, "307 s"), (911.0, "911 s"), (3671.0, "3671 s"),
                          (3600.0, "3600 s (daily harmonic)")):
            f0 = 1.0/Pd
            k = int(np.argmin(np.abs(f-f0))); sl = slice(max(0, k-2), k+3)
            lim = None; curve = []
            for a in (0.0, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3):
                hits = 0; N = 30
                for _ in range(N):
                    sur = np.roll(base, int(rng.integers(1000, n-1000)))
                    y = sur + a*sd*np.cos(2*np.pi*f0*t + rng.uniform(0, 2*np.pi))
                    y[~ok] = 0.0
                    _, P2 = pipeline(y, ok, mars)
                    if (P2[sl]/C[sl]).max() > thr: hits += 1
                pw = 100.0*hits/N
                curve.append((a, pw))
                if lim is None and a > 0 and pw >= 95: lim = a
            frac = (lim*sd) if lim else None
            print("  P=%-24s %s  95%% at %s" %
                  (label, " ".join("%.0e:%3.0f%%" % c for c in curve),
                   ("%.1e sigma = %.1e fractional" % (lim, frac)) if lim else ">5e-3 sigma"),
                  flush=True)
            rows[label] = {"curve": curve, "lim95_sigma": lim, "lim95_frac": frac}
        out[tag] = {"sigma": float(sd), "rows": rows}
    json.dump(out, open(os.path.join(D, "viewpoint_fast_limit.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_fast_limit.json", flush=True)

if __name__ == "__main__":
    main()
