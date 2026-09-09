"""
ITEM A, part 3 -- turning the viewpoint null into a limit.

Part 2 found no Earth-only peak that Mars had the power to see. An absence is
worth a number only if the number is measured, so: inject a coherent tone into a
surrogate Mars series across the accessible band and find the fractional
amplitude at which recovery reaches 95%, the convention used throughout.
"""
import os, json
import numpy as np
from scipy.ndimage import median_filter
import viewpoint_l2b as V
from viewpoint_search import spec, cont, thr_of

D = os.path.expanduser("~")
rng = np.random.default_rng(23)

def main():
    jm, vm = V.maven_daily(); jg, vg = V.goes_daily()
    dm, ym, _ = V.to_daily(jm, vm); dg, yg, _ = V.to_daily(jg, vg)
    lo = max(dm.min(), dg.min()); hi = min(dm.max(), dg.max())
    grid = np.arange(lo, hi+1)
    M = np.interp(grid, dm, ym)
    le, re, lm, rm = V.longitudes(grid)
    dlam = (lm-le) % 360.0
    tau = dlam/360.0*V.CARR
    tau = np.where(tau > V.CARR/2, tau-V.CARR, tau)
    Mal = np.interp(grid+tau, grid, M)
    Mal = Mal - median_filter(Mal, size=401, mode="nearest")
    f, P = spec(Mal); R = P/cont(P)
    keep = (1.0/f) <= 200.0
    tM, _ = thr_of(R[keep])
    sd = Mal.std(); n = len(grid); t = np.arange(n, dtype=float)
    print("Mars threshold %.2f   series sigma %.3e" % (tM, sd), flush=True)
    print("\n  %8s %s" % ("period", "recovery vs fractional amplitude"), flush=True)
    out = {}
    for Pd in (3.0, 7.0, 13.5, 27.0, 60.0, 120.0, 180.0):
        row = []
        lim = None
        for a in (0.0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35):
            hits = 0; N = 200
            for _ in range(N):
                sur = np.roll(Mal, int(rng.integers(50, n-50)))
                y = sur + a*sd*np.cos(2*np.pi*t/Pd + rng.uniform(0, 2*np.pi))
                f2, P2 = spec(y); R2 = P2/cont(P2)
                k = int(np.argmin(np.abs(f2-1.0/Pd)))
                if R2[max(0, k-2):k+3].max() > tM: hits += 1
            pw = 100.0*hits/N
            row.append((a, pw))
            if lim is None and a > 0 and pw >= 95: lim = a
        print("  %7.1fd  %s   95%% at %s" %
              (Pd, " ".join("%.2f:%3.0f%%" % r for r in row),
               ("%.2f" % lim) if lim else ">0.35"), flush=True)
        out["%.1f" % Pd] = {"curve": row, "lim95": lim}
    json.dump(out, open(os.path.join(D, "viewpoint_limit.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_limit.json", flush=True)

main()
