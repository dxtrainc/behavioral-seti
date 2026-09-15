"""
ROW 10, FAST BAND: 5-minute HARPS-N solar, periods 11 min to 2.8 h.

At daily cadence row 10 floored at 9.4e-05 because the daily series sits 25-43x
above the photon noise -- in the 4-40 day band the limit is solar rotation and
active-region evolution. Within an observing session the Sun is far quieter:
measured residual is 1.21x the 7.96e-05 photon floor on contrast and 1.44x on
equivalent width. That is the paper's own claim that the archives are more
informative at minutes than at days, and this is the search it implies.

SESSIONS, not days. Grouping by floor(BJD) gives a median span of 23.7 h per
bucket, which is impossible for a daytime telescope: the BJD boundary falls
mid-observation and each bucket straddles a night. A session is a run with no gap
over an hour -- 1,919 of them, median 5.53 h and 56 exposures.

THE NULL MUST PRESERVE THE WITHIN-SESSION CORRELATION. The residual is not white:
contrast shows lag-1 autocorrelation +0.415 decaying to +0.047 by lag 10, EW
+0.155 and white past lag 2. Permuting residuals would destroy that and
manufacture significance -- the same failure that gave 100% and 95% zero arms in
the daily search. So surrogates are drawn per session from that session's own
empirical autocovariance, at the OBSERVED times. Sessions are ~56 points, so the
factorisation is trivial where the daily search needed 2444x2444.

The surrogate carries no phase relationship BETWEEN sessions, which is the point:
the test asks whether a coherent signal persists across the decade, and the null
is the same noise without that coherence. The window function -- a 1/day comb
putting ~125 alias teeth in this band -- is identical for data and surrogate
because both live on the same timestamps.

CAVEAT, stated because it bounds the result: 5-minute p-modes sit at 288 c/d,
above the 133 c/d Nyquist of a 325 s cadence, and alias into this band near
22 c/d (~65 min). A candidate there is solar before it is anything else.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

SPLIT = 1.0 / 24.0


def build(chan, qc=True):
    d = pd.read_pickle("/home/dxtra/harpsn_sun.pkl")
    num = lambda c: pd.to_numeric(d[c], errors="coerce")
    m = d["spectro_drs_qc"].astype(bool) if qc else pd.Series(True, index=d.index)
    t = num("obj_date_bjd")
    y = (num("spectro_ccf_contrast") if chan == "contrast"
         else num("spectro_ccf_contrast") * num("spectro_ccf_fwhm"))
    ok = m & t.notna() & y.notna() & np.isfinite(y)
    tt, yy = t[ok].values, y[ok].values
    o = np.argsort(tt)
    return tt[o], yy[o]


def sessionise(t, y, minpts=30):
    s = np.concatenate([[0], np.cumsum(np.diff(t) > SPLIT)])
    out = []
    for k in np.unique(s):
        m = s == k
        if m.sum() >= minpts:
            out.append((t[m], y[m]))
    return out


def detrend(ts, ys):
    r = ys / np.median(ys)
    c = np.polyfit(ts - ts.mean(), r, 1)
    return r - np.polyval(c, ts - ts.mean())


class SessionNull:
    """Per-session Cholesky of the session's own empirical autocovariance."""

    def __init__(self, segs, maxlag=25):
        self.L = []
        for ts, x in segs:
            n = len(x)
            xc = x - x.mean()
            v = float(np.dot(xc, xc) / n)
            ac = [1.0]
            for k in range(1, min(maxlag, n - 2)):
                ac.append(float(np.dot(xc[:-k], xc[k:]) / (n * v)) if v > 0 else 0.0)
            ac = np.array(ac)
            idx = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
            C = v * np.where(idx < len(ac), ac[np.clip(idx, 0, len(ac) - 1)], 0.0)
            C[np.diag_indices(n)] = v
            Lk = None
            for j in (0.0, 1e-9, 1e-7, 1e-5):
                try:
                    Lk = np.linalg.cholesky(C + np.eye(n) * (j * v)); break
                except np.linalg.LinAlgError:
                    continue
            if Lk is None:
                w, V = np.linalg.eigh(C)
                Lk = V * np.sqrt(np.clip(w, 1e-12 * v, None))
            self.L.append(Lk)

    def draw(self, rng):
        return np.concatenate([Lk @ rng.standard_normal(Lk.shape[0]) for Lk in self.L])


def power(t, x, f):
    return float(LombScargle(t, x).power(f, method="fast").max())


CFG = {}


def _work(job):
    amp, trial = job
    rng = np.random.default_rng(31337 + trial * 71 + int(amp * 1e9))
    segs, null, f, f_inj = CFG["segs"], CFG["null"], CFG["f"], CFG["f_inj"]
    xs, ts = [], []
    ph0 = rng.uniform(0, 2 * np.pi)
    for tseg, yseg in segs:
        ys = yseg * (1.0 + amp * np.sin(2 * np.pi * f_inj * tseg + ph0)) if amp > 0 else yseg
        xs.append(detrend(tseg, ys)); ts.append(tseg)
    t_all, x_all = np.concatenate(ts), np.concatenate(xs)
    obs = power(t_all, x_all, f)
    ne = 0
    for _ in range(CFG["nsh"]):
        if power(t_all, null.draw(rng), f) >= obs:
            ne += 1
    return (1.0 + ne) / (1.0 + CFG["nsh"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chan", default="ew", choices=["contrast", "ew"])
    ap.add_argument("--period", type=float, default=0.75, help="injected period, hours")
    ap.add_argument("--shifts", type=int, default=120)
    ap.add_argument("--trials", type=int, default=12)
    ap.add_argument("--pool", type=int, default=20)
    ap.add_argument("--out", default="~/fastband.json")
    A = ap.parse_args()

    t, y = build(A.chan)
    segs = sessionise(t, y)
    segs_d = [(ts, detrend(ts, ys)) for ts, ys in segs]
    npts = sum(len(s[0]) for s in segs)
    print("\nROW 10 FAST BAND -- %s, HARPS-N solar 5-minute" % A.chan)
    print("  %d sessions, %d exposures, %.2f yr"
          % (len(segs), npts, (t.max() - t.min()) / 365.25))
    print("  residual scatter %.3e  (%.2fx the 7.96e-05 photon floor)"
          % (np.std(np.concatenate([s[1] for s in segs_d])),
             np.std(np.concatenate([s[1] for s in segs_d])) / 7.96e-5))
    fmin, fmax = 24.0 / 2.8, 24.0 / (11.0 / 60.0)      # 2.8 h down to 11 min
    f = np.linspace(fmin, fmax, 40000)
    print("  band %.1f .. %.1f c/d  (periods 11 min .. 2.8 h), injected at %.2f h"
          % (fmin, fmax, A.period))
    print("  %d surrogates per test, %d trials per amplitude\n" % (A.shifts, A.trials), flush=True)

    global CFG
    CFG = dict(segs=segs, null=SessionNull(segs_d), f=f,
               f_inj=24.0 / A.period, nsh=A.shifts)
    amps = [0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]

    from multiprocessing import Pool
    print("  %10s %9s %9s   %s" % ("amplitude", "p<0.05", "p<0.01", "median p"))
    print("  " + "-" * 50)
    curve = []
    pool = Pool(A.pool)
    try:
        for a in amps:
            ps = np.array(pool.map(_work, [(a, k) for k in range(A.trials)]))
            curve.append(dict(amp=a, rec05=float((ps < 0.05).mean()),
                              rec01=float((ps < 0.01).mean()), med_p=float(np.median(ps))))
            print("  %10.1e %8.1f%% %8.1f%%   %.4f"
                  % (a, 100 * (ps < 0.05).mean(), 100 * (ps < 0.01).mean(),
                     np.median(ps)), flush=True)
    finally:
        pool.close(); pool.join()

    r = [c["rec05"] for c in curve]
    a95 = None
    for i in range(1, len(r)):
        if r[i] >= 0.95 > r[i - 1]:
            fr = (0.95 - r[i - 1]) / max(r[i] - r[i - 1], 1e-9)
            a95 = float(np.exp(np.log(max(amps[i - 1], 1e-12)) +
                               fr * (np.log(amps[i]) - np.log(max(amps[i - 1], 1e-12)))))
            break
    print("\n  zero arm     : %.1f%%  %s" % (100 * r[0],
          "(correct)" if r[0] <= 0.10 else "*** FIRES ON NOTHING ***"))
    print("  95%% recovery : %s" % ("%.2e" % a95 if a95 else "not reached on this grid"))
    print("  daily-cadence result for comparison: 9.4e-05")
    json.dump(dict(chan=A.chan, sessions=len(segs), npts=npts, curve=curve, a95=a95),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("  saved %s" % A.out)


if __name__ == "__main__":
    main()
