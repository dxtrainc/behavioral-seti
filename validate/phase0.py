"""
PHASE 0: A ZERO ARM THAT REPRODUCES THE FAILURE.

The pairwise-preserving null passed its synthetic zero arm at median z -0.17 and
then failed on the real Sun at -2.41. The control was the problem, not just the
null: AR(1) x lognormal channels have no episodic clustering, which is precisely
what the surrogate destroys and therefore precisely what the test bed had to
contain. A control that cannot reproduce a failure cannot validate a fix.

CONSTRUCTION. Each real channel factors multiplicatively into a slow component
and a residual:

    channel = trend x residual,   trend = 181-day running median

  * the TREND is block-bootstrapped ONCE and shared by all three synthetic
    channels, which is what creates their pairwise correlation -- the common
    solar cycle, exactly as in the real data;
  * each RESIDUAL is block-bootstrapped INDEPENDENTLY from that channel's own
    real residual, so each synthetic channel inherits real episodic clustering,
    real marginals and real short-range spectrum;
  * blocks are contiguous and circular, so clustering INSIDE a block survives
    intact. Block length is the knob that controls how much clustering is kept.

Three-way structure beyond what the shared trend implies is therefore absent BY
CONSTRUCTION, while the clustering that breaks the null is present and real.
That is the definition of a zero arm for this test.

ACCEPTANCE (this file only builds and measures the bed; it does not fix anything):
the CURRENT null must FAIL on it the way it fails on the Sun --
    median p ~ 0.95-0.99, median z ~ -2, and median p RISING with series length.
If it does not, the bed is not reproducing the failure and no fix can be
validated against it.
"""
import argparse, json, os, sys
import numpy as np
from scipy.ndimage import median_filter
from scipy import stats as sps

sys.path.insert(0, os.path.expanduser("~/beacon-diag/validate"))
from tripsurr import Plan                                    # noqa: E402

W = 181


def prep(x):
    r = x / np.median(x)
    d = r - median_filter(r, size=W, mode="nearest")
    sd = d.std()
    return r, (d / sd if sd > 0 else d)


def stat(x):
    d = np.diff(x)
    s = 1.4826 * np.median(np.abs(d - np.median(d)))
    if s <= 0:
        return np.nan
    d = np.clip(d, -3.0 * s, 3.0 * s)
    y = np.concatenate([[x[0]], x[0] + np.cumsum(d)])
    sd = y.std()
    if sd <= 0:
        return np.nan
    z = (y - y.mean()) / sd
    d = np.diff(z)
    s = d.std()
    if s <= 0:
        return np.nan
    return float(sps.kurtosis(d, fisher=True) - 4.0 * np.median(np.abs(d)) / s)


def block_idx(n_out, n_src, L, rng):
    """Circular block bootstrap indices: contiguous runs of length L."""
    out = np.empty(n_out, dtype=np.int64)
    pos = 0
    while pos < n_out:
        start = int(rng.integers(0, n_src))
        take = min(L, n_out - pos)
        out[pos:pos + take] = (start + np.arange(take)) % n_src
        pos += take
    return out


def factor(x):
    """channel -> (trend, residual) with channel = trend * residual."""
    t = median_filter(x, size=W, mode="nearest")
    t = np.where(t <= 0, np.median(x[x > 0]) if np.any(x > 0) else 1.0, t)
    return t, x / t


def make_triple(chans, n, L_trend, L_res, rng):
    """Shared bootstrapped trend x independently bootstrapped real residuals."""
    tr0, _ = factor(chans[0])
    ish = block_idx(n, len(tr0), L_trend, rng)          # ONE schedule, shared
    trend = tr0[ish]
    out = []
    for c in chans:
        _, res = factor(c)
        idx = block_idx(n, len(res), L_res, rng)        # independent per channel
        out.append(trend * res[idx])
    return out


def evaluate(tri, nsh, iters, rng):
    """Run the CURRENT (raw-path) null on one synthetic triple."""
    P = [prep(x) for x in tri]
    obs = stat(P[0][1] * P[1][1] * P[2][1])
    if not np.isfinite(obs):
        return None
    plan = Plan(np.vstack(tri))
    null = []
    for _ in range(nsh):
        Q = [prep(r) for r in plan.draw(rng, iters)]
        v = stat(Q[0][1] * Q[1][1] * Q[2][1])
        if np.isfinite(v):
            null.append(v)
    if len(null) < nsh // 4:
        return None
    a = np.asarray(null)
    ne = int((a >= obs).sum())
    return dict(p=(1.0 + ne) / (1.0 + a.size), ne=ne,
                z=float((obs - a.mean()) / max(a.std(), 1e-12)), n=len(tri[0]))


CFG = {}


def _work(job):
    """Module-level so multiprocessing can pickle it."""
    t, chans, n, nm = job
    r = np.random.default_rng(1000 + t)
    tri = make_triple(chans, n, CFG["ltrend"], CFG["lres"], r)
    out = evaluate(tri, CFG["shifts"], CFG["iters"], r)
    if out:
        out["combo"] = " / ".join(nm)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntrip", type=int, default=60)
    ap.add_argument("--shifts", type=int, default=4000)
    ap.add_argument("--iters", type=int, default=6)
    ap.add_argument("--ltrend", type=int, default=1460, help="trend block, days")
    ap.add_argument("--lres", type=int, default=180, help="residual block, days")
    ap.add_argument("--pool", type=int, default=10)
    ap.add_argument("--out", default="~/phase0.json")
    A = ap.parse_args()

    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M, names = z["M"], [str(x) for x in z["names"]]
    fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    usable = [i for i in range(M.shape[0]) if fin[i].sum() >= 3000]

    print("\nPHASE 0 -- zero arm from block-bootstrapped REAL channels")
    print("  %d usable channels, trend block %d d, residual block %d d"
          % (len(usable), A.ltrend, A.lres))
    print("  %d synthetic triples x %d draws, iters=%d\n" % (A.ntrip, A.shifts, A.iters),
          flush=True)

    rng = np.random.default_rng(20260915)
    jobs = []
    for t in range(A.ntrip):
        pick = rng.choice(usable, size=3, replace=False)
        # length drawn across the real overlap range so the n-trend is testable
        n = int(rng.choice([3000, 4500, 6000, 9000, 13000]))
        chans = [M[i][fin[i]][:max(n, 3000)] for i in pick]
        if min(len(c) for c in chans) < n:
            n = min(len(c) for c in chans)
        jobs.append((t, [c[:n] for c in chans], n, [names[i] for i in pick]))

    from multiprocessing import Pool
    global CFG
    CFG = dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts, iters=A.iters)
    with Pool(A.pool) as p:
        R = [x for x in p.map(_work, jobs) if x]

    ps = np.array([r["p"] for r in R]); zs = np.array([r["z"] for r in R])
    ns = np.array([r["n"] for r in R])
    print("  %d triples completed\n" % len(R))
    print("  median p %.4f   median z %+.2f   min p %.4f" % (np.median(ps), np.median(zs), ps.min()))
    ks = sps.kstest(ps, "uniform")
    print("  KS vs uniform: D=%.3f p=%.2e" % (ks.statistic, ks.pvalue))
    print("\n  median p by length (the real run showed 0.952 / 0.984 / 0.998):")
    for lo, hi in ((0, 4000), (4000, 6000), (6000, 99999)):
        m = (ns >= lo) & (ns < hi)
        if m.sum():
            print("    n=%5d..%-5d  %3d  median p %.4f  median z %+.2f"
                  % (lo, hi, m.sum(), np.median(ps[m]), np.median(zs[m])))
    ok = (np.median(ps) > 0.90) and (np.median(zs) < -1.0)
    print("\n  REPRODUCES THE FAILURE: %s" % ("YES -- usable as a test bed" if ok
          else "NO -- this bed cannot validate a fix"))
    json.dump(dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts, iters=A.iters,
                   median_p=float(np.median(ps)), median_z=float(np.median(zs)),
                   reproduces=bool(ok), results=R),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("  saved %s" % A.out)


if __name__ == "__main__":
    main()
