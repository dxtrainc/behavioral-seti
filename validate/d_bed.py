"""rank3 against the phase-0 bed, with the ORIGINAL IAAFT null.

rank3 measured A/N = 0.14 against the incumbent's 3.85, so the whole point is
that the surrogate construction should no longer matter much. The plain IAAFT
null -- the best-understood of the three, and the one already in the paper -- is
therefore used unchanged. If rank3 is as construction-insensitive as A/N says,
this passes with the null we already have and no further surrogate engineering
is needed.

The bed is identical to phase 0, same seeds, same triples: real clustering, real
marginals, real spectra, and NO three-way structure by construction. A calibrated
statistic must return median p ~ 0.5, median z ~ 0, and no trend with length.
"""
import argparse, json, os, sys
import numpy as np
from scipy import stats as sps

sys.path.insert(0, "/home/dxtra/beacon-repo/validate")
from phase0 import prep, make_triple                        # noqa: E402
from tripsurr import Plan                                   # noqa: E402


def _cr(x):
    return sps.rankdata(x, method="average") / (len(x) + 1.0) - 0.5


def rank3(d):
    """Third-order rank co-moment of the residuals."""
    return float(np.mean(_cr(d[0]) * _cr(d[1]) * _cr(d[2])))


def rank3d(d):
    """Same, on the DIFFERENCED residuals: keeps sensitivity to fast structure
    while the rank transform keeps the construction artifact bounded."""
    e = [np.diff(x) for x in d]
    return float(np.mean(_cr(e[0]) * _cr(e[1]) * _cr(e[2])))


STATS = {"rank3": rank3, "rank3d": rank3d}


def evaluate(tri, nsh, iters, rng, sf=None):
    sf = sf or rank3
    D = [prep(x)[1] for x in tri]
    obs = sf(D)
    if not np.isfinite(obs):
        return None
    plan = Plan(np.vstack(tri))
    null = []
    for _ in range(nsh):
        Q = [prep(x)[1] for x in plan.draw(rng, iters)]
        v = sf(Q)
        if np.isfinite(v):
            null.append(v)
    if len(null) < nsh // 4:
        return None
    a = np.asarray(null)
    # two-sided: three-way structure can raise OR lower the co-moment
    ne = int((np.abs(a - a.mean()) >= abs(obs - a.mean())).sum())
    return dict(p=(1.0 + ne) / (1.0 + a.size), ne=ne, obs=obs,
                z=float((obs - a.mean()) / max(a.std(), 1e-12)), n=len(tri[0]))


CFG = {}


def _work(job):
    t, chans, n, nm = job
    r = np.random.default_rng(1000 + t)
    tri = make_triple(chans, n, CFG["ltrend"], CFG["lres"], r)
    out = evaluate(tri, CFG["shifts"], CFG["iters"], r,
                   STATS[CFG.get("stat", "rank3")])
    if out:
        out["combo"] = " / ".join(nm)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntrip", type=int, default=60)
    ap.add_argument("--shifts", type=int, default=4000)
    ap.add_argument("--iters", type=int, default=6)
    ap.add_argument("--ltrend", type=int, default=1460)
    ap.add_argument("--lres", type=int, default=180)
    ap.add_argument("--pool", type=int, default=60)
    ap.add_argument("--stat", default="rank3", choices=list(STATS))
    ap.add_argument("--out", default="~/candidate_d.json")
    A = ap.parse_args()

    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M, names = z["M"], [str(x) for x in z["names"]]
    fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    usable = [i for i in range(M.shape[0]) if fin[i].sum() >= 3000]

    print("\nCANDIDATE D -- %s on the phase-0 bed, IAAFT null unchanged" % A.stat)
    print("  %d triples x %d draws\n" % (A.ntrip, A.shifts), flush=True)

    rng = np.random.default_rng(20260915)      # identical bed to phase 0
    jobs = []
    for t in range(A.ntrip):
        pick = rng.choice(usable, size=3, replace=False)
        n = int(rng.choice([3000, 4500, 6000, 9000, 13000]))
        chans = [M[i][fin[i]][:max(n, 3000)] for i in pick]
        n = min(n, min(len(c) for c in chans))
        jobs.append((t, [c[:n] for c in chans], n, [names[i] for i in pick]))

    global CFG
    CFG = dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts,
               iters=A.iters, stat=A.stat)
    from multiprocessing import Pool
    with Pool(A.pool) as p:
        R = [x for x in p.map(_work, jobs) if x]

    ps = np.array([r["p"] for r in R]); zs = np.array([r["z"] for r in R])
    ns = np.array([r["n"] for r in R])
    ks = sps.kstest(ps, "uniform")
    print("  %d triples" % len(R))
    print("  median p %.4f   median z %+.3f   min p %.4f" % (np.median(ps), np.median(zs), ps.min()))
    print("  KS vs uniform: D=%.3f p=%.2e" % (ks.statistic, ks.pvalue))
    print("\n                   current stat   rank3")
    print("    median p         0.9738        %.4f" % np.median(ps))
    print("    median z        -1.87          %+.3f" % np.median(zs))
    print("    KS D             0.563         %.3f" % ks.statistic)
    print("\n  by length (flat = no length artifact):")
    for lo, hi in ((0, 4000), (4000, 6000), (6000, 99999)):
        m = (ns >= lo) & (ns < hi)
        if m.sum():
            print("    n=%5d..%-5d %3d  median p %.4f  median z %+.3f"
                  % (lo, hi, m.sum(), np.median(ps[m]), np.median(zs[m])))
    ok = (0.35 < np.median(ps) < 0.65) and abs(np.median(zs)) < 0.4 and ks.pvalue > 0.01
    print("\n  VERDICT: %s" % ("PASSES the phase-0 gate" if ok else "FAILS"))
    json.dump(dict(median_p=float(np.median(ps)), median_z=float(np.median(zs)),
                   ks_D=float(ks.statistic), ks_p=float(ks.pvalue),
                   passed=bool(ok), results=R),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("  saved %s" % A.out)


if __name__ == "__main__":
    main()
