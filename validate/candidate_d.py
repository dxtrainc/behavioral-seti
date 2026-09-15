"""
CANDIDATE D: find a three-way statistic that is not dominated by how the
surrogate was built.

THE FINDING THAT SENDS US HERE. On one triple, three different surrogate
constructions gave null distributions of essentially the same WIDTH (sd 0.0225 -
0.0240) but centred 0.05-0.06 apart -- roughly twice that width. The current
statistic, kurtosis of the differenced triple product, therefore reports the
choice of construction more loudly than it reports the sky. That is not a
tuning problem and no further surrogate engineering fixes it.

THE FIGURE OF MERIT is artifact-to-noise:

    A/N  =  (spread of null MEANS across constructions) / (typical null sd)

The current statistic sits near 2.4. Anything below ~0.3 is usable: the
construction choice would then move the answer by less than a third of the
null's own width. This script measures A/N for several candidates on the same
triples, which is far cheaper and more decisive than running each through the
full phase-0 bed.

CANDIDATES. All are computed on the DETRENDED residuals, since the trend carries
the pairwise structure the null is supposed to preserve.

  cur    kurtosis(diff(clean(d_a d_b d_c))) - 4 MAD/sd     -- the incumbent
  rank3  mean( R(d_a) R(d_b) R(d_c) )                      -- third-order rank
         co-moment: R is the centred normalised rank, so the statistic is
         invariant to any monotone transform of each channel, is a FIRST moment
         of bounded quantities, and still cannot factor into pairwise terms.
  rank3d same, on the differenced residuals -- keeps sensitivity to fast
         structure while staying rank-stabilised.
  sgn3   mean( sign(d_a) sign(d_b) sign(d_c) )             -- the extreme case:
         depends only on orthant, so it discards magnitude entirely.
"""
import os, sys
import numpy as np
from scipy import stats as sps

sys.path.insert(0, "/home/dxtra/beacon-repo/validate")
from phase0 import prep, stat as cur_stat, make_triple           # noqa: E402
from candidate_a import WaveletSignPlan                          # noqa: E402
from tripsurr import Plan                                        # noqa: E402


def _cr(x):
    """Centred, normalised rank in (-0.5, 0.5)."""
    r = sps.rankdata(x, method="average")
    return r / (len(x) + 1.0) - 0.5


def s_rank3(d):
    return float(np.mean(_cr(d[0]) * _cr(d[1]) * _cr(d[2])))


def s_rank3d(d):
    e = [np.diff(x) for x in d]
    return float(np.mean(_cr(e[0]) * _cr(e[1]) * _cr(e[2])))


def s_sgn3(d):
    return float(np.mean(np.sign(d[0]) * np.sign(d[1]) * np.sign(d[2])))


def s_cur(d):
    return cur_stat(d[0] * d[1] * d[2])


STATS = [("cur", s_cur), ("rank3", s_rank3), ("rank3d", s_rank3d), ("sgn3", s_sgn3)]


def main():
    ndraw = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    ntri = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M = z["M"]; fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    usable = [i for i in range(M.shape[0]) if fin[i].sum() >= 4600]
    rng = np.random.default_rng(7)

    acc = {k: {"means": [], "sds": []} for k, _ in STATS}
    for t in range(ntri):
        pick = rng.choice(usable, size=3, replace=False)
        n = 4500
        chans = [M[i][fin[i]][:n] for i in pick]
        tri = make_triple(chans, n, 1460, 180, np.random.default_rng(100 + t))
        D = [prep(x)[1] for x in tri]
        obs = {k: f(D) for k, f in STATS}

        cons = (("iaaft", Plan(np.vstack(tri))),
                ("sign-scale", WaveletSignPlan(np.vstack(tri), blk=4.0)),
                ("sign-shared", WaveletSignPlan(np.vstack(tri), blk=4.0, shared_scale=True)))
        per = {k: [] for k, _ in STATS}
        for cname, plan in cons:
            r = np.random.default_rng(3 + t)
            vals = {k: [] for k, _ in STATS}
            for _ in range(ndraw):
                Q = [prep(x)[1] for x in plan.draw(r, 6)]
                for k, f in STATS:
                    v = f(Q)
                    if np.isfinite(v):
                        vals[k].append(v)
            for k, _ in STATS:
                a = np.asarray(vals[k])
                per[k].append((a.mean(), a.std()))
        for k, _ in STATS:
            ms = np.array([m for m, _ in per[k]]); ss = np.array([s for _, s in per[k]])
            acc[k]["means"].append(ms.max() - ms.min())
            acc[k]["sds"].append(np.median(ss))
        print("  triple %d done" % (t + 1), flush=True)

    print("\n  %-8s %12s %12s %10s   %s" % ("stat", "artifact", "null sd", "A/N", "verdict"))
    print("  " + "-" * 62)
    for k, _ in STATS:
        art = float(np.median(acc[k]["means"])); sd = float(np.median(acc[k]["sds"]))
        an = art / sd if sd > 0 else float("inf")
        print("  %-8s %12.3e %12.3e %10.2f   %s"
              % (k, art, sd, an,
                 "usable" if an < 0.3 else ("marginal" if an < 1 else "dominated by construction")))


if __name__ == "__main__":
    main()
