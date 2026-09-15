"""
CANDIDATE A: a surrogate that preserves the clustering IAAFT destroys.

DIAGNOSIS. IAAFT keeps each channel's values and power spectrum but scrambles
global Fourier phase, and phase is what keeps solar extremes bunched into
episodes. De-clustered extremes jump around more than clustered ones, so the
surrogates out-score the data on a kurtosis statistic almost every time:
median z -2.41 on the Sun, and worsening with record length.

CONSTRUCTION. Decompose each channel with an a trous (undecimated, shift-
invariant) wavelet transform, which reconstructs by simple summation:

    x = smooth + sum_j detail_j

Then randomise the SIGNS of the detail coefficients. The arithmetic is the whole
argument:

  |coefficient| is untouched at every time and scale  -> the energy stays exactly
      where it was, so episodic CLUSTERING is preserved rather than smeared;
  the SAME sign pattern is applied to all three channels, so a pairwise product
      picks up s^2 = +1 -> every cross-relationship is preserved EXACTLY;
  a three-way product picks up s^3 = s              -> three-way alignment is
      randomised, which is the thing under test.

Second order preserved, third order destroyed, clustering intact. That is the
separation the Fourier surrogate could not achieve.

SIGNS FLIP IN BLOCKS, NOT POINTWISE. Detail coefficients at scale j are
correlated over ~2^j samples; flipping each independently would whiten them and
destroy the power spectrum the null must keep. Block length therefore scales with
the level, so the local waveform survives and only its polarity changes.

The SMOOTH component is never flipped: it carries the solar cycle, which is the
source of the pairwise correlation the null is required to preserve.

Tested against the phase-0 bed, where the current null gives median p 0.974 and
median z -1.87 on triples containing NO three-way structure. Acceptance is
median p ~ 0.5, median z ~ 0, and no trend with record length.
"""
import argparse, json, os, sys
import numpy as np
from scipy.ndimage import median_filter
from scipy import stats as sps

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase0 import prep, stat, make_triple, factor, block_idx   # noqa: E402

H = np.array([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0                  # B3 spline


def _smooth(c, step):
    n = len(c)
    out = np.zeros(n)
    base = np.arange(n)
    for k, w in zip((-2, -1, 0, 1, 2), H):
        out += w * c[np.clip(base + k * step, 0, n - 1)]
    return out


def atrous(x, J):
    """Undecimated wavelet transform. x == smooth + sum(details), exactly."""
    c = x.astype(float).copy()
    det = []
    for j in range(J):
        cj = _smooth(c, 1 << j)
        det.append(c - cj)
        c = cj
    return np.array(det), c


class WaveletSignPlan:
    """Drop-in for tripsurr.Plan: .draw(rng, iters) -> (k, n) surrogate."""

    def __init__(self, X, J=None, blk=4.0, rankmap=True, shared_scale=False):
        X = np.asarray(X, float)
        self.k, self.n = X.shape
        self.J = J or max(3, min(9, int(np.log2(self.n)) - 4))
        self.blk, self.rankmap = blk, rankmap
        self.shared_scale = shared_scale
        self.srt = np.sort(X, axis=1)
        d, s = [], []
        for row in X:
            dj, sj = atrous(row, self.J)
            d.append(dj); s.append(sj)
        self.det = np.array(d)          # (k, J, n)
        self.smo = np.array(s)          # (k, n)
        # block length per level: signs must be constant over ~the level's own
        # correlation time or the spectrum is whitened
        self.L = [max(2, int(self.blk * (1 << j))) for j in range(self.J)]

    def draw(self, rng, iters=None):
        n = self.n
        if self.shared_scale:
            # ONE sign per time block, applied to EVERY scale at once, so
            # cross-scale alignment -- what makes a transient sharp -- survives.
            L = max(2, int(self.blk * (1 << (self.J - 1))))
            nb = int(np.ceil(n / L))
            row = np.repeat(rng.choice((-1.0, 1.0), size=nb), L)[:n]
            S = np.broadcast_to(row, (self.J, n))
        else:
            S = np.empty((self.J, n))
            for j in range(self.J):
                L = self.L[j]
                nb = int(np.ceil(n / L))
                S[j] = np.repeat(rng.choice((-1.0, 1.0), size=nb), L)[:n]
        Y = (self.det * S[None, :, :]).sum(axis=1) + self.smo
        if self.rankmap:
            for a in range(self.k):
                Y[a, np.argsort(Y[a])] = self.srt[a]
        return Y


def evaluate(tri, nsh, rng, blk, rankmap, shared=False):
    P = [prep(x) for x in tri]
    obs = stat(P[0][1] * P[1][1] * P[2][1])
    if not np.isfinite(obs):
        return None
    plan = WaveletSignPlan(np.vstack(tri), blk=blk, rankmap=rankmap,
                           shared_scale=shared)
    null = []
    for _ in range(nsh):
        Q = [prep(r) for r in plan.draw(rng)]
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
    t, chans, n, nm = job
    r = np.random.default_rng(1000 + t)
    tri = make_triple(chans, n, CFG["ltrend"], CFG["lres"], r)   # SAME bed
    out = evaluate(tri, CFG["shifts"], r, CFG["blk"], CFG["rankmap"],
                   CFG.get("shared", False))
    if out:
        out["combo"] = " / ".join(nm)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntrip", type=int, default=60)
    ap.add_argument("--shifts", type=int, default=4000)
    ap.add_argument("--ltrend", type=int, default=1460)
    ap.add_argument("--lres", type=int, default=180)
    ap.add_argument("--blk", type=float, default=4.0)
    ap.add_argument("--norankmap", action="store_true")
    ap.add_argument("--shared", action="store_true")
    ap.add_argument("--pool", type=int, default=60)
    ap.add_argument("--out", default="~/candidate_a.json")
    A = ap.parse_args()

    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M, names = z["M"], [str(x) for x in z["names"]]
    fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    usable = [i for i in range(M.shape[0]) if fin[i].sum() >= 3000]

    print("\nCANDIDATE A -- wavelet block-sign surrogate, on the phase-0 bed")
    print("  sign block factor %.1f x 2^j, rank-map %s" % (A.blk, not A.norankmap))
    print("  %d triples x %d draws\n" % (A.ntrip, A.shifts), flush=True)

    rng = np.random.default_rng(20260915)          # SAME seed as phase 0
    jobs = []
    for t in range(A.ntrip):
        pick = rng.choice(usable, size=3, replace=False)
        n = int(rng.choice([3000, 4500, 6000, 9000, 13000]))
        chans = [M[i][fin[i]][:max(n, 3000)] for i in pick]
        n = min(n, min(len(c) for c in chans))
        jobs.append((t, [c[:n] for c in chans], n, [names[i] for i in pick]))

    global CFG
    CFG = dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts,
               blk=A.blk, rankmap=not A.norankmap, shared=A.shared)
    from multiprocessing import Pool
    with Pool(A.pool) as p:
        R = [x for x in p.map(_work, jobs) if x]

    ps = np.array([r["p"] for r in R]); zs = np.array([r["z"] for r in R])
    ns = np.array([r["n"] for r in R])
    ks = sps.kstest(ps, "uniform")
    print("  %d triples" % len(R))
    print("  median p %.4f   median z %+.2f   min p %.4f" % (np.median(ps), np.median(zs), ps.min()))
    print("  KS vs uniform: D=%.3f p=%.2e" % (ks.statistic, ks.pvalue))
    print("\n                        current null    candidate A")
    print("    median p             0.9738          %.4f" % np.median(ps))
    print("    median z            -1.87            %+.2f" % np.median(zs))
    print("\n  by length (flat = clustering no longer destroyed):")
    for lo, hi in ((0, 4000), (4000, 6000), (6000, 99999)):
        m = (ns >= lo) & (ns < hi)
        if m.sum():
            print("    n=%5d..%-5d %3d  median p %.4f  median z %+.2f"
                  % (lo, hi, m.sum(), np.median(ps[m]), np.median(zs[m])))
    passed = (0.3 < np.median(ps) < 0.7) and abs(np.median(zs)) < 0.6
    print("\n  VERDICT: %s" % ("PASSES the phase-0 gate" if passed else
                               "FAILS -- still miscalibrated"))
    json.dump(dict(blk=A.blk, rankmap=not A.norankmap, median_p=float(np.median(ps)),
                   median_z=float(np.median(zs)), ks_D=float(ks.statistic),
                   passed=bool(passed), results=R),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("  saved %s" % A.out)


if __name__ == "__main__":
    main()
