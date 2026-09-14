"""
A PAIRWISE-PRESERVING NULL FOR THE THREE-BODY STATISTIC, AND ITS ZERO ARM.

fptest.py voided both triple forms: log-2nd-diff fired 9x too often in the far
tail of synthetic triples containing no signal, and resid-triple could not fire
at all. Section 4.9 diagnosed the cause as the null rather than the statistic --
the circular shift destroys the PAIRWISE correlations along with the three-way
alignment, so the aligned value is compared against surrogates that have lost
structure it was never supposed to be tested against.

WHAT THE REPLACEMENT HAS TO DO. Destroy three-way structure; preserve everything
else. "Everything else" here means each channel's marginal distribution, each
channel's power spectrum, and every pairwise cross-spectrum.

THE CONSTRUCTION. Apply ONE common random phase sequence phi(f) to all three
channels at once. Relative phases between channels are untouched, so every
cross-spectrum -- hence every pairwise cross-correlation at every lag -- is
preserved exactly. The bispectrum is not: its phase picks up
phi(f1) + phi(f2) - phi(f1+f2), which is random. Second order survives, third
order does not. Rank-mapping each surrogate back onto the observed marginal then
restores the static nonlinearity, and iterating the two steps restores the power
spectrum that rank-mapping perturbs. This is multivariate IAAFT, and its null
hypothesis is the one we actually want to test against: "three channels that are
static nonlinear transforms of a jointly Gaussian linear process with the
observed cross-spectra."

A NOTE ON WHY THIS STATISTIC IS TESTABLE AT ALL. Phase randomisation is a no-op
for any statistic that depends only on the power spectrum -- such a control
cannot fail and is worthless. stat() is kurtosis-based, a higher-order moment,
so the randomisation bites. That is checked directly in --diag rather than
assumed.

BOTH NULLS RUN ON THE SAME SYNTHETIC TRIPLES IN THE SAME PROCESS, so the
comparison against the published circular-shift numbers (6.6% / 0.9% / -0.04 for
log-2nd-diff; 0.0% / 0.0% / -4.49 for resid-triple) is exact and does not depend
on reproducing an earlier run. The circular-shift arm is the positive control:
if it does not reproduce those numbers, the harness is wrong, not the surrogate.

The channel generator, prep/clean/stat and the two forms are VERBATIM from
fptest.py so the test bed is identical.
"""
import sys, time, argparse, json, os
import numpy as np
from multiprocessing import Pool
from scipy import stats
from scipy.ndimage import median_filter

W = 181

# ---------------------------------------------------------------- verbatim from fptest.py
def runmed(y, w): return median_filter(y, size=w, mode="nearest")

def clean(x):
    d = np.diff(x); s = 1.4826 * np.median(np.abs(d - np.median(d)))
    if s <= 0: return None
    d = np.clip(d, -3.0 * s, 3.0 * s); y = np.concatenate([[x[0]], x[0] + np.cumsum(d)])
    sd = y.std(); return (y - y.mean()) / sd if sd > 0 else None

def stat(x):
    z = clean(x)
    if z is None: return np.nan
    d = np.diff(z); s = d.std()
    if s <= 0: return np.nan
    return float(stats.kurtosis(d, fisher=True) - 4.0 * np.median(np.abs(d)) / s)

def prep(x):
    r = x / np.median(x); d = r - runmed(r, W); sd = d.std()
    return r, (d / sd if sd > 0 else d)

def rollp(p, s): return (np.roll(p[0], s), np.roll(p[1], s))

def l2d(pa, pb, pc):
    la = np.log(np.maximum(pa[0], 1e-12)); lb = np.log(np.maximum(pb[0], 1e-12))
    lc = np.log(np.maximum(pc[0], 1e-12))
    return la - 2.0 * lb + lc

def rtp(pa, pb, pc): return pa[1] * pb[1] * pc[1]

def red(rng, n, alpha=0.92):
    e = rng.standard_normal(n); x = np.zeros(n)
    for i in range(1, n): x[i] = alpha * x[i - 1] + e[i]
    return np.exp(0.35 * (x - x.mean()) / max(x.std(), 1e-9)) + 0.05

def channels(rng, n):
    """Three channels sharing a driver: pairwise correlation, NO three-way signal."""
    drv = red(rng, n, 0.97)
    w = rng.uniform(0.4, 0.9, size=3)
    return [drv ** w[k] * red(rng, n, 0.90) ** (1 - w[k]) for k in range(3)]

FORMS = (("log-2nd-diff", l2d), ("resid-triple", rtp))

# ---------------------------------------------------------------- the surrogate
class Plan(object):
    """Per-triple constants for the surrogate. Built once, drawn from many times.

    Recomputing the sort order and the FFT on every draw costs more than the
    draw itself; the observed series never changes within a triple.
    """
    __slots__ = ("n", "k", "srt", "Arel")

    def __init__(self, X):
        self.k, self.n = X.shape
        self.srt = np.sort(X, axis=1)                 # target marginals
        F0 = np.fft.rfft(X, axis=1)
        # amplitude spectra carrying the OBSERVED relative phases: applying one
        # common phase to this preserves every cross-spectrum exactly.
        self.Arel = np.abs(F0) * np.exp(1j * (np.angle(F0) - np.angle(F0[0])[None, :]))

    def draw(self, rng, iters=20):
        n, k, srt, Arel = self.n, self.k, self.srt, self.Arel
        ph = rng.uniform(0.0, 2.0 * np.pi, size=Arel.shape[1])
        ph[0] = 0.0
        if n % 2 == 0:
            ph[-1] = 0.0
        Y = np.fft.irfft(Arel * np.exp(1j * ph)[None, :], n=n, axis=1)
        for _ in range(iters):
            for j in range(k):                        # rank-map onto the marginal
                Y[j, np.argsort(Y[j])] = srt[j]
            # only the common phase evolves; relative phases stay locked to data
            phi0 = np.angle(np.fft.rfft(Y[0]))
            Y = np.fft.irfft(Arel * np.exp(1j * phi0)[None, :], n=n, axis=1)
        for j in range(k):                            # exact marginals on exit
            Y[j, np.argsort(Y[j])] = srt[j]
        return Y


def mviaaft(X, rng, iters=20):
    return Plan(X).draw(rng, iters)


# ---------------------------------------------------------------- one triple
def one(job):
    seed, nsh, iters, arms = job
    rng = np.random.default_rng(seed)
    n = NDAY
    ch = channels(rng, n)
    P = [prep(c) for c in ch]
    X = np.vstack(ch)
    plan = Plan(X) if "surr" in arms else None
    out = {}

    obs = {}
    for nm, fn in FORMS:
        v = stat(fn(*P))
        if np.isfinite(v):
            obs[nm] = v
    if not obs:
        return None

    for arm in arms:
        null = {nm: [] for nm in obs}
        for _ in range(nsh):
            if arm == "shift":
                s1, s2 = int(rng.integers(1, n)), int(rng.integers(1, n))
                Q = (P[0], rollp(P[1], s1), rollp(P[2], s2))
            else:
                Q = tuple(prep(r) for r in plan.draw(rng, iters))
            for nm, fn in FORMS:
                if nm not in obs:
                    continue
                v = stat(fn(*Q))
                if np.isfinite(v):
                    null[nm].append(v)
        for nm in obs:
            a = np.asarray(null[nm])
            if a.size == 0:
                continue
            ne = int((a >= obs[nm]).sum())
            out["%s|%s" % (nm, arm)] = dict(
                p=(1.0 + ne) / (1.0 + a.size), ne=ne, n=int(a.size),
                z=float((obs[nm] - a.mean()) / max(a.std(), 1e-12)))
    return out

# ---------------------------------------------------------------- surrogate fidelity
def diag(seed, iters, ndraw=40):
    """Does the surrogate preserve what it claims and destroy what it claims?"""
    rng = np.random.default_rng(seed)
    n = NDAY
    ch = channels(rng, n); X = np.vstack(ch)

    def xcorr(M, i, j, lags=(0, 1, 2, 5, 10, 20, 50)):
        a = (M[i] - M[i].mean()) / max(M[i].std(), 1e-12)
        b = (M[j] - M[j].mean()) / max(M[j].std(), 1e-12)
        return np.array([float(np.mean(a * np.roll(b, L))) for L in lags])

    A0 = np.abs(np.fft.rfft(X, axis=1))
    ref_x = {(i, j): xcorr(X, i, j) for i in range(3) for j in range(i + 1, 3)}
    ref_b = np.array([stat(f(*[prep(c) for c in ch])) for _, f in FORMS])

    plan = Plan(X)
    t0 = time.time()
    se, xe, bs, mar = [], [], [], []
    for _ in range(ndraw):
        Y = plan.draw(rng, iters)
        A1 = np.abs(np.fft.rfft(Y, axis=1))
        se.append(float(np.median(np.abs(A1 - A0) / np.maximum(A0, 1e-12))))
        xe.append(max(float(np.max(np.abs(xcorr(Y, i, j) - ref_x[(i, j)])))
                      for i in range(3) for j in range(i + 1, 3)))
        mar.append(float(np.max(np.abs(np.sort(Y, axis=1) - np.sort(X, axis=1)))))
        bs.append([stat(f(*[prep(r) for r in Y])) for _, f in FORMS])
    dt = (time.time() - t0) / ndraw
    bs = np.asarray(bs)

    print("  surrogate fidelity  (%d draws, %d iterations, n=%d)" % (ndraw, iters, n))
    print("    marginal, max abs deviation of sorted values : %.3e   (want 0)" % max(mar))
    print("    power spectrum, median relative error        : %.4f    (want << 1)" % np.median(se))
    print("    pairwise cross-corr, WORST abs error at lags : %.4f    (want << 1)" % max(xe))
    print("    cost per surrogate                           : %.1f ms" % (dt * 1e3))
    print()
    print("  does the randomisation actually bite? (statistic on data vs surrogates)")
    for k, (nm, _) in enumerate(FORMS):
        col = bs[:, k]; col = col[np.isfinite(col)]
        if col.size < 2:
            print("    %-13s : no finite surrogate values" % nm); continue
        sd = col.std()
        print("    %-13s data %+8.3f   surrogate %+8.3f +/- %.3f   -> %s"
              % (nm, ref_b[k], col.mean(), sd,
                 "MOVES, control can fail" if sd > 1e-6 else "*** NO-OP ***"))
    return 0

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nday", type=int, default=5000)
    ap.add_argument("--ntrip", type=int, default=320)
    ap.add_argument("--seed0", type=int, default=1000,
                    help="first seed; disjoint tranches pool by using disjoint ranges")
    ap.add_argument("--shifts", type=int, default=10000)
    ap.add_argument("--iters", type=int, default=20)
    ap.add_argument("--pool", type=int, default=72)
    ap.add_argument("--arms", default="shift,surr")
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--out", default="tripsurr.json")
    A = ap.parse_args()
    NDAY = A.nday

    if A.diag:
        sys.exit(diag(1000, A.iters))

    arms = tuple(A.arms.split(","))
    ALPHA = 0.05 / 10996
    print("pairwise-preserving null, zero arm")
    print("  %d synthetic triples x %d days, %d draws per null, arms=%s, iters=%d, pool=%d"
          % (A.ntrip, A.nday, A.shifts, "+".join(arms), A.iters, A.pool))
    print("  synthetic channels contain NO three-way signal; Bonferroni alpha = %.2e\n"
          % ALPHA, flush=True)

    t0 = time.time()
    jobs = [(s, A.shifts, A.iters, arms) for s in range(A.seed0, A.seed0 + A.ntrip)]
    with Pool(A.pool) as p:
        R = [r for r in p.map(one, jobs) if r]
    el = (time.time() - t0) / 60.0
    print("done in %.1f min\n" % el, flush=True)

    floor = 1.0 / (A.shifts + 1)
    print("  %-13s %-8s %7s %8s %9s %10s %9s"
          % ("form", "null", "p<0.05", "p<1e-3", "at floor", "median z", "n"))
    print("  %-13s %-8s %7s %8s %9s %10s %9s"
          % ("expected", "sound", "5.0%", "0.1%", "0.0%", "0.00", ""))
    print("  " + "-" * 70)
    res = {}
    for nm, _ in FORMS:
        for arm in arms:
            key = "%s|%s" % (nm, arm)
            v = [r[key] for r in R if key in r]
            if not v:
                continue
            ps = np.array([x["p"] for x in v]); zs = np.array([x["z"] for x in v])
            res[key] = dict(n=len(v), p=[float(x) for x in ps],
                            z=[float(x) for x in zs], p05=float((ps < 0.05).mean()),
                            p1e3=float((ps < 1e-3).mean()),
                            at_floor=float((ps <= floor + 1e-15).mean()),
                            med_z=float(np.median(zs)),
                            below_alpha=int((ps < ALPHA).sum()))
            print("  %-13s %-8s %6.1f%% %7.1f%% %8.1f%% %10.2f %9d"
                  % (nm, arm, 100 * res[key]["p05"], 100 * res[key]["p1e3"],
                     100 * res[key]["at_floor"], res[key]["med_z"], len(v)))
    print()
    for nm, _ in FORMS:
        a, b = "%s|shift" % nm, "%s|surr" % nm
        if a in res and b in res:
            print("  %-13s far-tail excess over nominal:  shift %5.1fx   surrogate %5.1fx"
                  % (nm, res[a]["p1e3"] / 1e-3, res[b]["p1e3"] / 1e-3))
    json.dump(dict(nday=A.nday, ntrip=A.ntrip, seed0=A.seed0, shifts=A.shifts, iters=A.iters,
                   alpha=ALPHA, floor=floor, elapsed_min=el, results=res),
              open(os.path.join(os.path.expanduser("~"), A.out), "w"), indent=1)
    print("\n  saved ~/%s" % A.out)
