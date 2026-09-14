"""
THE PAIR SPACE AGAINST A PAIRWISE-PRESERVING NULL.

The published pair sweep uses a circular shift, which is the correct null for
the question it asks: "is there any alignment-specific coupling between these
two channels?" It destroys the cross-alignment and preserves everything else,
and for a two-body statistic there is nothing else to preserve.

THIS ASKS A DIFFERENT AND SHARPER QUESTION. The multivariate IAAFT surrogate
preserves each channel's marginal, each channel's power spectrum AND the pairwise
cross-spectrum exactly -- so it preserves the entire LINEAR relationship between
the two channels -- while randomising the common phase, which destroys the
bispectrum. Its null hypothesis is therefore:

    "these two channels are static nonlinear transforms of a jointly Gaussian
     linear process with the observed cross-spectrum"

and rejecting it means there is coupling BEYOND the linear cross-correlation.
A beacon that modulates one observable by another produces exactly that, and a
linear search walks past it. Same 1,074 tests, same archives, new question.

WHY THE STATISTIC ADMITS THE TEST. stat() is kurtosis(diff) - 4*MAD(diff)/sd, a
higher-order moment. Common-phase randomisation is a NO-OP for any statistic
that sees only the power spectrum, and such a control cannot fail; this one is
not, and --diag measures that rather than assuming it.

EMPIRICAL COUNTING ONLY. No GPD extrapolation: p = (1+ne)/(1+M) with ne
reported, so the floor is explicit and every p-value is reproducible. The
surrogate budget must put that floor BELOW the family-wise threshold or the run
cannot fire -- the trap of section 5.2 -- so --shifts is checked against it and
the run refuses to start if it would be a no-op.

sweepT.py has no __main__ guard, so it is exec'd up to the job-building marker
rather than imported: that runs its channel loading and leaves its sweep alone.
"""
import argparse, json, os, sys, time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "validate"))
from tripsurr import Plan                      # noqa: E402  the surrogate


def load_sweepT(path):
    """Run sweepT.py's channel loading and nothing else."""
    src = open(path, encoding="utf-8").read()
    marker = "# ---- build the job list"
    if marker not in src:
        raise SystemExit("job-list marker not found in %s" % path)
    saved = sys.argv
    sys.argv = ["sweepT.py", "--mode", "pairs", "--shifts", "1", "--pool", "1"]
    ns = {"__name__": "sweepT_channels", "__file__": path}
    try:
        exec(compile(src[:src.index(marker)], path, "exec"), ns)
    finally:
        sys.argv = saved
    return ns


S = load_sweepT(os.path.join(HERE, "sweepT.py"))
M, names, N = S["M"], S["names"], S["N"]
prep, forms_pre, stat = S["prep"], S["forms_pre"], S["stat"]

# Each form consumes EITHER the normalised ratio r OR the detrended residual d,
# never both, so the surrogate is built on whichever the form actually reads.
# Surrogating the PREPARED array rather than the raw one takes prep()'s 181-wide
# median filter out of the inner loop -- the same economy the circular-shift null
# gets by rolling prepared arrays -- and states the null where the statistic
# actually lives: on the detrended residuals.
USES_R = {"ratio", "log-ratio"}


def form_one(fname, U, V):
    """One form, matching forms_pre() exactly."""
    if fname == "ratio":
        return U / V
    if fname == "log-ratio":
        return np.log(np.maximum(U, 1e-12)) - np.log(np.maximum(V, 1e-12))
    if fname == "resid-product":
        return U * V
    if fname == "resid-ratio":
        return U / np.where(np.abs(V) < 0.1, np.nan, V)
    raise KeyError(fname)


def one(arg):
    i, j, fname, nsh, iters, seed = arg
    ov = np.isfinite(M[i]) & np.isfinite(M[j])
    n = int(ov.sum())
    if n < MINOV:
        return None
    pa, pb = prep(M[i][ov]), prep(M[j][ov])
    f = forms_pre(pa[0], pa[1], pb[0], pb[1])
    if fname not in f:
        return None
    obs = stat(f[fname])
    if not np.isfinite(obs):
        return None

    U, V = (pa[0], pb[0]) if fname in USES_R else (pa[1], pb[1])
    # the fast path must reproduce the canonical one exactly
    if not np.allclose(form_one(fname, U, V), f[fname], equal_nan=True):
        return None

    plan = Plan(np.vstack([U, V]))
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(nsh):
        Y = plan.draw(rng, iters)
        g = form_one(fname, Y[0], Y[1])
        if not np.all(np.isfinite(g)):
            continue
        v = stat(g)
        if np.isfinite(v):
            null.append(v)
    if len(null) < nsh // 4:
        return None
    arr = np.asarray(null)
    ne = int((arr >= obs).sum())
    return dict(combo="%s / %s" % (names[i], names[j]), i=i, j=j, form=fname,
                obs=float(obs), p=(1.0 + ne) / (1.0 + arr.size), ne=ne,
                nnull=int(arr.size), ndays=n,
                z=float((obs - arr.mean()) / max(arr.std(), 1e-12)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shifts", type=int, default=100000)
    ap.add_argument("--iters", type=int, default=12)
    ap.add_argument("--pool", type=int, default=84)
    ap.add_argument("--minov", type=int, default=2000)
    ap.add_argument("--ntests", type=int, default=0, help="0 = all; >0 truncates, for timing")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--out", default="sweepP_surr.json")
    A = ap.parse_args()
    MINOV = A.minov

    import itertools
    jobs = []
    for i, j in itertools.combinations(range(N), 2):
        ov = np.isfinite(M[i]) & np.isfinite(M[j])
        if ov.sum() < MINOV:
            continue
        pa, pb = prep(M[i][ov]), prep(M[j][ov])
        for fn in forms_pre(pa[0], pa[1], pb[0], pb[1]):
            jobs.append((i, j, fn))
    ntest = len(jobs)
    thr = A.alpha / max(ntest, 1)
    floor = 1.0 / (A.shifts + 1)

    print("\npair space against the pairwise-preserving null")
    print("  %d tests, Bonferroni threshold %.3e" % (ntest, thr))
    print("  %d surrogate draws -> empirical floor %.3e" % (A.shifts, floor))
    if floor >= thr:
        raise SystemExit(
            "  REFUSING TO RUN: the floor is at or above the threshold, so no test\n"
            "  could fire whatever the data contain -- the section 5.2 trap.\n"
            "  Need --shifts > %d." % int(np.ceil(1.0 / thr)))
    print("  floor is %.1fx below the threshold -- the run can fire\n"
          % (thr / floor), flush=True)

    if A.diag:
        # does the randomisation bite on REAL channels, form by form?
        print("  does the surrogate move the statistic on real data?")
        seen = set()
        for i, j, fn in jobs:
            if fn in seen:
                continue
            seen.add(fn)
            ov = np.isfinite(M[i]) & np.isfinite(M[j])
            pa, pb = prep(M[i][ov]), prep(M[j][ov])
            o = stat(forms_pre(pa[0], pa[1], pb[0], pb[1])[fn])
            U, V = (pa[0], pb[0]) if fn in USES_R else (pa[1], pb[1])
            plan = Plan(np.vstack([U, V])); rng = np.random.default_rng(7)
            v, t0 = [], time.time()
            for _ in range(60):
                Y = plan.draw(rng, A.iters)
                g = form_one(fn, Y[0], Y[1])
                if np.all(np.isfinite(g)):
                    sv = stat(g)
                    if np.isfinite(sv):
                        v.append(sv)
            ms = (time.time() - t0) / 60 * 1e3
            v = np.asarray(v)
            print("    %-13s n=%5d  data %+8.3f  surrogate %+8.3f +/- %.3f  %6.1f ms/draw  -> %s"
                  % (fn, int(ov.sum()), o, v.mean(), v.std(), ms,
                     "MOVES" if v.std() > 1e-6 else "*** NO-OP ***"), flush=True)
        sys.exit(0)

    if A.ntests:
        jobs = jobs[:A.ntests]
    jobs = [(i, j, fn, A.shifts, A.iters, abs(hash((i, j, fn))) % (2 ** 32))
            for i, j, fn in jobs]
    print("  running %d of %d tests, pool=%d\n" % (len(jobs), ntest, A.pool), flush=True)

    t0 = time.time()
    with Pool(A.pool) as p:
        R = [r for r in p.map(one, jobs, chunksize=1) if r]
    el = time.time() - t0
    print("completed %d tests in %.1f min  (%.0f draws/s)"
          % (len(R), el / 60, len(R) * A.shifts / max(el, 1)), flush=True)

    json.dump(dict(shifts=A.shifts, iters=A.iters, ntests_total=ntest,
                   threshold=thr, floor=floor, elapsed_s=el, results=R),
              open(os.path.join(os.path.expanduser("~"), A.out), "w"))
    print("saved ~/%s" % A.out)

    R.sort(key=lambda r: r["p"])
    surv = [r for r in R if r["p"] < thr]
    at_floor = [r for r in R if r["ne"] == 0]
    print("\n  survivors below %.2e : %d" % (thr, len(surv)))
    print("  tests pinned at the floor (ne=0) : %d" % len(at_floor))
    print("\n  %-44s %-13s %9s %11s %5s" % ("combination", "form", "stat", "p", "ne"))
    for r in R[:20]:
        print("  %-44s %-13s %9.4f %11.3e %5d"
              % (r["combo"][:44], r["form"], r["obs"], r["p"], r["ne"]))
