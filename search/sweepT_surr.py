"""
THE TRIPLE SPACE UNDER THE PAIRWISE-PRESERVING NULL.

Section 4.9 reports the triple space as UNEXAMINED, not null: the circular shift
destroys the pairwise correlations along with the three-way alignment, so
log-2nd-diff separates from its null by construction (9x excess in the far tail
of signal-free synthetic triples) and resid-triple sits 4.5 sd below its own
null and cannot fire at all. Neither result is about the sky.

This runs the same enumeration against the multivariate IAAFT surrogate, which
preserves every marginal, every power spectrum and every pairwise cross-spectrum
while randomising the common phase -- so second-order structure survives and the
bispectrum does not. Rejecting it means genuine three-way structure.

BUILT FOR A RUN MEASURED IN DAYS, WHICH CHANGES THE ENGINEERING:

  * every completed test is appended to a JSONL checkpoint and flushed, so a
    crash costs the test in flight and nothing else. sweepT.py writes one
    json.dump at the very end; that cost 155 minutes of finished compute this
    morning when the output path was wrong, and the same bug at this scale
    would cost nine days.
  * --resume reads the checkpoint and skips what is already done, so an
    interrupted run continues rather than restarting.
  * --shard SI/SN splits the space; a dead shard loses its slice only.
  * progress and ETA are reported continuously, because a job this long needs
    to be observable while it runs, not only afterwards.

EMPIRICAL COUNTING ONLY. p = (1+ne)/(1+M), ne reported so the floor is explicit,
no GPD extrapolation anywhere. The run REFUSES to start if the draw budget puts
the floor at or above the family-wise threshold -- the section 5.2 trap, where a
run cannot fire whatever the data contain.
"""
import argparse, itertools, json, os, sys, time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "validate"))
from tripsurr import Plan                                      # noqa: E402


def load_sweepT(path):
    """sweepT.py has no __main__ guard; run its channel loading and stop."""
    src = open(path, encoding="utf-8").read()
    marker = "# ---- build the job list"
    saved = sys.argv
    sys.argv = ["sweepT.py", "--mode", "triples", "--shifts", "1", "--pool", "1"]
    ns = {"__name__": "sweepT_channels", "__file__": path}
    try:
        exec(compile(src[:src.index(marker)], path, "exec"), ns)
    finally:
        sys.argv = saved
    return ns


S = load_sweepT(os.path.join(HERE, "sweepT.py"))
M, names, N = S["M"], S["names"], S["N"]
prep, forms_tri, stat = S["prep"], S["forms_tri"], S["stat"]

# log-2nd-diff reads the normalised ratio r; resid-triple reads the detrended
# residual d. Neither reads both, so the surrogate is built on whichever the
# form consumes -- which also keeps prep()'s 181-wide median filter out of the
# inner loop, exactly as the circular-shift null avoids it by rolling.
USES_R = {"log-2nd-diff"}


def form_one(fname, U, V, W):
    """One three-body form, matching forms_tri() exactly. V is the middle term."""
    if fname == "resid-triple":
        return U * V * W
    if fname == "log-2nd-diff":
        return (np.log(np.maximum(U, 1e-12)) - 2.0 * np.log(np.maximum(V, 1e-12))
                + np.log(np.maximum(W, 1e-12)))
    raise KeyError(fname)


def one(arg):
    i, j, k, fname, nsh, iters, seed = arg
    key = "%d-%d-%d|%s" % (i, j, k, fname)
    try:
        ov = np.isfinite(M[i]) & np.isfinite(M[j]) & np.isfinite(M[k])
        n = int(ov.sum())
        if n < MINOV:
            return None
        P = [prep(M[c][ov]) for c in (i, j, k)]
        f = forms_tri(*P)
        if fname not in f:
            return None
        obs = stat(f[fname])
        if not np.isfinite(obs):
            return None

        idx = 0 if fname in USES_R else 1
        U, V, W = (P[0][idx], P[1][idx], P[2][idx])
        if not np.allclose(form_one(fname, U, V, W), f[fname], equal_nan=True):
            return None                     # fast path must match the canonical one

        plan = Plan(np.vstack([U, V, W]))
        rng = np.random.default_rng(seed)
        null = []
        for _ in range(nsh):
            Y = plan.draw(rng, iters)
            g = form_one(fname, Y[0], Y[1], Y[2])
            if not np.all(np.isfinite(g)):
                continue
            v = stat(g)
            if np.isfinite(v):
                null.append(v)
        if len(null) < nsh // 4:
            return None
        a = np.asarray(null)
        ne = int((a >= obs).sum())
        return dict(key=key, combo="%s / %s / %s" % (names[i], names[j], names[k]),
                    form=fname, obs=float(obs), p=(1.0 + ne) / (1.0 + a.size),
                    ne=ne, nnull=int(a.size), ndays=n,
                    z=float((obs - a.mean()) / max(a.std(), 1e-12)))
    except Exception as e:                  # one bad test must not kill nine days
        return dict(key=key, error="%s: %s" % (type(e).__name__, e))


def hms(s):
    s = int(s)
    return "%dd %02dh %02dm" % (s // 86400, (s % 86400) // 3600, (s % 3600) // 60)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shifts", type=int, default=250000)
    ap.add_argument("--iters", type=int, default=12)
    ap.add_argument("--pool", type=int, default=80)
    ap.add_argument("--minov", type=int, default=2000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--ntests", type=int, default=0, help="0 = all; >0 truncates, for testing")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--ckpt", default="")
    A = ap.parse_args()
    MINOV = A.minov
    SI, SN = (int(x) for x in A.shard.split("/"))

    jobs = []
    for i, j, k in itertools.combinations(range(N), 3):
        ov = np.isfinite(M[i]) & np.isfinite(M[j]) & np.isfinite(M[k])
        if ov.sum() < MINOV:
            continue
        jobs.append((i, j, k, "resid-triple"))
        # log-2nd-diff is not symmetric: each channel takes a turn as the middle
        jobs.append((i, j, k, "log-2nd-diff"))
        jobs.append((j, i, k, "log-2nd-diff"))
        jobs.append((i, k, j, "log-2nd-diff"))
    ntest = len(jobs)
    thr = A.alpha / max(ntest, 1)
    floor = 1.0 / (A.shifts + 1)

    print("\ntriple space against the pairwise-preserving null")
    print("  %d tests (%d triples x 4 forms)" % (ntest, ntest // 4))
    print("  Bonferroni threshold %.3e" % thr)
    print("  %d draws -> empirical floor %.3e" % (A.shifts, floor))
    if floor >= thr:
        raise SystemExit(
            "  REFUSING TO RUN: floor at or above threshold -- no test could fire\n"
            "  whatever the data contain (the section 5.2 trap).\n"
            "  Need --shifts > %d." % int(np.ceil(1.0 / thr)))
    print("  floor is %.1fx below the threshold -- the run can fire" % (thr / floor))

    jobs = [x for n_, x in enumerate(jobs) if n_ % SN == SI]
    ckpt = A.ckpt or os.path.join(os.path.expanduser("~"),
                                  "sweepT_surr_%d-%d.jsonl" % (SI, SN))
    done = {}
    if A.resume and os.path.exists(ckpt):
        for ln in open(ckpt, errors="replace"):
            try:
                r = json.loads(ln)
                done[r["key"]] = r
            except Exception:
                continue
        print("  resuming: %d tests already in %s" % (len(done), ckpt))

    todo = [(i, j, k, fn, A.shifts, A.iters, abs(hash((i, j, k, fn))) % (2 ** 32))
            for (i, j, k, fn) in jobs if "%d-%d-%d|%s" % (i, j, k, fn) not in done]
    if A.ntests:
        todo = todo[:A.ntests]
    print("  shard %d/%d: %d tests, %d to run, pool=%d" % (SI, SN, len(jobs), len(todo), A.pool))
    print("  checkpoint: %s\n" % ckpt, flush=True)

    t0 = time.time()
    nd = 0
    fh = open(ckpt, "a", buffering=1)
    with Pool(A.pool) as p:
        for r in p.imap_unordered(one, todo, chunksize=1):
            if r is None:
                nd += 1
                continue
            fh.write(json.dumps(r) + "\n")
            nd += 1
            if nd % 25 == 0:
                fh.flush(); os.fsync(fh.fileno())
                el = time.time() - t0
                rate = nd / max(el, 1e-9)
                print("  %6d/%-6d  %5.1f%%  elapsed %s  eta %s  (%.2f tests/min)"
                      % (nd, len(todo), 100.0 * nd / max(len(todo), 1), hms(el),
                         hms((len(todo) - nd) / max(rate, 1e-12)), rate * 60), flush=True)
    fh.flush(); os.fsync(fh.fileno()); fh.close()
    el = time.time() - t0
    print("\ncompleted %d tests in %s" % (nd, hms(el)), flush=True)

    R = list(done.values())
    for ln in open(ckpt, errors="replace"):
        try:
            r = json.loads(ln)
            if r["key"] not in done:
                R.append(r)
        except Exception:
            continue
    ok = [r for r in R if "error" not in r]
    bad = [r for r in R if "error" in r]
    out = os.path.join(os.path.expanduser("~"), "sweepT_surr_%d-%d.json" % (SI, SN))
    json.dump(dict(shifts=A.shifts, iters=A.iters, shard=A.shard, ntests_total=ntest,
                   threshold=thr, floor=floor, elapsed_s=el, errors=len(bad),
                   results=ok), open(out, "w"))
    print("saved %s  (%d results, %d errors)" % (out, len(ok), len(bad)))

    ok.sort(key=lambda r: r["p"])
    surv = [r for r in ok if r["p"] < thr]
    print("\n  survivors below %.2e : %d" % (thr, len(surv)))
    print("  pinned at the floor (ne=0) : %d" % len([r for r in ok if r["ne"] == 0]))
    byform = {}
    for r in surv:
        byform[r["form"]] = byform.get(r["form"], 0) + 1
    print("  survivors by form : %s" % (byform or "none"))
    print("\n  %-46s %-13s %10s %11s %5s" % ("combination", "form", "stat", "p", "ne"))
    for r in ok[:20]:
        print("  %-46s %-13s %10.4f %11.3e %5d"
              % (r["combo"][:46], r["form"], r["obs"], r["p"], r["ne"]))
