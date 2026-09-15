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


_NPZ = os.path.expanduser("~/beacon_channels.npz")
if os.path.exists(_NPZ):
    # Portable path: the assembled matrix, exported from zeus. Used on boxes
    # without the ingest, and identical to what the exec path produces.
    _z = np.load(_NPZ, allow_pickle=True)
    M = _z["M"]
    names = [str(x) for x in _z["names"]]
    N = M.shape[0]
    W = int(_z["W"])
    from scipy.ndimage import median_filter as _medfilt
    from scipy import stats as _sps

    def prep(x):
        r = x / np.median(x)
        d = r - _medfilt(r, size=W, mode="nearest")
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
        return float(_sps.kurtosis(d, fisher=True) - 4.0 * np.median(np.abs(d)) / s)

    def forms_tri(pa, pb, pc):
        ra, da = pa; rb, db = pb; rc, dc = pc
        out = {"resid-triple": da * db * dc}
        la = np.log(np.maximum(ra, 1e-12)); lb = np.log(np.maximum(rb, 1e-12))
        lc = np.log(np.maximum(rc, 1e-12))
        out["log-2nd-diff"] = la - 2.0 * lb + lc
        return {k: v for k, v in out.items() if np.all(np.isfinite(v))}
else:
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

        # RAW surrogate path. The surrogate is drawn from the undetrended
        # channels and prep() is applied to it afterwards, exactly as it is
        # applied to the data. Surrogating the prepared residuals instead would
        # be cheaper -- it keeps the 181-wide median filter out of this loop --
        # but IAAFT does not commute with detrending, so that is a DIFFERENT
        # null, and not the one the zero arm validated.
        plan = Plan(np.vstack([M[i][ov], M[j][ov], M[k][ov]]))
        rng = np.random.default_rng(seed)
        null = []
        for _ in range(nsh):
            Q = [prep(r) for r in plan.draw(rng, iters)]
            fb = forms_tri(*Q)
            if fname not in fb:
                continue
            v = stat(fb[fname])
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
    ap.add_argument("--shards", default="",
                    help="comma list of shard indices, e.g. 0,1,2 -- used with --shardn")
    ap.add_argument("--shardn", type=int, default=0,
                    help="total shards when --shards is given")
    ap.add_argument("--forms", default="resid-triple,log-2nd-diff",
                    help="comma-separated forms to enumerate. The zero arm shows\n                          log-2nd-diff cannot fire under this null, so the default\n                          for a real run is resid-triple alone.")
    ap.add_argument("--sample", type=int, default=0,
                    help="time N triples spread across the overlap range and stop")
    ap.add_argument("--ntests", type=int, default=0, help="0 = all; >0 truncates, for testing")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--ckpt", default="")
    A = ap.parse_args()
    MINOV = A.minov
    if A.shards:
        MINE = sorted(int(x) for x in A.shards.split(","))
        SN = A.shardn or (max(MINE) + 1)
        SI = MINE[0]
    else:
        SI, SN = (int(x) for x in A.shard.split("/"))
        MINE = [SI]

    want = tuple(f.strip() for f in A.forms.split(",") if f.strip())
    jobs = []
    for i, j, k in itertools.combinations(range(N), 3):
        ov = np.isfinite(M[i]) & np.isfinite(M[j]) & np.isfinite(M[k])
        if ov.sum() < MINOV:
            continue
        if "resid-triple" in want:
            jobs.append((i, j, k, "resid-triple"))
        if "log-2nd-diff" in want:
            # log-2nd-diff is not symmetric: each channel takes a turn as the middle
            jobs.append((i, j, k, "log-2nd-diff"))
            jobs.append((j, i, k, "log-2nd-diff"))
            jobs.append((i, k, j, "log-2nd-diff"))
    ntest = len(jobs)
    thr = A.alpha / max(ntest, 1)
    floor = 1.0 / (A.shifts + 1)

    print("\ntriple space against the pairwise-preserving null")
    print("  %d tests, forms: %s" % (ntest, ", ".join(want)))
    print("  Bonferroni threshold %.3e" % thr)
    print("  %d draws -> empirical floor %.3e" % (A.shifts, floor))
    if floor >= thr:
        raise SystemExit(
            "  REFUSING TO RUN: floor at or above threshold -- no test could fire\n"
            "  whatever the data contain (the section 5.2 trap).\n"
            "  Need --shifts > %d." % int(np.ceil(1.0 / thr)))
    print("  floor is %.1fx below the threshold -- the run can fire" % (thr / floor))

    jobs = [x for n_, x in enumerate(jobs) if (n_ % SN) in MINE]
    ckpt = A.ckpt or os.path.join(
        os.path.expanduser("~"),
        "sweepT_surr_%s-%d.jsonl" % ("-".join(str(x) for x in MINE), SN))
    done = {}
    if A.resume and os.path.exists(ckpt):
        for ln in open(ckpt, errors="replace"):
            try:
                r = json.loads(ln)
                done[r["key"]] = r
            except Exception:
                continue
        print("  resuming: %d tests already in %s" % (len(done), ckpt))

    if A.sample:
        # UNIFORM RANDOM sample: representative of the real overlap distribution
        # by construction. An even spread across the sorted range over-weights
        # long triples, and the count must exceed --pool or idle workers are
        # counted as though they were working.
        rsel = np.random.default_rng(20260914)
        pick = rsel.choice(len(jobs), size=min(A.sample, len(jobs)), replace=False)
        jobs = [jobs[int(p)] for p in pick]
        ns = [int((np.isfinite(M[a]) & np.isfinite(M[b]) & np.isfinite(M[c])).sum())
              for (a, b, c, _f) in jobs]
        print("  SAMPLE MODE: %d random triples, n %d..%d (mean %d), pool %d"
              % (len(jobs), min(ns), max(ns), int(np.mean(ns)), A.pool))
        if len(jobs) < A.pool:
            print("  WARNING: fewer tests than workers -- the rate will be understated")
    todo = [(i, j, k, fn, A.shifts, A.iters, abs(hash((i, j, k, fn))) % (2 ** 32))
            for (i, j, k, fn) in jobs if "%d-%d-%d|%s" % (i, j, k, fn) not in done]
    if A.ntests:
        todo = todo[:A.ntests]
    print("  shards %s of %d: %d tests, %d to run, pool=%d"
          % (",".join(str(x) for x in MINE), SN, len(jobs), len(todo), A.pool))
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
    out = os.path.join(os.path.expanduser("~"), "sweepT_surr_%s-%d.json"
                       % ("-".join(str(x) for x in MINE), SN))
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
