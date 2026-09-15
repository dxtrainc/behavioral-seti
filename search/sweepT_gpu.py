"""
THE TRIPLE SPACE ON THE GPU, AGAINST THE PAIRWISE-PRESERVING NULL.

The CPU version sustains ~580 draws/s across all 80 workers on zeus, which puts
the 2,934-test resid-triple enumeration at roughly six days. The A5000 does the
same work batched across draws -- thousands of independent surrogates of the
same three series at once -- and measured ~2,600 draws/s in float64 at
iters=12, ~5,200 at iters=6.

PRECISION. Everything stays in float64. float32 was measured and rejected: the
FFTs get 5.7x faster but they are only a quarter of the work, argsort barely
moves (19.7 -> 17.7 ms) because it is comparison- and bandwidth-bound rather
than FP64-ALU-bound, and the statistic shifts by 1.7e-05 -- real money in a
fourth moment.

THE STATISTIC. Computed from central moments rather than scipy.stats.kurtosis,
which has no batched form. Validated against the CPU implementation: agreement
to 6.7e-15, i.e. floating-point identical.

WHY NOT put_along_axis. cupy 13.6 has none, and the rank-map is a scatter.
Doing it by gather instead needs argsort(argsort(...)) -- two sorts, and the
sort is 73% of the runtime -- so it is done with flat indices, keeping it at one.

Checkpoint format, resume, sharding and the floor guard match sweepT_surr.py so
results from either are interchangeable.
"""
import argparse, json, hashlib, itertools, os, sys, time
import numpy as np
import cupy as cp
import cupyx.scipy.ndimage as cndi
from scipy.ndimage import median_filter

CH = os.path.expanduser("~/beacon_channels.npz")
_z = np.load(CH, allow_pickle=True)
M = _z["M"]
names = [str(x) for x in _z["names"]]
W = int(_z["W"])
N = M.shape[0]


# ---------------------------------------------------------------- CPU side, verbatim from sweepT
def prep(x):
    r = x / np.median(x)
    d = r - median_filter(r, size=W, mode="nearest")
    sd = d.std()
    return r, (d / sd if sd > 0 else d)


def gpu_prep(Y):
    """prep() over a batch of surrogates, on the GPU. (B,k,n) -> detrended residuals.

    Mirrors the CPU prep(): normalise by the median, subtract a 181-wide running
    median, divide by the standard deviation.
    """
    r = Y / cp.median(Y, axis=-1, keepdims=True)
    d = r - cndi.median_filter(r, size=(1, 1, W), mode="nearest")
    sd = d.std(axis=-1, keepdims=True)
    return d / cp.where(sd > 0, sd, 1.0)


def stable_seed(*parts):
    h = hashlib.blake2b(("|".join(str(p) for p in parts)).encode(), digest_size=8)
    return int.from_bytes(h.digest(), "big") % (2 ** 63)


# ---------------------------------------------------------------- GPU side
def gpu_stat(x):
    """Batched; x is (B, n). Matches sweepT.stat() to 6.7e-15."""
    d = cp.diff(x, axis=-1)
    med = cp.median(d, axis=-1, keepdims=True)
    s = 1.4826 * cp.median(cp.abs(d - med), axis=-1, keepdims=True)
    bad = (s <= 0).ravel()
    s = cp.where(s <= 0, 1.0, s)
    d = cp.clip(d, -3.0 * s, 3.0 * s)
    first = x[:, :1]
    y = cp.concatenate([first, first + cp.cumsum(d, axis=-1)], axis=-1)
    sd = y.std(axis=-1, keepdims=True)
    bad |= (sd <= 0).ravel()
    sd = cp.where(sd <= 0, 1.0, sd)
    y = (y - y.mean(axis=-1, keepdims=True)) / sd
    d2 = cp.diff(y, axis=-1)
    s2 = d2.std(axis=-1)
    bad |= (s2 <= 0)
    s2 = cp.where(s2 <= 0, 1.0, s2)
    dm = d2 - d2.mean(axis=-1, keepdims=True)
    m2 = (dm ** 2).mean(axis=-1)
    m4 = (dm ** 4).mean(axis=-1)
    out = m4 / (m2 ** 2) - 3.0 - 4.0 * cp.median(cp.abs(d2), axis=-1) / s2
    return cp.where(bad, cp.nan, out)


class GPUPlan:
    def __init__(self, X):
        self.k, self.n = X.shape
        self.srt = cp.asarray(np.sort(X, axis=1))
        F0 = np.fft.rfft(X, axis=1)
        self.Arel = cp.asarray(
            np.abs(F0) * np.exp(1j * (np.angle(F0) - np.angle(F0[0])[None, :])))
        self.nf = self.Arel.shape[1]

    def draw(self, B, rng, iters):
        n, k = self.n, self.k
        srtb = cp.ascontiguousarray(cp.broadcast_to(self.srt[None], (B, k, n)))
        srt_flat = srtb.ravel()
        off = (cp.arange(B * k, dtype=cp.int64) * n).reshape(B, k, 1)

        ph = rng.random((B, self.nf)) * (2.0 * np.pi)
        ph[:, 0] = 0.0
        if n % 2 == 0:
            ph[:, -1] = 0.0
        Y = cp.fft.irfft(self.Arel[None] * cp.exp(1j * ph)[:, None, :], n=n, axis=-1)

        def rankmap(Y):
            idx = cp.argsort(Y, axis=-1)
            flat = cp.empty(B * k * n, dtype=Y.dtype)
            flat[(idx + off).ravel()] = srt_flat
            return flat.reshape(B, k, n)

        for _ in range(iters):
            Y = rankmap(Y)
            phi0 = cp.angle(cp.fft.rfft(Y[:, 0, :], axis=-1))
            Y = cp.fft.irfft(self.Arel[None] * cp.exp(1j * phi0)[:, None, :], n=n, axis=-1)
        return rankmap(Y)


def run_test(i, j, k, nsh, iters, budget):
    """resid-triple: the product of the three detrended residuals."""
    ov = np.isfinite(M[i]) & np.isfinite(M[j]) & np.isfinite(M[k])
    n = int(ov.sum())
    raw = [M[c][ov] for c in (i, j, k)]
    P = [prep(x) for x in raw]
    U, V, Wv = P[0][1], P[1][1], P[2][1]
    obs = float(gpu_stat(cp.asarray((U * V * Wv)[None, :]))[0])
    if not np.isfinite(obs):
        return None

    # RAW path: the surrogate is drawn from the undetrended channels and prep()
    # is applied to it afterwards, exactly as it is applied to the data.
    plan = GPUPlan(np.vstack(raw))
    rng = cp.random.default_rng(stable_seed(i, j, k, "resid-triple"))
    B = max(32, min(2048, int(budget / (3 * n))))
    ne = nv = 0
    done = 0
    s1 = s2 = 0.0          # running sums for z; the null is never held whole
    while done < nsh:
        b = min(B, nsh - done)
        Y = plan.draw(b, rng, iters)
        D = gpu_prep(Y)
        g = D[:, 0, :] * D[:, 1, :] * D[:, 2, :]
        okg = cp.all(cp.isfinite(g), axis=-1)
        v = gpu_stat(g)
        ok = okg & cp.isfinite(v)
        vv = v[ok]
        ne += int(cp.count_nonzero(vv >= obs))
        nv += int(vv.size)
        s1 += float(cp.sum(vv))
        s2 += float(cp.sum(vv * vv))
        done += b
        del Y, D, g, v, ok, okg, vv
    if nv < nsh // 4:
        return None
    mean = s1 / nv
    var = max(s2 / nv - mean * mean, 0.0)      # population variance, as numpy std
    sd = var ** 0.5
    return dict(key="%d-%d-%d|resid-triple" % (i, j, k),
                combo="%s / %s / %s" % (names[i], names[j], names[k]),
                form="resid-triple", obs=obs, p=(1.0 + ne) / (1.0 + nv),
                ne=ne, nnull=nv, ndays=n,
                z=float((obs - mean) / max(sd, 1e-12)))


def hms(s):
    s = int(s)
    return "%dd %02dh %02dm" % (s // 86400, (s % 86400) // 3600, (s % 3600) // 60)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shifts", type=int, default=100000)
    ap.add_argument("--iters", type=int, default=6)
    ap.add_argument("--minov", type=int, default=2000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--shards", default="")
    ap.add_argument("--shardn", type=int, default=0)
    ap.add_argument("--budget", type=int, default=20_000_000,
                    help="elements per batch; batch size = budget/(3n)")
    ap.add_argument("--ntests", type=int, default=0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--ckpt", default="")
    ap.add_argument("--validate", default="",
                    help="jsonl of CPU results: rerun those exact tests and compare")
    A = ap.parse_args()
    if A.shards:
        MINE = sorted(int(x) for x in A.shards.split(","))
        SN = A.shardn or (max(MINE) + 1)
        SI = MINE[0]
    else:
        SI, SN = (int(x) for x in A.shard.split("/"))
        MINE = [SI]

    triples = []
    for i, j, k in itertools.combinations(range(N), 3):
        if (np.isfinite(M[i]) & np.isfinite(M[j]) & np.isfinite(M[k])).sum() >= A.minov:
            triples.append((i, j, k))
    ntest = len(triples)
    thr = A.alpha / max(ntest, 1)
    floor = 1.0 / (A.shifts + 1)

    print("\ntriple space on the GPU, resid-triple only")
    print("  %d tests, Bonferroni threshold %.3e" % (ntest, thr))
    print("  %d draws -> floor %.3e (%.1fx below threshold)" % (A.shifts, floor, thr / floor))
    if floor >= thr:
        raise SystemExit("  REFUSING: floor at or above threshold -- the run could not fire.")
    print("  iters=%d, float64, RAW surrogate path\n" % A.iters, flush=True)

    if A.validate:
        ref = {}
        for ln in open(A.validate):
            try:
                r = json.loads(ln)
                if r.get("form") == "resid-triple":
                    ref[r["key"]] = r
            except Exception:
                pass
        print("  validating %d CPU tests\n" % len(ref))
        print("  %-38s %11s %11s %9s | %10s %10s" %
              ("combination", "obs cpu", "obs gpu", "|d|", "p cpu", "p gpu"))
        dmax = 0.0
        for key, r in list(ref.items())[:A.ntests or len(ref)]:
            i, j, k = (int(x) for x in key.split("|")[0].split("-"))
            g = run_test(i, j, k, A.shifts, A.iters, A.budget)
            if g is None:
                print("  %-38s GPU returned None" % r["combo"][:38]); continue
            d = abs(g["obs"] - r["obs"])
            dmax = max(dmax, d)
            print("  %-38s %11.6f %11.6f %9.1e | %10.3e %10.3e"
                  % (r["combo"][:38], r["obs"], g["obs"], d, r["p"], g["p"]), flush=True)
        print("\n  worst |obs difference|: %.2e" % dmax)
        print("  %s" % ("OBS MATCHES" if dmax < 1e-9 else "*** OBS DIFFERS ***"))
        sys.exit(0)

    jobs = [t for n_, t in enumerate(triples) if (n_ % SN) in MINE]
    ckpt = A.ckpt or os.path.expanduser(
        "~/sweepT_gpu_%s-%d.jsonl" % ("-".join(str(x) for x in MINE), SN))
    done = set()
    if A.resume and os.path.exists(ckpt):
        for ln in open(ckpt, errors="replace"):
            try:
                done.add(json.loads(ln)["key"])
            except Exception:
                pass
        print("  resuming: %d already done" % len(done))
    todo = [t for t in jobs if "%d-%d-%d|resid-triple" % t not in done]
    if A.ntests:
        todo = todo[:A.ntests]
    print("  shards %s of %d: %d tests, %d to run"
          % (",".join(str(x) for x in MINE), SN, len(jobs), len(todo)))
    print("  checkpoint: %s\n" % ckpt, flush=True)

    t0 = time.time()
    fh = open(ckpt, "a", buffering=1)
    for c, (i, j, k) in enumerate(todo, 1):
        try:
            r = run_test(i, j, k, A.shifts, A.iters, A.budget)
        except Exception as e:
            r = dict(key="%d-%d-%d|resid-triple" % (i, j, k),
                     error="%s: %s" % (type(e).__name__, e))
            cp.get_default_memory_pool().free_all_blocks()
        if r:
            fh.write(json.dumps(r) + "\n")
        if c % 10 == 0:
            fh.flush(); os.fsync(fh.fileno())
            el = time.time() - t0
            rate = c / el
            print("  %5d/%-5d %5.1f%%  elapsed %s  eta %s  (%.0f draws/s)"
                  % (c, len(todo), 100.0 * c / len(todo), hms(el),
                     hms((len(todo) - c) / max(rate, 1e-12)), rate * A.shifts), flush=True)
    fh.flush(); os.fsync(fh.fileno()); fh.close()
    print("\ncompleted %d tests in %s" % (len(todo), hms(time.time() - t0)))
