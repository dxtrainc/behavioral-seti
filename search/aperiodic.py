"""
APERIODIC MORPHOLOGIES: the class the morphology suite could not bound.

Section 4.20 measures what the sinusoid assumption costs across six *strictly
periodic* shapes and finds the worst case is 3.15x. That bound rests on a
structural fact -- any periodic modulation deposits a comb of lines, and a
max-over-band search finds the strongest line -- and it therefore says nothing
about a signal that deposits no comb at all: a one-off transient, a
non-repeating code, or a carrier whose frequency wanders non-linearly.

This measures that class two ways, because the question has two halves:

  1. What does the paper's OWN search (Lomb-Scargle peak over the 4-40 d band)
     cost on an aperiodic signal? That is the size of the blind spot in every
     limit quoted so far.
  2. What would a search actually matched to such a signal achieve? A
     multi-scale matched filter -- max normalized response over arrival time
     AND duration -- is the morphology-agnostic detector for an isolated
     excursion of unknown length, and is what a burst search uses.

Running the 2x2 (and its aperiodic extensions) separates "the signal is
undetectable in this data" from "our statistic is the wrong one".

CONVENTION. Every waveform is scaled to the RMS of a unit-amplitude sine, the
same power-budget convention as the periodic suite, so the thresholds here are
directly comparable to the 4.20 table.

THE DETRENDING HAS A PASSBAND. The running median is W=181 d, so an excursion
much broader than that is removed before any statistic sees it, and one shorter
than the 1 d cadence is unresolved. Transient widths are reported against that
window rather than quoted as though the search were flat in duration.

THE FIRST AMPLITUDE IS ZERO, as everywhere else in this paper.
"""
import argparse, json, os, sys
import numpy as np
from scipy.signal import fftconvolve
from scipy.ndimage import median_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import row10 as R
import pandas as pd

SCALES = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)   # matched-filter durations, days
EDGE = 200                                    # grid points trimmed from each end


def _kernels():
    ks = []
    for s in SCALES:
        half = int(np.ceil(4 * s))
        u = np.arange(-half, half + 1, dtype=float)
        ks.append(np.exp(-0.5 * (u / s) ** 2))
    return ks


_K = _kernels()


def gridmap(t):
    """Daily grid indices. The cadence IS daily, so this is a relabelling."""
    idx = np.rint(t - np.floor(t.min())).astype(np.int64)
    return idx, int(idx.max()) + 1


def stat_transient(idx, n, x):
    """Max normalized matched-filter response over arrival time and duration.

    Applied identically to the data and to every surrogate, so the gridding is
    shared by both and cannot bias the calibration -- the null is still drawn
    at the observed times by GPNull.
    """
    xg = np.zeros(n); mg = np.zeros(n)
    xg[idx] = x; mg[idx] = 1.0
    sc = 1.4826 * float(np.median(np.abs(x - np.median(x))))
    if sc <= 0:
        sc = float(x.std())
    if sc <= 0:
        return 0.0
    lo, hi = (EDGE, n - EDGE) if n > 2 * EDGE + 10 else (0, n)
    best = 0.0
    for k in _K:
        num = fftconvolve(xg, k, mode="same")
        nrm = np.sqrt(np.clip(fftconvolve(mg, k * k, mode="same"), 0.0, None))
        good = nrm > 1e-12
        z = np.zeros(n)
        z[good] = num[good] / (sc * nrm[good])
        m = float(np.abs(z[lo:hi]).max())
        if m > best:
            best = m
    return best


def stat_transient_ratio(idx, n, x):
    """Largest excursion measured against the record's OWN population of them.

    max|z| alone is an extreme-value statistic, so against a Gaussian null it
    reports the Sun's non-Gaussianity rather than any injected transient: the
    residual carries excess kurtosis 1.58 and max|x|/MAD 6.74 against 3.41 for
    a Gaussian surrogate of the same autocovariance, putting the unmodified
    data 4-6.4 sd above that null at every scale. Matching the marginal (AAFT)
    removes most but not all of it, because real solar excursions also cluster
    in time.

    Dividing the peak by the 99th percentile of the same response makes the
    statistic scale-free and asks the question a single record can answer --
    is the biggest excursion bigger than this record's other big excursions?
    With AAFT draws that calibrates to within 1.2 sd at every scale.
    """
    xg = np.zeros(n); mg = np.zeros(n)
    xg[idx] = x; mg[idx] = 1.0
    sc = 1.4826 * float(np.median(np.abs(x - np.median(x))))
    if sc <= 0:
        sc = float(x.std())
    if sc <= 0:
        return 0.0
    lo, hi = (EDGE, n - EDGE) if n > 2 * EDGE + 10 else (0, n)
    best = 0.0
    for k in _K:
        num = fftconvolve(xg, k, mode="same")
        nrm = np.sqrt(np.clip(fftconvolve(mg, k * k, mode="same"), 0.0, None))
        good = nrm > 1e-12
        z = np.zeros(n)
        z[good] = num[good] / (sc * nrm[good])
        az = np.abs(z[lo:hi])
        q = float(np.quantile(az, 0.99))
        if q > 0:
            m = float(az.max()) / q
            if m > best:
                best = m
    return best


def modulation_ap(kind, t, rng, width=30.0, chip=9.7, f=1.0 / 9.7, wander=0.10):
    """Zero-mean modulation with the RMS of a unit-amplitude sine."""
    if kind == "sine":
        w = np.sin(2 * np.pi * f * (t - t.mean()) + rng.uniform(0, 2 * np.pi))
    elif kind == "transient":
        tc = rng.uniform(t.min() + 250.0, t.max() - 250.0)
        w = np.exp(-0.5 * ((t - tc) / width) ** 2)
    elif kind == "code":
        c = np.floor((t - t.min()) / chip).astype(np.int64)
        signs = rng.choice(np.array([-1.0, 1.0]), size=int(c.max()) + 1)
        w = signs[c]
    elif kind == "wander":
        o = np.argsort(t); ts = t[o]
        dt = np.diff(ts, prepend=ts[0])
        d = np.cumsum(rng.standard_normal(ts.size))
        sd = float(d.std())
        d = d / (sd if sd > 0 else 1.0) * wander
        ph = 2 * np.pi * f * np.cumsum(dt * (1.0 + d))
        ws = np.sin(ph + rng.uniform(0, 2 * np.pi))
        w = np.empty_like(ws); w[o] = ws
    else:
        raise ValueError("unknown morphology %r" % kind)
    w = w - w.mean()
    s = float(w.std())
    return w / (s * np.sqrt(2.0)) if s > 0 else w


CFG = {}
_GP = {}


def _gp(t, y):
    g = _GP.get("gp")
    if g is None:
        r0 = y / np.median(y)
        g = R.GPNull(t, r0 - median_filter(r0, size=R.W, mode="nearest"))
        _GP["gp"] = g
    return g


def _work(job):
    amp, trial = job
    rng = np.random.default_rng(9100 + trial * 131 + int(amp * 1e9))
    t, y = CFG["t"], CFG["y"]
    r = y / np.median(y)
    if amp > 0:
        r = r * (1.0 + amp * modulation_ap(CFG["wave"], t, rng,
                                           width=CFG["width"], f=CFG["f_inj"]))
    x = r - median_filter(r, size=R.W, mode="nearest")
    gp = _gp(t, y)
    idx, n = CFG["idx"], CFG["n"]
    if CFG["stat"] == "psd":
        fn = lambda v: R.stat(t, v, CFG["fmin"], CFG["fmax"])
    elif CFG["stat"] == "tratio":
        fn = lambda v: stat_transient_ratio(idx, n, v)
    else:
        fn = lambda v: stat_transient(idx, n, v)
    if CFG.get("null") == "aaft":
        xs = CFG["xs"]
        draw = lambda: xs[np.argsort(np.argsort(gp.draw(rng)))]
    else:
        draw = lambda: gp.draw(rng)
    obs = fn(x)
    null = np.asarray([fn(draw()) for _ in range(CFG["nsh"])])
    return (1.0 + int((null >= obs).sum())) / (1.0 + null.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chan", default="ew",
                    choices=["contrast", "ew", "bis_over_fwhm", "sn_colour"])
    ap.add_argument("--wave", default="transient",
                    choices=["sine", "transient", "code", "wander"])
    ap.add_argument("--stat", default="tratio", choices=["psd", "trans", "tratio"])
    ap.add_argument("--null", default="aaft", choices=["gp", "aaft"])
    ap.add_argument("--width", type=float, default=30.0,
                    help="transient duration (Gaussian sigma), days")
    ap.add_argument("--period", type=float, default=9.7)
    ap.add_argument("--shifts", type=int, default=300)
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--pool", type=int, default=12)
    ap.add_argument("--amps", default="")
    ap.add_argument("--out", default="~/aperiodic.json")
    A = ap.parse_args()

    d = pd.read_pickle("/home/dxtra/harpsn_sun.pkl")
    t, y = R.build(d, A.chan)
    idx, n = gridmap(t)
    fmin, fmax = 1.0 / 40.0, 1.0 / 4.0

    print("\nAPERIODIC -- %s, HARPS-N solar" % A.chan)
    print("  %d daily points, %.2f yr; grid %d d, %d trimmed each end"
          % (len(t), (t.max() - t.min()) / 365.25, n, EDGE))
    print("  morphology %s   statistic %s"
          % (A.wave, "Lomb-Scargle peak, 4-40 d (the paper's own search)"
             if A.stat == "psd" else "matched filter %s, %s d"
             % ("peak/q99" if A.stat == "tratio" else "peak", list(map(int, SCALES)))))
    print("  null: %s" % ("rank-matched (AAFT)" if A.null == "aaft" else "Gaussian (GP)"))
    if A.wave == "transient":
        print("  transient sigma %.0f d  (detrend window W=%d d)" % (A.width, R.W))
    print("  %d surrogates per test, %d trials per amplitude" % (A.shifts, A.trials), flush=True)

    global CFG
    r0 = y / np.median(y)
    x0 = r0 - median_filter(r0, size=R.W, mode="nearest")
    CFG = dict(t=t, y=y, idx=idx, n=n, fmin=fmin, fmax=fmax, nsh=A.shifts,
               wave=A.wave, stat=A.stat, width=A.width, f_inj=1.0 / A.period,
               null=A.null, xs=np.sort(x0))

    amps = ([0.0] + [float(v) for v in A.amps.split(",") if v.strip()]
            if A.amps else [0.0, 4e-5, 7e-5, 1.2e-4, 2e-4, 3.5e-4, 6e-4,
                            1e-3, 1.7e-3, 3e-3, 5e-3, 8.5e-3, 1.4e-2])

    from multiprocessing import Pool
    print("\n  %10s %9s %9s   %s" % ("amplitude", "p<0.05", "p<0.01", "median p"))
    print("  " + "-" * 50)
    curve = []
    pool = Pool(A.pool)
    try:
        for a in amps:
            ps = np.array(pool.map(_work, [(a, k) for k in range(A.trials)]))
            rec = float((ps < 0.05).mean())
            curve.append(dict(amp=a, rec05=rec, rec01=float((ps < 0.01).mean()),
                              med_p=float(np.median(ps))))
            print("  %10.1e %8.1f%% %8.1f%%   %.4f"
                  % (a, 100 * rec, 100 * (ps < 0.01).mean(), np.median(ps)), flush=True)
    finally:
        pool.close(); pool.join()

    r = [c["rec05"] for c in curve]
    mono = all(r[i] >= r[i - 1] - 0.15 for i in range(1, len(r)))
    a95 = None
    for i in range(1, len(r)):
        if r[i] >= 0.95 > r[i - 1]:
            f = (0.95 - r[i - 1]) / max(r[i] - r[i - 1], 1e-9)
            a95 = float(np.exp(np.log(max(amps[i - 1], 1e-12)) +
                               f * (np.log(amps[i]) - np.log(max(amps[i - 1], 1e-12)))))
            break
    print("\n  zero arm      : %.1f%%  %s" % (100 * r[0],
          "(correct)" if r[0] <= 0.10 else "*** FIRES ON NOTHING ***"))
    print("  monotonic     : %s" % mono)
    print("  95%% recovery  : %s" % ("%.2e" % a95 if a95 else "not reached on this grid"))
    json.dump(dict(chan=A.chan, wave=A.wave, stat=A.stat, null=A.null, width=A.width,
                   n=len(t), curve=curve, a95=a95, zero_arm=r[0],
                   monotonic=bool(mono)),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("\n  saved %s" % A.out)


if __name__ == "__main__":
    main()
