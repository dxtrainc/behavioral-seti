"""
ROW 10: disc-integrated line-profile ratios, HARPS-N solar telescope.

Section 3.6 lists this as the designer's tenth channel -- level 1e-07, our reach
estimated at ~1e-06, margin ~0.1, NOT SEARCHED. The estimate has never been
measured, and every row in that table whose estimate was later measured came in
worse than the estimate. This measures it.

DATA. The HARPS-N solar telescope public release: 173,793 disc-integrated
spectra, 2015-07-15 to 2025-05-23, 9.86 years, 3,239 observing days, 300 s
exposures at 325 s cadence. Pulled through the DACE public API; no key needed.

CHANNELS, all dimensionless so instrument gain divides out -- the property that
makes this the designer's deepest channel:

  contrast     CCF depth as a fraction of the continuum. Already a ratio.
  bispan/fwhm  bisector span over line width, both velocities.
  sn_colour    S/N in a blue order over a red one: a flux ratio, so throughput
               and exposure time cancel. Physically a colour temperature.

Per-exposure relative precision is 7.96e-05 on contrast and on fwhm; bispan is
two orders worse (8.6e-03) and is expected to be the weak channel.

METHOD. Daily means (3,239 points) to escape the daytime-only window function
and its one-day aliases. Long-term variation -- the solar cycle -- is removed
with the same running median the rest of the paper uses. A sinusoid of
fractional amplitude d is injected; detection is the Lomb-Scargle peak over a
search band, against a null of circular shifts of the same series, which
preserves the sampling and the red noise exactly.

THE FIRST AMPLITUDE IS ZERO. Row 7's original injection reported 100% recovery
at its smallest amplitude because it had no zero arm and was measuring its own
residual. Not repeated.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from astropy.timeseries import LombScargle

W = 181            # the running-median window used throughout the paper


def build(d, chan, qc=True):
    num = lambda c: pd.to_numeric(d[c], errors="coerce")
    m = pd.Series(True, index=d.index)
    if qc and "spectro_drs_qc" in d:
        m &= d["spectro_drs_qc"].astype(bool)
    t = num("obj_date_bjd")
    if chan == "contrast":
        y = num("spectro_ccf_contrast")
    elif chan == "bis_over_fwhm":
        y = num("spectro_ccf_bispan") / num("spectro_ccf_fwhm")
    elif chan == "ew":
        # contrast x fwhm ~ equivalent width. Depth and width anticorrelate
        # under activity, so the product cancels much of the activity signal.
        y = num("spectro_ccf_contrast") * num("spectro_ccf_fwhm")
    elif chan == "sn_colour":
        y = num("spectro_flux_sn10") / num("spectro_flux_sn50")
    else:
        raise KeyError(chan)
    m &= t.notna() & y.notna() & np.isfinite(y)
    t, y = t[m].values, y[m].values
    day = np.floor(t).astype(np.int64)
    df = pd.DataFrame({"d": day, "y": y}).groupby("d")["y"].agg(["mean", "count"])
    df = df[df["count"] >= 5]                       # a day needs real coverage
    return df.index.values.astype(float), df["mean"].values


def prep(y):
    """Fractional variation with the slow component removed."""
    r = y / np.median(y)
    return r - median_filter(r, size=W, mode="nearest")


def stat(t, x, fmin, fmax, nf=4000):
    f = np.linspace(fmin, fmax, nf)
    return float(LombScargle(t, x).power(f).max())


class GPNull:
    """Gaussian surrogates with the data's empirical autocovariance.

    Built once per series: the Cholesky factor is reused for every draw, which
    is what makes a few hundred surrogates per test affordable.
    """

    def __init__(self, t, x, max_lag=120.0, nbin=60, jitter=1e-10):
        n = len(t)
        d = np.abs(t[:, None] - t[None, :])
        # empirical autocovariance, binned in lag
        edges = np.linspace(0.0, max_lag, nbin + 1)
        iu = np.triu_indices(n, 1)
        lag, prod = d[iu], (x[:, None] * x[None, :])[iu]
        acf = np.empty(nbin)
        for b in range(nbin):
            m = (lag >= edges[b]) & (lag < edges[b + 1])
            acf[b] = prod[m].mean() if m.any() else 0.0
        ctr = 0.5 * (edges[:-1] + edges[1:])
        var = float(np.mean(x * x))
        C = np.interp(d, ctr, acf, left=var, right=0.0)
        C[np.diag_indices(n)] = var
        # nearest positive-definite: clip negative eigenvalues
        w, V = np.linalg.eigh(C)
        w = np.clip(w, jitter * var, None)
        self.L = V * np.sqrt(w)
        self.n = n

    def draw(self, rng):
        return self.L @ rng.standard_normal(self.n)


def run(t, y, amp, f_inj, nsh, fmin, fmax, rng, drift=0.0, ks=None):
    r = y / np.median(y)
    if amp > 0:
        t0 = t - t.mean()
        # drift: fractional change in frequency across the whole span
        ph = 2 * np.pi * f_inj * (t0 + 0.5 * (drift / (t0.max() - t0.min())) * t0 ** 2)
        r = r * (1.0 + amp * np.sin(ph + rng.uniform(0, 2 * np.pi)))
    x = r - median_filter(r, size=W, mode="nearest")
    key = CFG.get("gpkey")
    gp = _GP_CACHE.get(key)
    if gp is None:
        # built from the UNINJECTED detrended series, once per worker
        r0 = y / np.median(y)
        gp = GPNull(t, r0 - median_filter(r0, size=W, mode="nearest"))
        _GP_CACHE[key] = gp
    if ks is None:
        obs = stat(t, x, fmin, fmax)
        null = [stat(t, gp.draw(rng), fmin, fmax) for _ in range(nsh)]
    else:
        obs = stat_chirp(t, x, fmin, fmax, ks)
        null = [stat_chirp(t, gp.draw(rng), fmin, fmax, ks) for _ in range(nsh)]
    a = np.asarray(null)
    ne = int((a >= obs).sum())
    return (1.0 + ne) / (1.0 + a.size)


CFG = {}
_GP_CACHE = {}


def _work(job):
    amp, trial = job
    rng = np.random.default_rng(7000 + trial * 97 + int(amp * 1e9))
    return run(CFG["t"], CFG["y"], amp, CFG["f_inj"], CFG["nsh"],
               CFG["fmin"], CFG["fmax"], rng, CFG.get("drift", 0.0), CFG.get("ks"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chan", default="contrast",
                    choices=["contrast", "ew", "bis_over_fwhm", "sn_colour"])
    ap.add_argument("--period", type=float, default=9.7, help="injected period, days")
    ap.add_argument("--shifts", type=int, default=400)
    ap.add_argument("--trials", type=int, default=24)
    ap.add_argument("--pool", type=int, default=60)
    ap.add_argument("--drift", type=float, default=0.0,
                    help="fractional frequency drift of the INJECTED signal across the span")
    ap.add_argument("--chirp", action="store_true",
                    help="search over frequency drift as well as frequency")
    ap.add_argument("--out", default="~/row10.json")
    A = ap.parse_args()

    d = pd.read_pickle("/home/dxtra/harpsn_sun.pkl")
    t, y = build(d, A.chan)
    print("\nROW 10 -- %s, HARPS-N solar" % A.chan)
    print("  %d daily points, %.2f yr, median gap %.1f d"
          % (len(t), (t.max() - t.min()) / 365.25, np.median(np.diff(np.sort(t)))))
    rel = float(np.std(prep(y)))
    print("  residual scatter after detrend: %.3e (fractional)" % rel)
    # search band: 4 to 40 days, clear of the 27 d rotation harmonics at the edges
    fmin, fmax = 1.0 / 40.0, 1.0 / 4.0
    print("  search band %.3f .. %.3f /d  (periods 4 .. 40 d), injected at %.1f d"
          % (fmin, fmax, A.period))
    print("  %d shifts per test, %d trials per amplitude\n" % (A.shifts, A.trials), flush=True)

    global CFG
    span = t.max() - t.min()
    ks = None
    if A.chirp:
        # drift rates spanning +/-3% of the frequency across the record
        ks = [d / span for d in np.linspace(-0.03, 0.03, 7)]
    CFG = dict(t=t, y=y, f_inj=1.0 / A.period, nsh=A.shifts, fmin=fmin, fmax=fmax,
               drift=A.drift, ks=ks, gpkey=A.chan)
    print("  injected drift %.1f%% across the span; search = %s"
          % (100 * A.drift, "de-chirped (7 drift rates)" if A.chirp else "fixed frequency"))
    amps = [0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]

    from multiprocessing import Pool
    print("  %10s %9s %9s   %s" % ("amplitude", "p<0.05", "p<0.01", "median p"))
    print("  " + "-" * 50)
    curve = []
    for a in amps:
        jobs = [(a, k) for k in range(A.trials)]
        with Pool(A.pool) as p:
            ps = np.array(p.map(_work, jobs))
        rec = float((ps < 0.05).mean())
        curve.append(dict(amp=a, rec05=rec, rec01=float((ps < 0.01).mean()),
                          med_p=float(np.median(ps))))
        print("  %10.1e %8.1f%% %8.1f%%   %.4f"
              % (a, 100 * rec, 100 * (ps < 0.01).mean(), np.median(ps)), flush=True)

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
    print("  designer level: 1e-07   ->  margin %s"
          % ("%.3f" % (1e-7 / a95) if a95 else "n/a"))
    print("  paper's tabulated estimate for row 10: reach ~1e-06, margin ~0.1")
    json.dump(dict(chan=A.chan, n=len(t), curve=curve, a95=a95,
                   zero_arm=r[0], monotonic=bool(mono)),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("\n  saved %s" % A.out)


if __name__ == "__main__":
    main()
