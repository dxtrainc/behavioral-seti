"""
ENVELOPE-PRESERVING SURROGATE.

WHY. The common-phase null of search/wl_sweep.py calibrates perfectly (1.02x, 0/45
uniformity failures) and is still rejected by real solar data at median z = 45.3.
The rejection is the Sun's own higher-order structure. But "higher-order structure"
is not one thing, and the dominant part of it is cheap to name: solar activity is
AMPLITUDE MODULATED. Variance rises and falls with the cycle, so the series is not
stationary, and a surrogate that assumes stationarity destroys that modulation and
looks nothing like the data -- for reasons that have nothing to do with four-way
alignment.

So before concluding the space is untestable, separate the two hypotheses:
    (a) the Sun has irreducible four-way structure  -> the space really is void
    (b) the null wrongly assumed stationarity       -> the space reopens

CONSTRUCTION. Estimate a slowly varying robust scale s_k(t) per channel. Standardise
u_k = x_k / s_k, which is approximately stationary. Apply the common phase screen to
u -- preserving every pairwise cross-spectrum of the STANDARDISED series and
destroying higher-order alignment as before -- then multiply the envelope back.
The surrogate then carries the data's own amplitude modulation exactly.

THE WINDOW IS A SENSITIVITY LIMIT, NOT A FREE PARAMETER. Any modulation slower than
the envelope window is absorbed into the envelope and becomes invisible. W is
therefore scanned rather than chosen, and the blind band is stated with any result.

CALENDAR SYMMETRY, after the lesson in validate/wl_cal4.py: everything is done on
the full daily grid and subsampled through the real 1,064-break mask at the end, so
data and surrogate share a calendar and not merely a spectrum.
"""
import io, os, sys, time
import numpy as np
from scipy.stats import norm, kstest
from scipy.ndimage import median_filter
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "search"))
import wl_sweep as W

FULL = np.arange(len(W.OV)); SUB = np.where(W.OV)[0]; NF = len(FULL)

def _fill(ch):
    v = W.M[W.ix[ch]].astype(float); m = np.isfinite(v)
    return np.interp(FULL, FULL[m], v[m])

def _ns(x):
    return norm.ppf((np.argsort(np.argsort(x)) + 0.5) / len(x))

XF = np.vstack([_ns(_fill(c)) for c in W.WL])          # full grid, normal score

def envelope(Y, win):
    """robust slowly varying scale, per channel"""
    d = Y - median_filter(Y, size=(1, win), mode="nearest")
    s = 1.4826 * median_filter(np.abs(d), size=(1, win), mode="nearest")
    return np.maximum(s, 1e-6), d

def env_surrogate(Y, win, rng, cache=None):
    s, d = cache if cache is not None else envelope(Y, win)
    u = d / s
    F = np.fft.rfft(u, axis=1)
    ph = rng.uniform(0, 2*np.pi, F.shape[1]); ph[0] = 0.0
    if NF % 2 == 0: ph[-1] = 0.0
    us = np.fft.irfft(F * np.exp(1j*ph)[None, :], n=NF, axis=1)
    return us * s                                       # envelope re-imposed

def zscores(win, nsur, seed=5150):
    rng = np.random.default_rng(seed)
    cache = envelope(XF, win)
    base = (cache[1] / cache[0]) * cache[0]             # detrended data, same pipeline
    obs = W.all_stats(base[:, SUB])
    S = np.array([W.all_stats(env_surrogate(XF, win, rng, cache)[:, SUB])
                  for _ in range(nsur)])
    return (obs - np.nanmean(S, 0)) / np.nanstd(S, 0)

def _cal_one(arg):
    win, seed = arg
    rng = np.random.default_rng(seed)
    cache = envelope(XF, win)
    YF = env_surrogate(XF, win, rng, cache)             # no 4-way structure, real envelope
    c2 = envelope(YF, win)
    obs = W.all_stats(((c2[1]/c2[0])*c2[0])[:, SUB])
    S = np.array([W.all_stats(env_surrogate(YF, win, rng, c2)[:, SUB]) for _ in range(200)])
    return np.array([(1.0+(S[:, t] >= obs[t]).sum())/201.0 for t in range(len(W.TESTS))])

if __name__ == "__main__":
    WINS = [91, 365, 1461]
    print("envelope-preserving surrogate; blind to modulation slower than the window\n")
    print("=== 1. does real data still reject the null? ===")
    print("  %6s %8s %9s %9s %9s %9s" % ("window", "median z", "max z", "|z|>20", "3rd-diff", "quad-prod"))
    for win in WINS:
        z = zscores(win, 300)
        zq = np.array([z[t] for t, (q, f) in enumerate(W.TESTS) if f == "quad-prod"])
        z3 = np.array([z[t] for t, (q, f) in enumerate(W.TESTS) if f == "3rd-diff"])
        print("  %5dd %8.1f %9.1f %9d %9.1f %9.1f"
              % (win, np.nanmedian(z), np.nanmax(z), (np.abs(z) > 20).sum(),
                 np.nanmedian(z3), np.nanmedian(zq)))
    print("\n  (common-phase, no envelope, for comparison: median z 45.3, max 426.9, 30/45)")

    print("\n=== 2. is the envelope null itself calibrated? ===")
    NT = 150
    for win in WINS:
        t0 = time.time()
        with Pool(80) as p:
            P = np.array(p.map(_cal_one, [(win, 6000+i) for i in range(NT)], chunksize=1))
        bad = sum(1 for t in range(len(W.TESTS))
                  if kstest(P[:, t][np.isfinite(P[:, t])], "uniform").pvalue < 0.05/len(W.TESTS))
        allp = P[np.isfinite(P)]
        print("  %5dd  mean p %.4f  frac<0.05 %.4f  ratio %.2fx  failing %d/45  (%.0fs)"
              % (win, allp.mean(), (allp < 0.05).mean(), (allp < 0.05).mean()/0.05,
                 bad, time.time()-t0))
