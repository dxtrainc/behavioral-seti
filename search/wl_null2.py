"""Exact-preservation variant: work in NORMAL-SCORE space.

Rank-transform each channel to exact standard normal ONCE, up front. Then a common
phase screen preserves: the marginal (a phase-randomised Gaussian is still Gaussian),
each power spectrum, and every pairwise cross-spectrum -- all three EXACTLY, by
construction rather than by annealing. No IAAFT iteration, no compromise.
Higher-order phase alignment is still destroyed. Cost: the per-channel marginal
non-Gaussianity is discarded, which single- and pair-searches already cover.
"""
import numpy as np, os, time
from scipy.stats import norm
import sys as _s, os as _o
_s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
from wl_null import xcorr_profile

def normal_score(x):
    r = np.argsort(np.argsort(x))
    return norm.ppf((r + 0.5) / len(x))

def common_phase(X, rng):
    k, n = X.shape
    F = np.fft.rfft(X, axis=1)
    ph = rng.uniform(0, 2*np.pi, F.shape[1]); ph[0] = 0.0
    if n % 2 == 0: ph[-1] = 0.0
    return np.fft.irfft(F * np.exp(1j*ph)[None, :], n=n, axis=1)

if __name__ == "__main__":
    M = np.load(os.path.expanduser("~/wl_M.npy"))
    names = open(os.path.expanduser("~/wl_names.txt")).read().split("\n")
    ix = {n: i for i, n in enumerate(names)}
    q = ["TSI", "F10.7", "cosmic ray", "X-ray bg"]
    ov = np.all(np.isfinite(M[[ix[n] for n in q]]), axis=0)
    X = np.vstack([normal_score(M[ix[n]][ov]) for n in q])
    print("quadruple:", " / ".join(q), " n =", X.shape[1])

    rng = np.random.default_rng(20260912)
    lags = [0,1,2,3,5,8,13,21,34,55]
    obs = xcorr_profile(X, lags)
    t0=time.time(); devs=[]; kurt=[]
    for t in range(200):
        S = common_phase(X, rng)
        devs.append(np.abs(xcorr_profile(S, lags) - obs))
        kurt.append([float(((s-s.mean())**4).mean()/s.var()**2 - 3.0) for s in S])
    el=time.time()-t0; devs=np.array(devs); kurt=np.array(kurt)

    print("\n=== EXACTNESS CHECK, 200 surrogates ===")
    print("  observed |xcorr| range      : %.4f .. %.4f"%(np.abs(obs).min(),np.abs(obs).max()))
    print("  mean |surrogate - observed| : %.3e"%devs.mean())
    print("  worst deviation             : %.3e"%devs.max())
    print("  as %% of observed rms        : %.2e %%"%(100*devs.mean()/np.sqrt((obs**2).mean())))
    print("\n  marginal excess kurtosis, observed  : %s"%np.round([float(((x-x.mean())**4).mean()/x.var()**2-3.0) for x in X],4))
    print("  marginal excess kurtosis, surrogate : %s +/- %s"%(np.round(kurt.mean(0),4),np.round(kurt.std(0),4)))
    print("\n  %.4f s per surrogate -> 15 keys x 10000 = %.2f core-hours (%.1f min on 88 cores)"
          %(el/200, el/200*15*10000/3600, el/200*15*10000/60/88))
