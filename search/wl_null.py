"""
THE MATCHED HIGHER-ORDER NULL.

Section 4.9 voided the triple sweep because the circular shift destroys pairwise
structure as well as three-way structure. The surrogate must hold every pairwise
cross-spectrum FIXED and destroy only the higher-order alignment.

CONSTRUCTION. Give every channel the SAME random phase screen phi(f):
    X_k(f) -> X_k(f) exp(i phi(f))
The cross-spectrum S_ij(f) = X_i X_j* picks up exp(i phi) exp(-i phi) = 1, so every
pairwise cross-spectrum -- hence every pairwise cross-correlation at every lag -- is
preserved EXACTLY. Each power spectrum |X_k|^2 is untouched. But the bispectrum
picks up exp(i[phi(f1)+phi(f2)-phi(f1+f2)]), which is random unless phi is linear in
f, so three- and four-way phase alignment is destroyed. That is precisely the
decomposition the void called for.

THE CATCH. Phase randomisation Gaussianises the marginal, and our statistic is a
KURTOSIS of increments -- a non-Gaussianity measure. A Gaussianised surrogate would
have a wildly wrong null. So the phase screen is applied inside an IAAFT loop that
rank-remaps each channel back to its own empirical marginal each iteration. The
remap perturbs the cross-spectra slightly; how much is measured below, not assumed.
"""
import numpy as np, os, time

def common_phase_iaaft(X, rng, iters=60):
    """X: (k, n) real. Returns surrogate with per-channel marginal + spectrum and
    ALL pairwise cross-spectra preserved, higher-order alignment destroyed."""
    k, n = X.shape
    amp  = np.abs(np.fft.rfft(X, axis=1))          # target per-channel amplitude
    srt  = np.sort(X, axis=1)                      # target per-channel marginal
    # common random phase screen
    ph = rng.uniform(0, 2*np.pi, amp.shape[1]); ph[0] = 0.0
    if n % 2 == 0: ph[-1] = 0.0
    Y = np.fft.irfft(np.fft.rfft(X, axis=1) * np.exp(1j*ph)[None, :], n=n, axis=1)
    for _ in range(iters):
        # 1. impose each channel marginal by rank remap
        for r in range(k):
            Y[r, np.argsort(Y[r])] = srt[r]
        # 2. impose each channel amplitude, keeping the phases Y now has
        F = np.fft.rfft(Y, axis=1)
        F = amp * np.exp(1j*np.angle(F))
        Y = np.fft.irfft(F, n=n, axis=1)
    for r in range(k):
        Y[r, np.argsort(Y[r])] = srt[r]
    return Y

def xcorr_profile(X, lags):
    """all pairwise cross-correlations at the given lags -> flat vector"""
    k, n = X.shape
    Z = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)
    out = []
    for i in range(k):
        for j in range(i+1, k):
            for L in lags:
                out.append(float(np.mean(Z[i] * np.roll(Z[j], L))))
    return np.array(out)

if __name__ == "__main__":
    M = np.load(os.path.expanduser("~/wl_M.npy"))
    names = open(os.path.expanduser("~/wl_names.txt")).read().split("\n")
    ix = {n: i for i, n in enumerate(names)}
    q = ["TSI", "F10.7", "cosmic ray", "X-ray bg"]
    rows = [ix[n] for n in q]
    ov = np.all(np.isfinite(M[rows]), axis=0)
    X = M[rows][:, ov]
    print("quadruple:", " / ".join(q), " n =", X.shape[1])

    # gap structure on the overlap -- FFT assumes even sampling
    idx = np.where(ov)[0]; gaps = np.diff(idx)
    print("gaps: %d breaks, max %d d, %.2f%% of steps > 1 d"
          % ((gaps > 1).sum(), gaps.max(), 100.0*(gaps > 1).mean()))

    X = np.log(np.maximum(X / np.median(X, axis=1, keepdims=True), 1e-12))
    rng = np.random.default_rng(20260912)
    lags = [0, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    obs = xcorr_profile(X, lags)

    t0 = time.time()
    devs = []
    for t in range(20):
        S = common_phase_iaaft(X, rng)
        devs.append(np.abs(xcorr_profile(S, lags) - obs))
    el = time.time() - t0
    devs = np.array(devs)
    print("\n=== pairwise cross-correlation preservation, 20 surrogates ===")
    print("  observed |xcorr| range      : %.4f .. %.4f" % (np.abs(obs).min(), np.abs(obs).max()))
    print("  mean |surrogate - observed| : %.5f" % devs.mean())
    print("  worst deviation             : %.5f" % devs.max())
    print("  as %% of observed rms        : %.2f%%" % (100*devs.mean()/np.sqrt((obs**2).mean())))
    print("\n  %.3f s per surrogate  ->  15 keys x 10000 surrogates = %.1f core-hours"
          % (el/20, el/20*15*10000/3600))
