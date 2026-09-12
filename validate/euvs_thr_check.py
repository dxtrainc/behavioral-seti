"""
THE DECISIVE ONE: does the threshold hold on the channel behind the headline limit?

Section 4.2's strongest claim is 1.4e-6 in fractional Lyman-alpha at two minutes, from
GOES-16 EUVS. The band-wide false-alarm control has now failed on two channels -- MAVEN's
Mars line at 15.2x in threshold, GOES XRS-B at 7.10x -- but NEITHER IS THE CHANNEL THAT
PRODUCED THE HEADLINE. XRS-B is flare-dominated and clipped at 4 sigma before the
spectrum, which truncates the very tail that causes the problem; EUVS is not.

So this measures the family-wise rate on EUVS Lyman-alpha itself, with the same empirical
mu = median(R)/ln 2 the paper uses, and reads the empirical threshold off the
distribution of the maximum of R across surrogates.

If the rate is near 0.05 the headline stands and only the X-ray rows move. If it is in
the tens or hundreds, the paper's strongest single number requotes.
"""
import os, sys, json, time
import numpy as np
from scipy.ndimage import median_filter

NSUR = int(os.environ.get("NSUR", "200"))
D = os.path.expanduser("~")

def prep(x, wmin=1440):
    """EUVS Lyman-alpha: log irradiance, slow trend removed, gaps zeroed.
    No flare clip -- unlike XRS-B, which is why this channel is the one to test."""
    y = np.full(len(x), np.nan)
    ok = np.isfinite(x) & (x > 0)
    y[ok] = np.log(x[ok])
    med = np.nanmedian(y)
    v = np.nan_to_num(y, nan=med); m = ok.astype(np.float64)
    cv = np.concatenate([[0], np.cumsum(v)]); cm = np.concatenate([[0], np.cumsum(m)])
    h = wmin//2; i = np.arange(len(y))
    lo = np.maximum(0, i-h); hi = np.minimum(len(y), i+h)
    den = cm[hi]-cm[lo]
    trend = np.where(den > 10, (cv[hi]-cv[lo])/np.maximum(den, 1), med)
    return np.where(ok, y-trend, 0.0), ok

def spectrum(r):
    n = len(r)//2*2
    F = np.fft.rfft(r[:n]*np.hanning(n))
    return np.fft.rfftfreq(n, d=60.0)[1:], (np.abs(F)**2)[1:]

def cont(P, w=801):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

if __name__ == "__main__":
    t0 = time.time()
    z = np.load(os.path.join(D, "euvs1m_g16.npz"))
    x = z["irr_1216"].astype(float)                 # Lyman-alpha, the headline channel
    r, ok = prep(x)
    f, P = spectrum(r)
    R = P/cont(P)
    mu = np.median(R)/np.log(2.0)
    N = len(R)
    thr = mu*np.log(N/0.05)
    print("GOES-16 EUVS Lyman-alpha (121.6 nm): %d samples, %.1f%% present"
          % (len(x), 100*ok.mean()))
    print("  %d bins, mu = %.4f, analytic threshold = %.2f\n" % (N, mu, thr), flush=True)

    n = len(r); rng = np.random.default_rng(2718)
    counts, maxes = [], []
    for i in range(NSUR):
        sur = np.roll(r, int(rng.integers(1000, n-1000))).copy(); sur[~ok] = 0.0
        _, P2 = spectrum(sur)
        R2 = P2/cont(P2)
        counts.append(int((R2 > thr).sum())); maxes.append(float(R2.max()))
    counts = np.array(counts); maxes = np.array(maxes)
    emp = float(np.quantile(maxes, 0.95))
    print("  BAND-WIDE FALSE ALARMS per surrogate over %d bins" % N)
    print("    measured  %.2f      expected  0.05      ratio %.0fx"
          % (counts.mean(), counts.mean()/0.05))
    print("\n  EMPIRICAL threshold, 95th pct of max-R   %8.2f" % emp)
    print("  analytic threshold                       %8.2f" % thr)
    fac = np.sqrt(max(emp/thr, 1.0))
    print("  analytic is %.2fx too permissive; limits %.2fx looser" % (emp/thr, fac))
    print("\n  HEADLINE LIMIT")
    print("    section 4.2 as printed   1.4e-06")
    print("    requoted                 %.1e" % (1.4e-6*fac))
    print("\n  %s" % ("--> the headline STANDS" if counts.mean() <= 0.5 else
                      "--> the headline REQUOTES"))
    json.dump(dict(mu=float(mu), n_bins=int(N), analytic=float(thr), empirical=emp,
                   mean_false_alarms=float(counts.mean()), factor=float(fac),
                   requoted=float(1.4e-6*fac)),
              open(os.path.join(D, "euvs_thr_check.json"), "w"), indent=1)
    print("  %.0f s" % (time.time()-t0))
