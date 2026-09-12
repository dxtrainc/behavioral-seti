"""
DOES SECTION 4.2's THRESHOLD HOLD ITS NOMINAL FALSE-ALARM RATE ON GOES?

Section 4.2 states that continuum-normalised power exceeds an Exp(1) threshold about
4e3 times more often than that distribution predicts, and answers: "That would matter if
the searches had assumed Exp(1); they do not -- each estimates the scale empirically as
mu = median(R)/ln 2."

THAT ADDRESSES THE SCALE AND NOT THE TAIL SHAPE, and they are different failures. If R
were exponential with a mis-estimated mean, an empirical mu fixes it entirely. If R's
tail is HEAVIER than exponential, mu*ln(N/alpha) under-covers even with mu exactly right,
because the formula's derivation assumes the exponential form and not merely its scale.

On MAVEN, with the empirical mu in use, the measured band-wide rate was 197 exceedances
per surrogate against an expected 0.05 -- so the empirical scale did not fix it there.
This asks the same question of GOES XRS-B, which is the channel behind section 4.2's
limits, using the same surrogate the fast search uses.

IF THE RATE COMES BACK NEAR 0.05 the defence holds on GOES and only section 4.12 moves.
IF IT COMES BACK IN THE HUNDREDS the paper's headline 1.4e-6 is quoted against a
threshold that is not family-wise 0.05 either, and every fast-band limit requotes.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import fastsearch as F

NSUR = int(os.environ.get("NSUR", "200"))

if __name__ == "__main__":
    t0 = time.time()
    A16, B16, s16 = F.load("g16")
    r, ok = F.prep(B16)
    f, P = F.spectrum(r)
    C = F.continuum(P); R = P/C
    mu = np.median(R)/np.log(2.0)
    N = len(R)
    thr = mu*np.log(N/0.05)
    print("GOES-16 XRS-B: %d bins, mu = %.4f (section 4.2 quotes 1.44)" % (N, mu))
    print("  analytic threshold mu*ln(N/alpha) = %.2f\n" % thr, flush=True)

    n = len(r)
    rng = np.random.default_rng(777)
    counts, maxes = [], []
    for i in range(NSUR):
        sur = np.roll(r, int(rng.integers(1000, n-1000))).copy()
        sur[~ok] = 0.0
        _, P2 = F.spectrum(sur)
        R2 = P2/F.continuum(P2)
        counts.append(int((R2 > thr).sum()))
        maxes.append(float(R2.max()))
    counts = np.array(counts); maxes = np.array(maxes)
    emp = float(np.quantile(maxes, 0.95))
    print("  BAND-WIDE FALSE ALARMS per surrogate over %d bins" % N)
    print("    measured  %.2f      expected  0.05" % counts.mean())
    print("    ratio     %.0fx" % (counts.mean()/0.05) if counts.mean() > 0 else "    ratio     <1")
    print()
    print("  EMPIRICAL threshold, 95th pct of max-R   %8.2f" % emp)
    print("  analytic threshold                       %8.2f" % thr)
    print("  the analytic value is %.2fx too permissive" % (emp/thr))
    print("  limits scale as sqrt(threshold): %.2fx looser" % np.sqrt(max(emp/thr, 1.0)))
    print()
    if counts.mean() <= 0.5:
        print("  --> section 4.2's defence HOLDS on GOES. Only section 4.12 moves.")
    else:
        print("  --> section 4.2's defence FAILS on GOES TOO.")
        print("      1.4e-6 requotes to %.1e" % (1.4e-6*np.sqrt(emp/thr)))
    json.dump(dict(mu=float(mu), n_bins=int(N), analytic=float(thr), empirical=emp,
                   mean_false_alarms=float(counts.mean()), nsur=NSUR),
              open(os.path.expanduser("~/goes_thr_check.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
