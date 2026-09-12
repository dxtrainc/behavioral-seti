"""
THE TWIN-INSTRUMENT TEST.

WHAT IT IS FOR, PRECISELY. It does NOT see through solar nonlinearity -- the Sun's
higher-order structure is present in both spacecraft, so requiring a candidate in
both vetoes instrumental artefacts and nothing else. What it settles is the question
left open by validate/wl_env.py: whether the median z = 45.3 of the quadruple sweep
is the Sun or the pipeline.

GOES-16 and GOES-17 are both geostationary and view the same Sun to within 0.04 s of
light travel. The log difference of the same band on the two spacecraft therefore
cancels the solar signal AND the solar nonlinearity, leaving instrument only. It is a
DARK CHANNEL MADE OF REAL DATA -- real cadence, real gaps, real noise, no Sun --
which is a far better control than synthetic data can be.

  arm A   G16 bands            solar present     expect large z if z measures the Sun
  arm B   G17 bands            solar present     expect the same
  arm C   log G16 - log G17    solar cancelled   z ~ 0 if the pipeline is clean

If arm C returns z ~ 0 while A and B return large z, the void of the quadruple sweep
is solar and the pipeline is exonerated. If arm C also returns large z, the pipeline
manufactures structure and the void is an artefact of it.

MgII_EXIS IS EXCLUDED. It is computed from the bands already in the set, so including
it would repeat the error that produced sixteen spurious survivors in the first pair
sweep: dimensionless is necessary but not sufficient, and a quantity must also be
independently measured.
"""
import os, sys, time, json
import numpy as np
from scipy.stats import norm, kurtosis, kstest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "search"))
import wl_sweep as W          # reuse stat(), FORMS, the form algebra

BANDS = ["irr_256", "irr_284", "irr_304", "irr_1175", "irr_1216", "irr_1335", "irr_1405"]
import itertools
QUADS = list(itertools.combinations(range(len(BANDS)), 4))
TESTS = [(q, fn) for q in QUADS for fn in W.FORMS]

def load(tag):
    z = np.load(os.path.expanduser("~/%s.npz" % tag))
    return float(z["start"]), {b: z[b] for b in BANDS}

def build():
    s16, d16 = load("euvs1m_g16")
    s17, d17 = load("euvs1m_g17")
    o16 = int(round((max(s16, s17) - s16) * 1440))
    o17 = int(round((max(s16, s17) - s17) * 1440))
    n = min(len(d16[BANDS[0]]) - o16, len(d17[BANDS[0]]) - o17)
    A = np.vstack([d16[b][o16:o16+n] for b in BANDS]).astype(float)
    B = np.vstack([d17[b][o17:o17+n] for b in BANDS]).astype(float)
    A[A <= 0] = np.nan; B[B <= 0] = np.nan
    ok = np.all(np.isfinite(A), 0) & np.all(np.isfinite(B), 0)
    return A[:, ok], B[:, ok], ok, n

def ns(x):
    return norm.ppf((np.argsort(np.argsort(x)) + 0.5) / len(x))

def prep(Z, keep):
    """normal-score each channel; supply the gap mask used by the statistic"""
    W.KEEP = keep                      # stat() drops increments spanning a break
    return np.vstack([ns(r) for r in Z])

def stats_of(Z):
    return np.array([W.stat(W.FORMS[fn](*[Z[i] for i in q])) for q, fn in TESTS])

def screen(Z, rng):
    n = Z.shape[1]
    F = np.fft.rfft(Z, axis=1)
    ph = rng.uniform(0, 2*np.pi, F.shape[1]); ph[0] = 0.0
    if n % 2 == 0: ph[-1] = 0.0
    return np.fft.irfft(F * np.exp(1j*ph)[None, :], n=n, axis=1)

def run(Z, keep, nsur, label, seed):
    rng = np.random.default_rng(seed)
    X = prep(Z, keep)
    obs = stats_of(X)
    S = np.array([stats_of(screen(X, rng)) for _ in range(nsur)])
    z = (obs - np.nanmean(S, 0)) / np.nanstd(S, 0)
    zq = np.array([z[t] for t, (q, f) in enumerate(TESTS) if f == "quad-prod"])
    print("  %-28s median z %8.1f   max z %8.1f   |z|>20 %3d/%d   quad-prod med %6.1f"
          % (label, np.nanmedian(z), np.nanmax(z), (np.abs(z) > 20).sum(), len(TESTS),
             np.nanmedian(zq)))
    return z

if __name__ == "__main__":
    NSUR = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    A, B, ok, ntot = build()
    idx = np.where(ok)[0]
    keep = np.diff(idx) == 1
    print("GOES-16 / GOES-17 EUVS, %d bands, 1-min cadence" % len(BANDS))
    print("common finite samples: %d of %d (%.1f%%), %.1f days" % (A.shape[1], ntot,
          100*A.shape[1]/ntot, A.shape[1]/1440))
    print("gap breaks: %d (%.2f%% of steps)\n" % ((~keep).sum(), 100*(~keep).mean()))
    print("tests: %d quadruples x %d forms = %d, %d surrogates\n"
          % (len(QUADS), len(W.FORMS), len(TESTS), NSUR))

    D = np.log(A) - np.log(B)                    # dark channel: solar cancelled
    D = D - np.median(D, axis=1, keepdims=True)

    t0 = time.time()
    zA = run(A, keep, NSUR, "arm A  G16 (solar)", 101)
    zB = run(B, keep, NSUR, "arm B  G17 (solar)", 202)
    zC = run(D, keep, NSUR, "arm C  log G16 - log G17", 303)
    print("\n  (daily six-channel sweep, for comparison: median z 45.3)")
    print("  %.1f s" % (time.time()-t0))

    mA, mC = float(np.nanmedian(zA)), float(np.nanmedian(zC))
    print("\n  VERDICT:")
    if abs(mC) < 5 and abs(mA) > 20:
        print("    arm C near zero, arms A/B large -> pipeline CLEAN, the void is SOLAR")
    elif abs(mC) > 20:
        print("    arm C also large -> pipeline MANUFACTURES structure; void is an artefact")
    else:
        print("    intermediate: %.1f vs %.1f -- neither conclusion is supported" % (mA, mC))
    json.dump({"bands": BANDS, "n": int(A.shape[1]), "nsur": NSUR,
               "z_G16": zA.tolist(), "z_G17": zB.tolist(), "z_diff": zC.tolist()},
              open(os.path.expanduser("~/wl_twin.json"), "w"), indent=1)
