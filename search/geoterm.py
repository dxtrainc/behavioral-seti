"""
ITEM G -- THE ONE-YEAR GEOMETRIC TERM.

A modulator fixed in inertial space, viewed from a moving Earth, does not present
a stationary carrier. Earth's distance to a point near the Sun varies by the
orbital eccentricity, +/-0.0167 AU = +/-8.3 light-seconds, and its projected
velocity varies likewise; both impose an ANNUAL modulation on the carrier's
phase and amplitude. Every fast-band search in section 4.2 assumes a stationary
carrier and integrates coherently for 1,945 days. If an annual term is present,
that integration is wrong and the search loses power to sidebands.

WHERE THE SIDEBANDS ARE. Over T = 1,945 d the frequency resolution is
df = 5.95e-9 Hz and f_yr = 3.169e-8 Hz, so the sidebands sit 5.33 bins from the
carrier -- resolvable, but NOT on an integer bin, so each is split across two
bins and further spread by the Blackman window's three-bin mainlobe. The
statistic must therefore sum a block, not sample a bin:

    R_geo(i) = R(i) + R(i-6) + R(i-5) + R(i+5) + R(i+6)

THRESHOLD. The searches estimate the scale empirically as mu = median(R)/ln 2
rather than assuming Exp(1). A five-bin sum is Gamma(5, mu) IF the bins are
independent, which the window makes false. So the threshold is set from the
Gamma form and then the false-alarm rate is MEASURED by injection at depth zero,
which is what decides whether the threshold is honest.

WHAT IS BEING ASKED. Two questions, and they are different:
  1 does the geometric statistic recover an annually-modulated signal that the
    coherent statistic misses?  (injection -- a statement about method)
  2 does any real candidate that was sub-threshold under the coherent model rise
    above threshold under the geometric one?  (the data -- a statement about
    the Sun)
"""
import os, sys, json
import numpy as np
from scipy.ndimage import median_filter
from scipy.special import gammaincinv

sys.path.insert(0, os.path.expanduser("~/beacon-repo/search"))
D = os.path.expanduser("~")
CH = ["irr_256","irr_284","irr_304","irr_1175","irr_1216","irr_1335","irr_1405","MgII_EXIS"]
NICE = {"irr_256":"25.6nm","irr_284":"28.4nm","irr_304":"30.4nm","irr_1175":"117.5nm",
        "irr_1216":"Ly-alpha","irr_1335":"133.5nm","irr_1405":"140.5nm","MgII_EXIS":"MgII"}
STEP = 60.0
YEAR = 365.25*86400.0
rng = np.random.default_rng(7)

# ---- pipeline, identical in arithmetic to euvsfast2.py -------------------
def daily_profile(r, ok, day=1440):
    n = (len(r)//day)*day
    v = r[:n].reshape(-1, day); m = ok[:n].reshape(-1, day)
    prof = np.where(m.sum(axis=0) > 10, (v*m).sum(axis=0)/np.maximum(m.sum(axis=0), 1), 0.0)
    out = r.copy(); out[:n] = (v-prof).ravel(); out[~ok] = 0.0
    return out

def prep(x, wmin=1440):
    ok = np.isfinite(x) & (x > 0)
    y = np.full(len(x), np.nan); y[ok] = np.log(x[ok])
    med = np.nanmedian(y); s = 1.4826*np.nanmedian(np.abs(y-med))
    y = np.clip(y, med-5*s, med+5*s)
    a = np.nan_to_num(y); m = ok.astype(np.float64)
    ca = np.concatenate([[0], np.cumsum(a)]); cm = np.concatenate([[0], np.cumsum(m)])
    h = wmin//2; i = np.arange(len(y))
    lo = np.maximum(0, i-h); hi = np.minimum(len(y), i+h)
    num = ca[hi]-ca[lo]; den = cm[hi]-cm[lo]
    trend = np.where(den > 10, num/np.maximum(den, 1), med)
    return np.where(ok, y-trend, 0.0), ok

def notch(r, ok, periods):
    t = np.arange(len(r), dtype=np.float64)*STEP
    A = np.vstack(sum(([np.cos(2*np.pi*t/P), np.sin(2*np.pi*t/P)] for P in periods), [])).T
    coef, *_ = np.linalg.lstsq(A[ok], r[ok], rcond=None)
    out = r - A@coef; out[~ok] = 0.0
    return out

def spec(r):
    n = len(r)//2*2; w = np.blackman(n)
    P = np.abs(np.fft.rfft(r[:n]*w))**2
    return np.fft.rfftfreq(n, d=STEP)[1:], P[1:]

def cont(P, w=801):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

LINES = [86400.0/k for k in (1, 2, 3, 4, 5, 6)]

# ---------------------------------------------------------------------------
# THE NULL. Summing five correlated bins of a heavy-tailed R does NOT give
# Gamma(5): measured against the Gamma(5) threshold the statistic fires 3-5x
# more often than the single-bin statistic does against its own. So the Gamma
# threshold cannot be used to declare a candidate, and an earlier version of
# this script that did so produced 63,856 "new candidates" in one channel, all
# of them threshold inflation.
#
# The matched null is a DECOY OFFSET. The annual sideband sits 5.33 bins from
# the carrier; a sideband at any other spacing has identical statistical
# structure -- same number of bins, same correlations, same heavy tail -- and no
# geometric meaning. Counting exceedances at the true offset against the
# distribution of counts over decoy offsets is a control that CAN fail, and it
# needs no distributional assumption at all.
DECOYS = [(2, 3), (3, 4), (8, 9), (9, 10), (12, 13), (15, 16), (19, 20), (24, 25)]

def geo(R, off=(5, 6)):
    """carrier plus both annual sidebands, each split across two bins"""
    g = R.copy()
    for d in off:
        g[d:] += R[:-d]          # lower sideband contribution
        g[:-d] += R[d:]          # upper
    return g

def thresholds(R, k, alpha=0.05):
    """R is always the SINGLE-BIN ratio: the scale mu belongs to one bin, and
    a k-bin sum is Gamma(k, mu) with the SAME mu. Estimating mu from the summed
    series instead inflates the threshold by ~6.5x and silently disables the
    test -- which is what the first version of this script did."""
    mu = np.median(R)/np.log(2.0)
    N = len(R)
    return mu*gammaincinv(k, 1.0 - alpha/N)

def run_channel(x, tag, inject=None):
    r, ok = prep(x)
    if ok.mean() < 0.3: return None
    if inject is not None:
        f0, amp, annual = inject
        t = np.arange(len(r), dtype=np.float64)*STEP
        env = 1.0 + annual*np.cos(2*np.pi*t/YEAR)
        r = r + amp*env*np.cos(2*np.pi*f0*t)
    r = daily_profile(r, ok)
    r = notch(r, ok, LINES)
    f, P = spec(r); C = cont(P); R = P/C
    Rg = geo(R)
    thr1 = thresholds(R, 1); thrG = thresholds(R, 5)
    # The five summed bins are NOT independent -- the Blackman mainlobe is three
    # bins wide -- so Gamma(5) is an approximation. An earlier version tried to
    # replace it with the empirical (1 - alpha/N) quantile of Rg; that is
    # degenerate, because with N = 1.4e6 samples the (1 - 3.6e-8) quantile is
    # just the maximum, which forces exactly one exceedance by construction.
    # The Gamma form is kept and the excess over it is REPORTED rather than
    # absorbed, exactly as section 4.5 does for the single-bin case.
    q = float(np.quantile(Rg, 1.0 - 1e-4))
    return dict(f=f, R=R, Rg=Rg, thr1=thr1, thrG=thrG,
                thrG_gamma=thrG, thrG_emp=q, ok=float(ok.mean()))

def main():
    z = np.load(os.path.join(D, "euvs1m_g16.npz"))
    out = {}

    print("=== 1. THE DATA: true annual offset against decoy offsets ===\n", flush=True)
    print("  the count at the true offset is compared with the counts at eight", flush=True)
    print("  decoy offsets of identical structure and no geometric meaning\n", flush=True)
    print("  %-10s %6s %7s %7s %8s %8s %7s %6s" %
          ("channel", "%pres", "n>thr1", "n_true", "decoy med", "decoy max", "p", "verdict"),
          flush=True)
    for c in CH:
        res = run_channel(z[c], c)
        if res is None:
            print("  %-10s skipped (coverage)" % NICE[c], flush=True); continue
        R, thr1 = res["R"], res["thr1"]
        thrG = thresholds(R, 5)
        n1 = int((R > thr1).sum())
        nt = int((geo(R) > thrG).sum())
        nd = np.array([int((geo(R, off=d) > thrG).sum()) for d in DECOYS])
        p = (1 + (nd >= nt).sum())/(1 + len(nd))
        print("  %-10s %5.1f%% %7d %7d %8.0f %8d %7.3f  %s" %
              (NICE[c], 100*res["ok"], n1, nt, np.median(nd), nd.max(), p,
               "excess" if p < 0.05 else "none"), flush=True)
        out[NICE[c]] = {"n1": n1, "n_true": nt, "decoy_med": float(np.median(nd)),
                        "decoy_max": int(nd.max()), "decoys": nd.tolist(),
                        "p": float(p), "thrG": float(thrG), "thr1": float(thr1)}

    print("\n=== 2. THE METHOD: can the geometric statistic see what the", flush=True)
    print("       coherent one misses? ===", flush=True)
    print("  Injected at P = 617.3 s (NOT a daily harmonic -- 600 s is the 144th,", flush=True)
    print("  and the daily-profile subtraction annihilated it in the first run).", flush=True)
    print("  Annual envelope depth 1.0 = full amplitude swing over the year.\n", flush=True)
    x = z["irr_1216"]; f0 = 1/617.3
    base = run_channel(x, "cal")
    df = base["f"][1]-base["f"][0]
    j = int(round(f0/df))
    print("  %-20s %10s %8s %10s %8s   %s" %
          ("injection", "R(f0)", "thr1", "Rgeo(f0)", "thrG", "coherent / geometric"), flush=True)
    rec = []
    for amp in (0.0, 1e-5, 2e-5, 5e-5, 1e-4):
        for annual in (0.0, 1.0):
            res = run_channel(x, "inj", inject=(f0, amp, annual))
            Rr = res["R"]; thrG = thresholds(Rr, 5); Rg = geo(Rr)
            lo, hi = max(0, j-8), j+9
            k = int(np.argmax(Rr[lo:hi])) + lo
            kg = int(np.argmax(Rg[lo:hi])) + lo
            d1 = bool(Rr[k] > res["thr1"]); dg = bool(Rg[kg] > thrG)
            print("  a=%.0e annual=%.0f   %10.1f %8.1f %10.1f %8.1f   %s / %s" %
                  (amp, annual, Rr[k], res["thr1"], Rg[kg], thrG,
                   "HIT " if d1 else "miss", "HIT " if dg else "miss"), flush=True)
            rec.append({"amp": amp, "annual": annual, "R": float(Rr[k]),
                        "thr1": float(res["thr1"]), "Rgeo": float(Rg[kg]),
                        "thrG": float(thrG), "coh_hit": d1, "geo_hit": dg})
    out["injection"] = rec
    json.dump(out, open(os.path.join(D, "geoterm.json"), "w"), indent=1)
    print("\nsaved ~/geoterm.json", flush=True)

main()
