"""
LUNAR LASER RANGING -- a MODEL-FREE search for non-random periodic structure.

The point of doing it this way. LLR residuals as usually published are what is
left after fitting a lunar ephemeris with hundreds of parameters, and anything
periodic in there has had every chance to be absorbed by the fit. But a laser
round trip is self-describing: over one observing session the range is a smooth
function of time, so a low-order polynomial removes it with NO ephemeris, no
station coordinates, no libration model, no tidal model. What is left is
measurement residual, and periodic structure in that needs no dynamics at all.

Three tests, none of which requires a physical model:

 1 FIRE-EPOCH QUANTISATION. The laser fires on a clock. Returns can only occur
   on that grid, so the fractional part of the epoch is not uniform. This is a
   real non-random periodic artifact and it is EXPECTED; the test is whether the
   grid is the one the hardware implies, and whether anything else rides on it.

 2 RESIDUAL PERIODICITY. Lomb-Scargle on the pooled within-session residuals,
   since sampling is irregular. Null: shuffle residuals WITHIN each session,
   which preserves the sampling pattern and the residual distribution exactly
   and destroys only the time-ordering. That is the surrogate matched to this
   statistic in the sense of section 3.4.

 3 SLOT OCCUPANCY. Given the fire grid, which slots actually return photons?
   A detector or timing artifact shows as non-uniform occupancy of the grid.
"""
import os, glob, sys
import numpy as np

C = 2.99792458e8
D = os.path.expanduser("~/llrfr")

def sessions():
    """(epoch seconds-of-day, two-way ToF) per observing session"""
    out = []
    for path in sorted(glob.glob(os.path.join(D, "*.frd"))):
        cur_t, cur_v = [], []
        for line in open(path, errors="ignore"):
            f = line.split()
            if not f: continue
            if f[0] == "h4":
                if len(cur_t) >= 25: out.append((np.array(cur_t), np.array(cur_v)))
                cur_t, cur_v = [], []
            elif f[0] == "10" and len(f) >= 3:
                try:
                    t = float(f[1]); v = float(f[2])
                except ValueError:
                    continue
                if v > 2.0 and v < 3.0: cur_t.append(t); cur_v.append(v)
        if len(cur_t) >= 25: out.append((np.array(cur_t), np.array(cur_v)))
    return out

def unwrap_day(t):
    """seconds-of-day resets at midnight; a session that crosses it would
    otherwise be fitted as if time ran backwards"""
    t = t.copy()
    j = np.where(np.diff(t) < -43200)[0]
    for k in j: t[k+1:] += 86400.0
    return t

def detrend(t, v, half=12, deg=2):
    """LOCAL quadratic over a sliding window of neighbours. A single global
    polynomial was wrong: over a long session the range change is far from
    polynomial and the first version left 876 km of 'residual', which is the
    lunar orbit, not noise. Local fitting makes no model assumption at all
    beyond smoothness over a few tens of seconds."""
    o = np.argsort(t); t, v = t[o], v[o]
    n = len(t); out = np.full(n, np.nan)
    for i in range(n):
        lo, hi = max(0, i-half), min(n, i+half+1)
        if hi-lo < deg+2: continue
        tt = t[lo:hi]-t[i]
        try:
            c = np.polyfit(tt, v[lo:hi], deg)
        except Exception:
            continue
        out[i] = v[i] - np.polyval(c, 0.0)
    return out

def main():
    S = sessions()
    n = sum(len(t) for t, _ in S)
    print("sessions %d   returns %d   median per session %d"
          % (len(S), n, int(np.median([len(t) for t, _ in S]))), flush=True)

    # ---------------- 1. fire-epoch quantisation ----------------
    print("\n=== 1. fire-epoch grid ===", flush=True)
    allt = np.concatenate([t for t, _ in S])
    # A coarse scan over "nice" quanta is ambiguous: if epochs sit on a 0.05 s
    # grid then 0.1, 0.01 and 0.001 all score too, because each is commensurate.
    # The unambiguous question is which FREQUENCY the epoch series is combed at,
    # so scan the comb frequency continuously and read off the fundamental.
    fs = np.linspace(1.0, 2000.0, 200000)         # Hz -- the first scan
                                                  # stopped at 400 Hz and the
                                                  # peak sat on that ceiling
    z = np.array([np.abs(np.mean(np.exp(2j*np.pi*f*allt))) for f in fs]) \
        if False else None
    # vectorised: chunk to keep memory sane
    z = np.empty(len(fs))
    CH = 2000
    for i in range(0, len(fs), CH):
        ph = np.exp(2j*np.pi*np.outer(fs[i:i+CH], allt))
        z[i:i+CH] = np.abs(ph.mean(axis=1))
    order = np.argsort(-z)[:400]
    peaks = []
    for k in order:
        if all(abs(fs[k]-p0) > 0.5 for p0, _ in peaks):
            peaks.append((float(fs[k]), float(z[k])))
        if len(peaks) >= 8: break
    peaks.sort()
    print("   comb frequencies in the fire epochs (concentration 0 = uniform):", flush=True)
    for f0, zz in peaks:
        print("     %8.3f Hz  = %7.4f s   concentration %.4f" % (f0, 1.0/f0, zz), flush=True)
    fund = max(peaks, key=lambda pz: pz[1])[0]   # strongest, not lowest
    print("   -> strongest comb %.3f Hz = %.6f s" % (fund, 1.0/fund), flush=True)
    # a uniform series of the same length gives this much concentration by chance
    print("   chance level for %d epochs: %.4f" % (len(allt), 1.0/np.sqrt(len(allt))), flush=True)
    best = (1.0/fund, max(z))

    # ---------------- 2. residual periodicity ----------------
    print("\n=== 2. residual periodicity, model-free ===", flush=True)
    from scipy.signal import lombscargle
    R, T = [], []
    for t, v in S:
        t = unwrap_day(t)
        if np.ptp(t) < 60: continue
        r = detrend(t, v)
        ok = np.isfinite(r)
        if ok.sum() < 25: continue
        t, r = t[ok], r[ok]
        s = r.std()
        if not np.isfinite(s) or s <= 0: continue
        keep = np.abs(r) < 5*s              # drop obvious outliers
        R.append(r[keep]/s); T.append(t[keep]-t[keep].min())
    print("   usable sessions %d, residuals %d" % (len(R), sum(len(x) for x in R)), flush=True)
    rms_ps = np.median([np.std(x) for x in R])
    allr = np.concatenate(R)
    pooled = []
    for t, v in S:
        t = unwrap_day(t)
        if np.ptp(t) < 60: continue
        r = detrend(t, v); r = r[np.isfinite(r)]
        if len(r): pooled.append(r)
    pooled = np.concatenate(pooled) if pooled else np.array([0.0])
    rms_s = float(np.std(pooled))
    print("   pooled residual: %.1f ps  =  %.1f mm one-way"
          % (rms_s*1e12, rms_s*C/2*1e3), flush=True)
    print("   (LLR single-shot residuals are tens to hundreds of ps; if this is"
          " far outside that, the detrend is wrong, not the Moon)", flush=True)

    periods = np.logspace(np.log10(0.5), np.log10(600.0), 400)
    w = 2*np.pi/periods
    def power(Rs, Ts):
        acc = np.zeros(len(w))
        for r, t in zip(Rs, Ts):
            if len(r) < 25 or np.ptp(t) < 2*periods[0]: continue
            acc += lombscargle(t.astype(float), r.astype(float), w, normalize=True)
        return acc/max(len(Rs), 1)
    obs = power(R, T)
    rng = np.random.default_rng(5)
    NULL = 200
    nulls = np.empty((NULL, len(w)))
    for k in range(NULL):
        Rs = [rng.permutation(r) for r in R]     # shuffle WITHIN session
        nulls[k] = power(Rs, T)
    hi = np.percentile(nulls, 99.9, axis=0)
    exc = obs > hi
    print("   frequencies above the 99.9th percentile of the shuffled null: %d of %d"
          % (exc.sum(), len(w)), flush=True)
    gmax = float((obs/np.maximum(nulls.max(axis=0), 1e-12)).max())
    k = int(np.argmax(obs/np.maximum(hi, 1e-12)))
    pglob = (1 + (nulls.max(axis=1) >= obs.max()).sum())/(1+NULL)
    print("   strongest: P = %.2f s, power %.4g, null 99.9%% %.4g" % (periods[k], obs[k], hi[k]), flush=True)
    print("   global p (max over band vs null max): %.4f  -> %s"
          % (pglob, "STRUCTURE" if pglob < 0.05 else "no structure"), flush=True)

    # ---------------- 3. slot occupancy ----------------
    print("\n=== 3. occupancy of the fire grid ===", flush=True)
    q = best[0]
    slot = np.round((allt/q) % 1.0 * 20).astype(int) % 20
    cnt = np.bincount(slot, minlength=20).astype(float)
    chi = ((cnt-cnt.mean())**2/cnt.mean()).sum()
    print("   20 phase bins, chi2 = %.1f on 19 dof" % chi, flush=True)
    print("   occupancy: %s" % " ".join("%3.0f" % c for c in cnt), flush=True)

if __name__ == "__main__":
    main()
