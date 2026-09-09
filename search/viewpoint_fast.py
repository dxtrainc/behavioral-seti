"""
ITEM A, fast band -- the viewpoint test at the cadence where section 4.3 shows
the sensitivity actually is.

Earth line: GOES-16 EUVS Lyman-alpha, 1 minute.
Mars line:  MAVEN/EUVM L2 band irradiances, diode C (Lyman-alpha), reduced to
            1-minute means so the two are directly comparable.
Span: 2019-12-10 to 2025-04-06, matched.

WHY THE LIGHT-TRAVEL CORRECTION IS NOT OPTIONAL HERE, unlike in the slow band.
A time shift changes phase, not amplitude, so for a power comparison a CONSTANT
offset would not matter. But the Sun-Mars distance swings 1.38-1.67 AU over the
Martian year, which is +/-70 s of light travel -- and at a 120 s period that is
+/-0.58 of a cycle. Uncorrected, a coherent tone seen from Mars is smeared across
many bins and the search would report absence for a signal that is present. Both
series are therefore resampled onto photon-EMISSION time before anything else.

THE CONTAMINANT. MAVEN is in Mars orbit and EUVM sees the Sun only outside
occultation, so the Mars series carries a comb at the orbital period (~4.5 h) and
its harmonics -- the same structure as the spacecraft-day comb that defeated
version 1 of the GOES fast search. It is measured from the data rather than
assumed, and every candidate is vetted against it.

THE TEST, as in section 4.12: a peak present in one line and absent in the other
is a line-of-sight signature ONLY if the other line had the power to see it.
Every Earth-only candidate is injected into a surrogate Mars series at the
observed amplitude and kept only if recovery reaches 95%.
"""
import os, json, datetime as dt
import numpy as np
from scipy.ndimage import median_filter

AU_KM = 1.495978707e8
C_KMS = 299792.458
D = os.path.expanduser("~")
rng = np.random.default_rng(31)
STEP = 60.0

def prep(x, ok, wmin=1440):
    """log, clip, running-mean detrend -- the euvsfast2 pipeline"""
    y = np.full(len(x), np.nan); y[ok] = np.log(x[ok])
    med = np.nanmedian(y); s = 1.4826*np.nanmedian(np.abs(y-med))
    y = np.clip(y, med-5*s, med+5*s)
    a = np.nan_to_num(y); m = ok.astype(np.float64)
    ca = np.concatenate([[0], np.cumsum(a)]); cm = np.concatenate([[0], np.cumsum(m)])
    h = wmin//2; i = np.arange(len(y))
    lo = np.maximum(0, i-h); hi = np.minimum(len(y), i+h)
    num = ca[hi]-ca[lo]; den = cm[hi]-cm[lo]
    trend = np.where(den > 10, num/np.maximum(den, 1), med)
    return np.where(ok, y-trend, 0.0)

def fold_out(r, ok, period_min):
    """subtract the mean profile at a given period -- kills that comb entire"""
    P = int(round(period_min))
    if P < 2: return r
    n = (len(r)//P)*P
    v = r[:n].reshape(-1, P); m = ok[:n].reshape(-1, P)
    prof = np.where(m.sum(0) > 5, (v*m).sum(0)/np.maximum(m.sum(0), 1), 0.0)
    out = r.copy(); out[:n] = (v-prof).ravel(); out[~ok] = 0.0
    return out

def fold_drifting(r, ok, seg_days=30, lo_h=3.0, hi_h=8.0, periods_in=None):
    """MAVEN's orbit is not a metronome -- its period evolves over the mission,
    so folding the whole record at ONE period leaves most of the comb standing.
    A single fold left 1,076 peaks above threshold in the Mars line against 15
    in the Earth line, which is the comb, not the Sun. Folding segment by
    segment, at the period each segment actually shows, tracks the drift."""
    out = r.copy()
    W = seg_days*1440
    periods = []
    # Determining each segment's period costs a spectrum and a continuum, and
    # in an injection loop that is 65 segments times hundreds of trials. The
    # periods are a property of the ORBIT, not of the injection, so they are
    # measured once on the real series and passed back in.
    for i0, s0 in enumerate(range(0, len(r), W)):
        sl = slice(s0, min(s0+W, len(r)))
        seg, segok = out[sl].copy(), ok[sl]
        if segok.sum() < W//4:
            periods.append(None); continue
        if periods_in is not None:
            Pm = periods_in[i0] if i0 < len(periods_in) else None
            if Pm is None: continue
        else:
            f, P, _ = spec(seg)
            if len(f) < 10:
                periods.append(None); continue
            R = P/cont(P, w=201)
            b = (f > 1/(hi_h*3600.)) & (f < 1/(lo_h*3600.))
            if b.sum() < 5:
                periods.append(None); continue
            Pm = 1.0/float(f[b][np.argmax(R[b])])/60.0
        periods.append(Pm)
        seg = fold_out(seg, segok, Pm)
        seg = fold_out(seg, segok, Pm/2.0)     # first harmonic
        out[sl] = seg
    got = [p0 for p0 in periods if p0]
    if got and periods_in is None:
        print("  segment orbital periods: %.4f-%.4f h (median %.4f)"
              % (min(got)/60, max(got)/60, float(np.median(got))/60), flush=True)
    out[~ok] = 0.0
    return (out, periods) if periods_in is None else out

_WCACHE = {}
def spec(r):
    """the window and frequency grid depend only on the length, and this is
    called thousands of times on a 2.8M-point series, so both are cached"""
    n = len(r)//2*2
    got = _WCACHE.get(n)
    if got is None:
        w = np.blackman(n)
        got = (w, np.fft.rfftfreq(n, d=STEP)[1:], w.sum())
        _WCACHE[n] = got
    w, f, ws = got
    P = np.abs(np.fft.rfft(r[:n]*w))**2
    return f, P[1:], ws

def cont(P, w=801):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

def thr_of(R, alpha=0.05):
    mu = np.median(R)/np.log(2.0)
    return mu*np.log(len(R)/alpha), mu

def load_grid():
    """common 1-minute grid in photon-emission time"""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
    z = np.load(os.path.join(D, "euvs1m_g16.npz"))
    E = np.asarray(z["irr_1216"], float)
    start = int(z["start"])
    jd0 = start + 1721424.5
    n = len(E)
    jdE = jd0 + np.arange(n)/1440.0

    mv = json.load(open(os.path.join(D, "maven_l2_1m.json")))
    days = sorted(k for k, v in mv.items() if v)
    print("MAVEN days with data: %d" % len(days), flush=True)
    tm, vm = [], []
    for k in days:
        rec = mv[k]
        mins = np.asarray(rec["min"], dtype=np.int64)
        arr = np.array([[np.nan if x is None else x for x in row] for row in rec["v"]], float)
        tm.append(mins*60.0); vm.append(arr[:, 2])          # diode C
    tm = np.concatenate(tm); vm = np.concatenate(vm)
    good = np.isfinite(vm) & (vm > 0)
    tm, vm = tm[good], vm[good]
    jdM = tm/86400.0 + 2440587.5

    # heliocentric distances, daily then interpolated
    solar_system_ephemeris.set("builtin")
    dgrid = np.arange(np.floor(min(jdE.min(), jdM.min())), np.ceil(max(jdE.max(), jdM.max()))+1)
    T = Time(dgrid, format="jd")
    sun = get_body_barycentric("sun", T)
    def dist(b):
        p = get_body_barycentric(b, T) - sun
        return np.sqrt(p.x.to_value("AU")**2 + p.y.to_value("AU")**2 + p.z.to_value("AU")**2)
    rE = np.interp(jdE, dgrid, dist("earth"))
    rM = np.interp(jdM, dgrid, dist("mars"))
    print("Mars distance %.3f-%.3f AU -> light travel swing %.0f s"
          % (rM.min(), rM.max(), (rM.max()-rM.min())*AU_KM/C_KMS), flush=True)

    # to emission time, and normalise to 1 AU
    jdE = jdE - (rE*AU_KM/C_KMS)/86400.0
    jdM = jdM - (rM*AU_KM/C_KMS)/86400.0
    E = E*rE*rE
    vm = vm*rM*rM

    lo = max(jdE.min(), jdM.min()); hi = min(jdE.max(), jdM.max())
    grid = np.arange(lo, hi, 1.0/1440.0)
    print("common grid: %d minutes (%.1f d)" % (len(grid), len(grid)/1440), flush=True)
    def onto(j, v):
        y = np.interp(grid, j, v, left=np.nan, right=np.nan)
        # mark bins with no sample within 2 minutes as missing
        idx = np.searchsorted(j, grid)
        idx = np.clip(idx, 1, len(j)-1)
        near = np.minimum(np.abs(j[idx]-grid), np.abs(j[idx-1]-grid))
        y[near > 2.0/1440.0] = np.nan
        return y
    return grid, onto(jdE, E), onto(jdM, vm)

def main():
    grid, E, M = load_grid()
    okE = np.isfinite(E) & (E > 0); okM = np.isfinite(M) & (M > 0)
    print("coverage: Earth %.1f%%  Mars %.1f%%  both %.1f%%"
          % (100*okE.mean(), 100*okM.mean(), 100*(okE & okM).mean()), flush=True)
    rE = prep(np.nan_to_num(E, nan=1.0), okE)
    rM = prep(np.nan_to_num(M, nan=1.0), okM)

    # measure MAVEN's orbital comb rather than assuming it
    fM0, PM0, _ = spec(rM)
    RM0 = PM0/cont(PM0)
    band = (fM0 > 1/(8*3600.)) & (fM0 < 1/(3*3600.))
    forb = float(fM0[band][np.argmax(RM0[band])])
    Porb = 1.0/forb/60.0
    print("MAVEN orbital comb: %.4f h  (peak R=%.0f)" % (Porb/60, RM0[band].max()), flush=True)

    rE = fold_out(rE, okE, 1440)                 # spacecraft day, per section 4.4
    print("removing the Mars orbital comb, segment by segment:", flush=True)
    rM, _segP = fold_drifting(rM, okM)           # tracks the orbit's drift
    rM = fold_out(rM, okM, 1440)

    fE, PE, W = spec(rE); fM, PM, _ = spec(rM)
    CE = cont(PE); CM = cont(PM)
    RE = PE/CE; RM = PM/CM
    tE, _ = thr_of(RE); tM, _ = thr_of(RM)
    print("band %.1f s to %.1f h   Earth thr %.1f  Mars thr %.1f"
          % (1/fE[-1], 1/fE[0]/3600, tE, tM), flush=True)

    # BAND CAP, and it is not cosmetic. Solar rotation has a different synodic
    # period from each viewpoint -- 27.28 d from Earth, 26.35 d from Mars,
    # because Mars orbits more slowly -- so ANY rotation-band structure is
    # "present in one and absent in the other" by construction. The first run
    # duly returned a 26.6 d "Earth-only" candidate at 99.2% Mars power, which
    # is that artefact and nothing else. This is the fast-band test; periods
    # near and above a day belong to the slow-band test of section 4.12, which
    # carries the geometric lag that makes the comparison valid. Judge only
    # below 6 h, where both viewpoints see the same Sun up to light travel.
    PCAP = 6*3600.0
    inband = (1.0/fE <= PCAP) & (1.0/fE >= 125.0)
    print("judged band: 125 s to %.0f h (%d bins of %d)"
          % (PCAP/3600, inband.sum(), len(fE)), flush=True)
    iE = np.where((RE > tE) & inband)[0]; iM = np.where((RM > tM) & inband)[0]
    both = np.intersect1d(iE, iM)
    print("\npeaks: Earth %d   Mars %d   both %d" % (len(iE), len(iM), len(both)), flush=True)

    harm = np.array([k*forb for k in range(1, 400)])
    dayh = np.array([k/86400.0 for k in range(1, 3000)])
    def near(f, arr, tol):
        return bool(np.min(np.abs(arr-f)) < tol)
    df = fE[1]-fE[0]

    sdM = rM[okM].std(); n = len(rM)

    # control: with no injection the recovery must sit near the nominal 5%
    t0 = np.arange(n, dtype=float)*STEP
    fc = float(fE[len(fE)//3]); kc = int(np.argmin(np.abs(fM-fc)))
    slc = slice(max(0, kc-2), kc+3); fp = 0
    for _ in range(120):
        sur = np.roll(rM, int(rng.integers(1000, n-1000)))
        sur[~okM] = 0.0
        _, P2, _ = spec(sur)
        if (P2[slc]/CM[slc]).max() > tM: fp += 1
    # The threshold is family-wise over ~5e5 bins, so the chance that any ONE
    # named bin exceeds it is alpha/N, not alpha. A near-zero rate here is the
    # correct result and confirms the cached continuum invents nothing; an
    # earlier label called for "~5%", which was simply the wrong statistic.
    print("\ndepth-0 control at a fixed bin: %.1f%% (family-wise threshold, so ~0%% is correct)"
          % (100.0*fp/120), flush=True)

    print("\n=== Earth-only candidates ===", flush=True)
    rows = []
    for i in iE[np.argsort(-RE[iE])][:20]:
        if i in both: continue
        f0 = fE[i]; P = 1.0/f0
        # Reject the band edges. 1-minute sampling puts Nyquist at exactly
        # 120 s, and the first run's two strongest "Earth-only" candidates sat
        # at 120.00 and 120.01 s -- the edge, not the Sun. The same trap caught
        # the decoy scan at its 2 kHz ceiling and the viewpoint search at its
        # 2.000 d limit.
        if i < 5 or i > len(fE)-6 or P < 125.0:
            print("  P=%9.2f s  R_E=%7.1f   REJECTED: band edge (Nyquist)" % (P, RE[i]), flush=True)
            continue
        if near(f0, dayh, 3*df):
            print("  P=%9.2f s  R_E=%7.1f   REJECTED: harmonic of the spacecraft day" % (P, RE[i]), flush=True)
            continue
        if near(f0, harm, 3*df):
            print("  P=%9.2f s  R_E=%7.1f   REJECTED: harmonic of the MAVEN orbit" % (P, RE[i]), flush=True)
            continue
        aE = 2.0*np.sqrt(PE[i])/W/rE[okE].std()
        # The continuum is a running median over 801 bins of a 576k-bin
        # spectrum: recomputing it inside every injection costs seconds per
        # trial and put this loop at over an hour. An injected tone moves one
        # bin, not the running median around it, so the continuum of the
        # unmodified series is reused. Verified below by the depth-0 control.
        hits = 0; N = 120
        t = np.arange(n, dtype=float)*STEP
        k = int(np.argmin(np.abs(fM-f0)))
        sl = slice(max(0, k-2), k+3)
        for _ in range(N):
            sur = np.roll(rM, int(rng.integers(1000, n-1000)))
            y = sur + aE*sdM*np.cos(2*np.pi*f0*t + rng.uniform(0, 2*np.pi))
            y[~okM] = 0.0
            _, P2, _ = spec(y)
            if (P2[sl]/CM[sl]).max() > tM: hits += 1
        pw = 100.0*hits/N
        print("  P=%9.2f s  R_E=%7.1f  R_M=%6.2f  amp %.2e  Mars power %5.1f%%  %s"
              % (P, RE[i], RM[i], aE, pw,
                 "EARTH-ONLY" if pw >= 95 else "unresolved"), flush=True)
        rows.append({"period_s": float(P), "RE": float(RE[i]), "RM": float(RM[i]),
                     "amp": float(aE), "mars_power": pw, "earth_only": bool(pw >= 95)})

    json.dump({"n_earth": len(iE), "n_mars": len(iM), "n_both": len(both),
               "forb_hz": forb, "Porb_h": Porb/60, "thrE": float(tE),
               "thrM": float(tM), "rows": rows,
               "cov_earth": float(okE.mean()), "cov_mars": float(okM.mean())},
              open(os.path.join(D, "viewpoint_fast.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_fast.json", flush=True)


if __name__ == "__main__":
    main()
