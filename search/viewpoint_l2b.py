"""
ITEM A -- THE LINE-OF-SIGHT TEST, slow band, on measured data.

A stationary occulter or redistributor modulates flux along ONE line of sight.
Observers on different lines see either nothing or something different. A
modulation intrinsic to the Sun is seen by everyone, offset only by geometry.
So: Earth line (GOES-16 EUVS Lyman-alpha, 1 min) against Mars line (MAVEN EUVM
diode C, Lyman-alpha, orbit-averaged), 2019-12-10 to 2025-04-06.

WHY THIS PRODUCT. L3B is FISM-M model output partly driven by Earth data and
would contaminate the comparison with the line it is being compared against.
L2B is measurement. Its cadence is 3.66 h, so this reaches periods above ~7.3 h
only; the two-minute band waits on the L2 daily files.

THE CONTROL COMES FIRST, AND IT IS NOT OPTIONAL.
Solar rotation already provides a signal both platforms must see, with a lag
fixed by geometry: a feature facing Earth faces Mars after the Sun has turned
through the heliocentric longitude difference, tau = (lam_M - lam_E)/360 x
27.2753 d. Over this span the Earth-Mars angle sweeps through more than two
synodic periods, so tau ranges over the whole rotation. If the pipeline cannot
recover THAT -- a known signal at a predicted, time-varying lag -- then a null
on anything else means nothing at all. Measured lag against predicted lag is
therefore the first result, and it can fail.

TWO NULLS, per section 3.4 and Fable's note:
  anti-geometric  apply the light-travel and lag corrections with the WRONG
                  SIGN; coherence must vanish
  shuffled        break the pairing between window and its predicted lag
"""
import os, json, datetime as dt
import numpy as np
from scipy.io import readsav
from scipy.ndimage import median_filter

AU_KM = 1.495978707e8
C_KMS = 299792.458
CARR = 27.2753
D = os.path.expanduser("~")
rng = np.random.default_rng(3)

# ---------------------------------------------------------------- ephemeris
def longitudes(jds):
    """heliocentric ecliptic longitude of Earth and Mars, degrees"""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
    solar_system_ephemeris.set("builtin")
    T = Time(jds, format="jd")
    sun = get_body_barycentric("sun", T)
    out = []
    for body in ("earth", "mars"):
        p = get_body_barycentric(body, T) - sun
        x = p.x.to_value("AU"); y = p.y.to_value("AU"); z = p.z.to_value("AU")
        out.append((np.degrees(np.arctan2(y, x)) % 360.0, np.sqrt(x*x+y*y+z*z)))
    (le, re), (lm, rm) = out
    return le, re, lm, rm

# ---------------------------------------------------------------- series
def maven_daily():
    d = readsav(os.path.join(D, "mvn_euv_l2b.sav"), verbose=False)["mvn_euv_l2_orbit"][0]
    t = np.asarray(d["TIME_UNIX"], float)
    c = np.asarray(d["DIODE_MEAN"], float)[2]          # diode C = Ly-alpha
    r = np.asarray(d["MEAN_SUN_DISTANCE"], float)/AU_KM
    ok = np.isfinite(t) & (t > 0) & np.isfinite(c) & (c > 0) & np.isfinite(r) & (r > 0.5)
    t, c, r = t[ok], c[ok], r[ok]
    c = c*r*r                                          # normalise to 1 AU
    t = t - r*AU_KM/C_KMS                              # to photon-emission time
    jd = t/86400.0 + 2440587.5
    return jd, c

def goes_daily():
    z = np.load(os.path.join(D, "euvs1m_g16.npz"))
    x = np.asarray(z["irr_1216"], float)
    start = int(z["start"])
    jd0 = start + 1721424.5                            # ordinal -> JD at 00:00
    jd = jd0 + np.arange(len(x))/1440.0
    ok = np.isfinite(x) & (x > 0)
    return jd[ok], x[ok]

def to_daily(jd, v):
    day = np.floor(jd).astype(np.int64)
    u, inv = np.unique(day, return_inverse=True)
    s = np.bincount(inv, weights=v); n = np.bincount(inv)
    return u.astype(float), s/np.maximum(n, 1), n

def band(y, lo=20.0, hi=35.0):
    """keep the rotational band: difference of two median filters"""
    a = median_filter(y, size=int(lo), mode="nearest")
    b = median_filter(y, size=int(hi*2+1), mode="nearest")
    return a - b

def xcorr_lag(a, b, lags):
    a = a - a.mean(); b = b - b.mean()
    sa, sb = a.std(), b.std()
    if sa <= 0 or sb <= 0: return np.zeros(len(lags))
    out = []
    n = len(a)
    for L in lags:
        if L >= 0: u, v = a[:n-L], b[L:]
        else:      u, v = a[-L:], b[:n+L]
        out.append(float(np.dot(u-u.mean(), v-v.mean())/(len(u)*u.std()*v.std()))
                   if len(u) > 30 and u.std() > 0 and v.std() > 0 else 0.0)
    return np.array(out)

def circdiff(a, b, per=CARR):
    d = (a - b) % per
    return np.where(d > per/2, d - per, d)

def main():
    jm, vm = maven_daily()
    jg, vg = goes_daily()
    dm, ym, nm = to_daily(jm, vm)
    dg, yg, ng = to_daily(jg, vg)
    lo = max(dm.min(), dg.min()); hi = min(dm.max(), dg.max())
    grid = np.arange(lo, hi+1)
    print("overlap: JD %.0f to %.0f  (%d days, %.2f yr)"
          % (lo, hi, len(grid), len(grid)/365.25), flush=True)

    E = np.interp(grid, dg, yg)
    M = np.interp(grid, dm, ym)
    # gaps: MAVEN has real ones; mark days with no orbit within 1 day
    okM = np.array([np.min(np.abs(dm-g)) <= 1.0 for g in grid])
    print("MAVEN days with data: %d/%d (%.1f%%)" % (okM.sum(), len(grid), 100*okM.mean()), flush=True)

    le, re, lm, rm = longitudes(grid)
    dlam = (lm - le) % 360.0
    tau_pred = dlam/360.0*CARR
    tau_pred = np.where(tau_pred > CARR/2, tau_pred-CARR, tau_pred)
    print("Earth-Mars longitude difference sweeps %.0f to %.0f deg"
          % (dlam.min(), dlam.max()), flush=True)

    Eb, Mb = band(E), band(M)
    lags = np.arange(-13, 14)

    print("\n=== POSITIVE CONTROL: does the measured lag track the geometry? ===",
          flush=True)
    W, S = 121, 20
    meas, pred, amp, mid = [], [], [], []
    for s in range(0, len(grid)-W, S):
        sl = slice(s, s+W)
        if okM[sl].mean() < 0.7: continue
        cc = xcorr_lag(Eb[sl], Mb[sl], lags)
        k = int(np.argmax(cc))
        meas.append(float(lags[k])); pred.append(float(np.median(tau_pred[sl])))
        amp.append(float(cc[k])); mid.append(float(grid[s+W//2]))
    meas = np.array(meas); pred = np.array(pred); amp = np.array(amp)
    print("  windows used: %d   median peak correlation %.3f" % (len(meas), np.median(amp)),
          flush=True)
    resid = circdiff(meas, pred)
    rms = float(np.sqrt(np.mean(resid**2)))
    # null: break the window-to-lag pairing
    nulls = []
    for _ in range(2000):
        p = rng.permutation(len(pred))
        nulls.append(float(np.sqrt(np.mean(circdiff(meas, pred[p])**2))))
    nulls = np.array(nulls)
    p = (1 + (nulls <= rms).sum())/(1 + len(nulls))
    print("  RMS(measured - predicted lag) = %.2f d" % rms, flush=True)
    print("  shuffled-pairing null: median %.2f d, 5th pct %.2f d" %
          (np.median(nulls), np.percentile(nulls, 5)), flush=True)
    print("  p = %.4f  -> %s" % (p, "GEOMETRY RECOVERED" if p < 0.05
                                 else "geometry NOT recovered"), flush=True)

    # anti-geometric control: wrong-sign prediction must NOT match
    rms_anti = float(np.sqrt(np.mean(circdiff(meas, -pred)**2)))
    print("  anti-geometric (wrong sign): RMS %.2f d  %s" %
          (rms_anti, "(correctly worse)" if rms_anti > rms else "(!! not worse)"), flush=True)

    json.dump({"n_win": len(meas), "rms": rms, "p": float(p),
               "rms_anti": rms_anti, "null_med": float(np.median(nulls)),
               "peak_corr_med": float(np.median(amp)),
               "meas": meas.tolist(), "pred": pred.tolist(), "mid": mid,
               "overlap_days": int(len(grid)),
               "maven_coverage": float(okM.mean())},
              open(os.path.join(D, "viewpoint_l2b.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_l2b.json", flush=True)

if __name__=="__main__": main()
