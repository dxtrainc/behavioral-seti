"""
ITEM A, part 2 -- the search itself, once the control of part 1 has earned it.

Part 1 established that the pipeline recovers solar rotation at the predicted,
time-varying geometric lag: RMS 1.97 d against a shuffled-pairing null of 7.65 d,
p = 0.0005, with the wrong-sign control correctly worse. That is outcome 2 in
Fable's scheme -- a real signal seen from both lines with the geometric phase --
and it is what makes a null on anything else meaningful.

THE SEARCH. Align the Mars series onto Earth-equivalent time using the
time-varying lag, then compare spectra. Three outcomes are distinguishable:

  Earth-only peak   present in E above threshold, absent in M WHERE M HAD THE
                    POWER TO SEE IT -> line-of-sight modulator
  both              intrinsic to the Sun
  neither           null

THE TRAP, and it is the same one section 4.10 fell into. MAVEN is orbit-averaged
and noisier than GOES, so "absent at Mars" is worthless unless Mars could have
seen it. Every Earth-only candidate is therefore injected into a surrogate Mars
series at the amplitude observed at Earth, and kept only if recovery exceeds 95%.
An Earth-only peak that Mars had no power to see is reported as unresolved, not
as a detection.
"""
import os, json
import numpy as np
from scipy.ndimage import median_filter
import viewpoint_l2b as V

D = os.path.expanduser("~")
rng = np.random.default_rng(19)

def spec(y):
    n = len(y)//2*2
    w = np.blackman(n)
    P = np.abs(np.fft.rfft((y[:n]-y[:n].mean())*w))**2
    f = np.fft.rfftfreq(n, d=1.0)          # cycles per day
    return f[1:], P[1:]

def cont(P, w=101):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

def thr_of(R, alpha=0.05):
    mu = np.median(R)/np.log(2.0)
    return mu*np.log(len(R)/alpha), mu

def main():
    jm, vm = V.maven_daily(); jg, vg = V.goes_daily()
    dm, ym, _ = V.to_daily(jm, vm); dg, yg, _ = V.to_daily(jg, vg)
    lo = max(dm.min(), dg.min()); hi = min(dm.max(), dg.max())
    grid = np.arange(lo, hi+1)
    E = np.interp(grid, dg, yg); M = np.interp(grid, dm, ym)
    okM = np.array([np.min(np.abs(dm-g)) <= 1.0 for g in grid])

    le, re, lm, rm = V.longitudes(grid)
    # THE CORRECTION THAT WAS MISSING. MAVEN was normalised to 1 AU inside
    # maven_daily(); GOES was not. Earth's orbital eccentricity puts a 6.69%
    # peak-to-peak annual term into any irradiance measured at Earth, purely
    # geometric. Left in, it appears as a large "Earth-only" peak near one year
    # plus a comb of window sidelobes at 1.07, 0.89, 0.76, 0.67 and 0.59 yr --
    # which is exactly what the first run of this search reported, as ten
    # line-of-sight modulator candidates. It was the missing 1/r^2.
    E = E*re*re
    dlam = (lm-le) % 360.0
    tau = dlam/360.0*V.CARR
    tau = np.where(tau > V.CARR/2, tau-V.CARR, tau)

    # Mars onto Earth-equivalent time
    Mal = np.interp(grid+tau, grid, M)
    # okM marks days MAVEN actually observed. They were computed and never used: the
    # 5.7% missing days were linearly interpolated and then treated as data, in the
    # spectrum, the continuum, the threshold and the injections. Interpolated stretches
    # are smoother than data, which lowers Mars's apparent noise and raises Mars power.
    okAl = np.interp(grid+tau, grid, okM.astype(float)) > 0.999
    # A median filter of width W flattens everything slower than about W/2, so
    # the detrend scale sets the longest period the search can honestly claim.
    DETREND = 401
    PMAX = DETREND/2.0
    E = E - median_filter(E, size=DETREND, mode="nearest")
    Mal = Mal - median_filter(Mal, size=DETREND, mode="nearest")

    fE, PE = spec(E); fM, PM = spec(Mal)
    RE = PE/cont(PE); RM = PM/cont(PM)
    keep = (1.0/fE) <= PMAX
    fE, PE, RE = fE[keep], PE[keep], RE[keep]
    fM, PM, RM = fM[keep], PM[keep], RM[keep]
    tE, muE = thr_of(RE); tM, muM = thr_of(RM)
    print("band: %.1f to %.1f day periods, %d bins  (capped at detrend/2)"
          % (1/fE[-1], 1/fE[0], len(fE)), flush=True)
    print("Earth thr %.2f (mu %.3f)   Mars thr %.2f (mu %.3f)" % (tE, muE, tM, muM), flush=True)

    iE = np.where(RE > tE)[0]; iM = np.where(RM > tM)[0]
    both = np.intersect1d(iE, iM)
    print("\npeaks above threshold:  Earth %d   Mars %d   both %d"
          % (len(iE), len(iM), len(both)), flush=True)

    print("\n=== Earth-only candidates, and whether Mars could have seen them ===",
          flush=True)
    Mal = np.where(okAl, Mal, 0.0)
    okM = okAl
    sdM = Mal[okM].std(); n = len(grid)
    rows = []
    for i in iE[np.argsort(-RE[iE])][:25]:
        if i in both: continue
        P = 1.0/fE[i]
        # amplitude at Earth, as a fraction of Earth's own scatter
        # For a coherent tone of amplitude a through window w, the peak power
        # is (a*sum(w)/2)^2, so a = 2*sqrt(P)/sum(w). The first version divided
        # by sqrt(sum(w^2)) instead, which inflated every amplitude by 33.6x and
        # made "could Mars have seen it?" trivially yes for everything.
        WSUM = np.blackman(len(E)//2*2).sum()
        # UNITS. A line-of-sight modulator imposes the same FRACTIONAL amplitude at
        # both viewpoints, so the Earth amplitude must be injected into Mars in the
        # units the series is in -- not rescaled from Earth's scatter to Mars's. The
        # previous form divided by E.std() here and multiplied by sdM below, inflating
        # the injected tone by sdM/sdE and so overstating Mars power, in the direction
        # of declaring Earth-only detections.
        aE = 2.0*np.sqrt(PE[i])/WSUM
        # can Mars see a signal of that fractional amplitude?
        hits = 0; N = 200
        t = np.arange(n, dtype=float)
        for _ in range(N):
            sur = np.roll(Mal, int(rng.integers(50, n-50)))
            sur = sur.copy(); sur[~okM] = 0.0
            y = sur + aE*np.cos(2*np.pi*t/P + rng.uniform(0, 2*np.pi))
            y[~okM] = 0.0
            f2, P2 = spec(y); R2 = P2/cont(P2)
            k = int(np.argmin(np.abs(f2-fE[i])))
            if R2[max(0,k-2):k+3].max() > tM: hits += 1
        pw = 100.0*hits/N
        verdict = ("EARTH-ONLY, Mars had power" if pw >= 95 else
                   "unresolved (Mars power %.0f%%)" % pw)
        print("  P=%8.2f d  R_E=%6.1f  R_M=%6.2f  ampl %.2e  Mars power %5.1f%%  %s"
              % (P, RE[i], RM[i], aE, pw, verdict), flush=True)
        rows.append({"period_d": float(P), "RE": float(RE[i]), "RM": float(RM[i]),
                     "amp_frac": float(aE), "mars_power": pw,
                     "earth_only": bool(pw >= 95)})

    json.dump({"n_earth": len(iE), "n_mars": len(iM), "n_both": len(both),
               "thrE": float(tE), "thrM": float(tM), "rows": rows},
              open(os.path.join(D, "viewpoint_search.json"), "w"), indent=1)
    print("\nsaved ~/viewpoint_search.json", flush=True)

if __name__=="__main__": main()
