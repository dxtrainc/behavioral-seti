"""
TWIN-INSTRUMENT TEST, REDONE WITH CALENDAR-SYMMETRIC SURROGATES.

The first version (validate/wl_twin.py) printed "pipeline MANUFACTURES structure" and
that verdict was not supportable. It used INDEX-AXIS surrogates -- the construction
validate/wl_cal3.py had already measured as over-firing 9.39x on gappy data -- and arm
C carries 12,762 breaks. A large arm-C z is exactly what that miscalibration produces,
so the test could not distinguish a pipeline artefact from its own null being wrong.

validate/wl_cal4.py fixed this for the daily sweep: randomise on the FULL grid and
subsample through the real mask, so data and surrogate share a calendar and not merely
a spectrum. That fix was never propagated here. It is propagated now.

THE SECOND CAVEAT STANDS AND IS NOT FIXABLE BY BETTER SURROGATES. GOES-16 sits at 75W
and GOES-17 at 137W. Different longitudes mean different local times, different eclipse
seasons and different thermal cycles, so the DIFFERENCE of the two instruments carries
genuine spacecraft-specific structure of its own. Arm C is therefore a partial dark
channel, not a clean one, and a residual z in it does not prove the pipeline is at
fault. That is stated rather than argued away.
"""
import os, sys, json, time, itertools
import numpy as np
from scipy.stats import norm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "search"))
import wl_sweep as W
import wl_twin as T

def ns(x): return norm.ppf((np.argsort(np.argsort(x))+0.5)/len(x))

if __name__ == "__main__":
    NSUR = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    A, B, ok, ntot = T.build()
    idx = np.where(ok)[0]
    keep = np.diff(idx) == 1
    NF = ntot                                   # full one-minute grid
    SUB = idx
    print("GOES-16/17 EUVS, %d bands; full grid %d, common %d (%.1f%%), %d breaks"
          % (len(T.BANDS), NF, len(SUB), 100*len(SUB)/NF, (~keep).sum()), flush=True)

    def fullgrid(Z2d, mask):
        """put a common-sample series back on the full grid by interpolation, so the
        surrogate can be generated the way the data actually arose: a smooth process
        sampled on an irregular calendar."""
        out = np.empty((Z2d.shape[0], NF))
        xf = np.arange(NF)
        for r in range(Z2d.shape[0]):
            out[r] = np.interp(xf, mask, Z2d[r])
        return out

    def screen_cal(F_amp, F_ang, rng):
        ph = rng.uniform(0, 2*np.pi, F_amp.shape[1]); ph[0] = 0.0
        if NF % 2 == 0: ph[-1] = 0.0
        Y = np.fft.irfft(F_amp*np.exp(1j*(F_ang+ph[None, :])), n=NF, axis=1)
        return Y[:, SUB]

    def run(Z, label, seed):
        W.KEEP = keep
        X = np.vstack([ns(r) for r in Z])
        XF = fullgrid(X, SUB)
        Fx = np.fft.rfft(XF, axis=1)
        amp, ang = np.abs(Fx), np.angle(Fx)
        rng = np.random.default_rng(seed)
        obs = T.stats_of(X)
        S = np.array([T.stats_of(screen_cal(amp, ang, rng)) for _ in range(NSUR)])
        z = (obs-np.nanmean(S, 0))/np.nanstd(S, 0)
        print("  %-26s median z %9.1f   max z %10.1f   |z|>20 %3d/%d"
              % (label, np.nanmedian(z), np.nanmax(z), (np.abs(z) > 20).sum(), len(T.TESTS)), flush=True)
        return z

    t0 = time.time()
    D = np.log(A)-np.log(B); D = D - np.median(D, axis=1, keepdims=True)
    zA = run(A, "arm A  G16 (solar)", 11)
    zB = run(B, "arm B  G17 (solar)", 22)
    zC = run(D, "arm C  log G16 - log G17", 33)
    print("\n  index-axis version gave: A 2929.6, B 1061.8, C 613.9")
    print("  daily six-channel sweep, calendar-symmetric: median z 45.3")
    mA, mC = float(np.nanmedian(zA)), float(np.nanmedian(zC))
    print("\n  arm C / arm A ratio: %.3f" % (mC/mA if mA else np.nan))
    print("  NOTE: arm C is a PARTIAL dark channel -- the two spacecraft differ in")
    print("  longitude, eclipse season and thermal cycle -- so residual z in C does")
    print("  not by itself convict the pipeline.")
    json.dump({"median_z": {"A": mA, "B": float(np.nanmedian(zB)), "C": mC}, "nsur": NSUR},
              open(os.path.expanduser("~/wl_twin2.json"), "w"), indent=1)
    print("  %.0f s" % (time.time()-t0))
