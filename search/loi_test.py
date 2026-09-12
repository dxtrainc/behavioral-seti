"""
VIRGO/LOI -- THE ONE SPATIALLY RESOLVED CHANNEL IN THIS PROGRAMME.

Every other search here is disc-integrated. LOI gives 16 channels -- 12 science pixels
across the solar disc plus 4 guiding references -- at 60 s over 29.0 years with a 97.19%
duty cycle, which is the best-sampled record in the paper.

WHY IT IS WORTH A SEARCH THE OTHERS ARE NOT. A feature covering fraction f of the disc
at depth d contributes f*d to a disc-integrated series and d to the pixel that contains
it. For a feature the size of one pixel that is a signal gain of 12 against a noise cost
of sqrt(12) -- net 3.46x -- which the pre-acquisition bound put at 5.8e-8 for a localised
source against the 2.0e-7 measured on SPM. It is also 3.46x WORSE than SPM for a global
modulation, because the photons have been divided for nothing. The gain is conditional
on localisation and the write-up must say so.

AND ONE TEST NO DISC-INTEGRATED CHANNEL CAN DO AT ANY SENSITIVITY. A localised feature
ROTATES: it crosses the disc in about 13 days, entering and leaving pixels in a sequence
fixed by geometry. Every pixel pair therefore has a PREDICTED lag. The control is free,
because solar rotation is already in this data -- active regions do exactly this -- so
the known rotational lag pattern is the positive control, in the same way the recovered
Earth-Mars geometry validated section 4.12.

CAVEAT CARRIED FROM THE HEADER. L2 states correction for orbit, outliers, attractors,
roll sensitivity and degradation. The outlier step is the same one section 4.15 flags for
SPM: it may remove exactly the impulsive structure a kurtosis statistic looks for, and
testing that needs L1 data the mission bundles do not carry.
"""
import os, sys, json, time
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
sys.path.insert(0, os.path.expanduser("~/beacon-repo"))
from beacon.controls import denominator_safe, recovers_known

VD = os.path.expanduser("~/virgo")
FN = "VIRGO-LOI-ALL-PIXELS-LEVEL2-19960401-20250331_V01.fits"
DT = 60.0
NU_MAX = 3.09e-3

def load():
    """FITS stores NAXIS1 fastest, so the array arrives as (NAXIS2, NAXIS1) =
    (16, 15252480) -- CHANNELS FIRST. The first version of this script wrote
    N, NCH = D.shape and got 16 samples of 15 million channels, then began looping
    over them. Transposed once here, at the boundary, so nothing downstream has to
    remember which way round it is."""
    d = np.asarray(fits.open(os.path.join(VD, FN))[0].data, float)
    if d.shape[0] < d.shape[1]:
        d = d.T                                   # -> (N samples, 16 channels)
    return d

def prep(x, win=1441):
    ok = np.isfinite(x)
    y = np.where(ok, x, 0.0)
    tr = median_filter(y, size=win, mode="nearest")
    return np.where(ok, y - tr, 0.0), ok             # ppm already relative: SUBTRACT

def spec(r):
    n = len(r)//2*2
    F = np.fft.rfft(r[:n]*np.hanning(n))
    return np.fft.rfftfreq(n, d=DT)[1:], (np.abs(F)**2)[1:]

def cont(P, w=801):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

def rebin(f, P, width):
    idx = ((f-f[0])/width).astype(np.int64); nb = int(idx.max())+1
    c = np.bincount(idx, minlength=nb).astype(float)
    v = np.bincount(idx, weights=P, minlength=nb)
    ok = c > 0
    return f[0]+(np.arange(nb)[ok]+0.5)*width, v[ok]/c[ok]

if __name__ == "__main__":
    t0 = time.time()
    D = load()
    N, NCH = D.shape
    print("VIRGO/LOI: %d samples x %d channels at %.0f s = %.1f yr"
          % (N, NCH, DT, N*DT/3.156e7), flush=True)
    finite = np.isfinite(D).mean(axis=0)
    sd = np.array([np.nanstd(D[:, k]) for k in range(NCH)])
    print("\n  ch  finite   rms(ppm)   role")
    for k in range(NCH):
        print("  %2d  %5.1f%%  %9.1f" % (k, 100*finite[k], sd[k]), flush=True)

    # science pixels: the twelve with comparable scatter. The four guiding channels
    # differ in level or variance; identified from the data rather than assumed, and
    # the split is printed so it can be checked against Appourchaux et al. 1997.
    med_sd = np.median(sd)
    sci = [k for k in range(NCH) if 0.2*med_sd < sd[k] < 5*med_sd and finite[k] > 0.5]
    print("\n  science-pixel candidates by scatter: %s (%d of %d)"
          % (sci, len(sci), NCH), flush=True)

    R = {}
    for k in sci:
        R[k], _ = prep(D[:, k])
    print("  prepped %d channels in %.0f s" % (len(sci), time.time()-t0), flush=True)

    print("\nGATE 1 -- LOI IS A HELIOSEISMOLOGY IMAGER: p-modes must be present")
    env = []
    for k in sci:
        f, P = spec(R[k])
        b = (f > 2.0e-3) & (f < 4.5e-3)
        fb, Pb = rebin(f[b], P[b], 0.5e-6)
        Cb = median_filter(Pb, size=61, mode="nearest")
        fe, Re = rebin(fb, Pb/Cb, 50e-6)
        j = int(np.argmax(Re)); env.append(fe[j])
        print("    ch %2d  envelope %.4f mHz (R=%.2f)  %s"
              % (k, fe[j]*1e3, Re[j], "pass" if 2.7e-3 < fe[j] < 3.6e-3 else "FAIL"), flush=True)
    npass = sum(1 for e in env if 2.7e-3 < e < 3.6e-3)
    print("  %d of %d pixels show the five-minute oscillation" % (npass, len(sci)))

    print("\nROTATION -- the control that is free, and the test no disc-integrated channel has")
    # daily means: rotation lives at 27 d, far above the 1-day detrend, so use raw ppm
    day = N//1440
    dm = np.array([[np.nanmean(D[i*1440:(i+1)*1440, k]) for i in range(day)] for k in sci])
    dm = np.where(np.isfinite(dm), dm, 0.0)
    dm = dm - median_filter(dm, size=(1, 181), mode="nearest")
    lags = {}
    for a in range(len(sci)):
        for b in range(a+1, len(sci)):
            x, y = dm[a], dm[b]
            x = (x-x.mean())/(x.std() or 1); y = (y-y.mean())/(y.std() or 1)
            c = np.correlate(x, y, "full")/len(x)
            m = len(x)-1; w = 20
            seg = c[m-w:m+w+1]
            lags[(sci[a], sci[b])] = (int(np.argmax(seg))-w, float(seg.max()))
    good = {k: v for k, v in lags.items() if v[1] > 0.3}
    print("  %d of %d pixel pairs correlate above r = 0.3" % (len(good), len(lags)))
    if good:
        ls = np.array([v[0] for v in good.values()])
        print("  lag distribution: median %.1f d, range %d to %d d"
              % (np.median(ls), ls.min(), ls.max()))
        print("  a rotating localised source gives NON-ZERO lags set by pixel separation;")
        print("  a global modulation gives lag 0 for every pair.")
        print("  pairs with |lag| >= 1 d: %d of %d" % (int((np.abs(ls) >= 1).sum()), len(ls)))
    json.dump({"n": int(N), "nch": int(NCH), "sci": sci,
               "envelope_mHz": [float(e*1e3) for e in env], "gate1_pass": int(npass),
               "rot_pairs": len(good), "rot_total": len(lags)},
              open(os.path.expanduser("~/loi_test.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
