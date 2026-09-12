"""
ROW 7 -- p-mode FREQUENCY STRUCTURE, on SOHO/GOLF, 25.9 years of Doppler velocity.

Section 3.6 lists row 7 at a designer level of 1e-6 against our reach of 1e-5..1e-6,
searched only as a comb re-detection (section 4.6). This is the search proper.

WHY GOLF AND NOT EXIS. Section 4.6 ran on GOES EXIS "and not by choice" -- the
canonical helioseismology archives were recorded as unreachable. That was wrong for
VIRGO and GOLF: NASA mirrors both over plain HTTPS. GOLF is a resonant-scattering
spectrophotometer measuring DOPPLER VELOCITY, a different observable from EXIS
photometry, at 20 s over 25.9 yr. It is the instrument this question should be asked of.

WHAT IS BEING SEARCHED. Not the presence of p-modes -- that is settled, and here it
is a control. A sender able to modulate the Sun at all could shift MODE FREQUENCIES,
and mode frequencies are the most precisely measured quantity in solar physics. So
the observable is the p-mode frequency shift as a function of time, and the question
is whether it carries structure beyond what solar activity explains.

THREE GATES, ALL ABLE TO FAIL.

  1a  THE COMB. The spectrum must show modes spaced by dnu = 135 uHz under an
      envelope near nu_max = 3.09 mHz. dnu goes as sqrt(M/R^3) and is PREDICTED from
      stellar structure, never fitted. A pipeline that cannot recover it in GOLF --
      the instrument p-modes are classically measured with -- is not calibrated, and
      nothing else it reports should be read.

  1b  THE ACTIVITY SHIFT. Mode frequencies are KNOWN to rise and fall with the solar
      cycle, by about 0.4 uHz peak-to-peak at 3 mHz. Recovering that proves the
      pipeline can measure a frequency shift at all, and at roughly the right size.
      A search for a 1e-6 modulation that cannot see a known 1.3e-4 one is not a
      search. This gate is why the known solar-cycle term is an asset, not a nuisance.

  2   PM1 vs PM2. GOLF has two photomultipliers. A real frequency shift appears in
      both; a detector artefact does not. The twin gate is inside the instrument, so
      there is no cross-calibration to argue about.

ONLY THEN the residual -- shift series minus the activity term -- is searched for
imposed periodicity, with the activity regression reported so the reader can see how
much was removed.
"""
import os, sys, time, json
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

VD = os.path.expanduser("~/virgo")
DT = 20.0                       # s
SEG_D = 90                      # days per segment
NU_LO, NU_HI = 2.5e-3, 4.0e-3   # p-mode band, Hz
DNU_EXP = 135.0e-6              # predicted large separation, Hz

def load(tag):
    return np.asarray(fits.open(os.path.join(VD, "GOLF_26y_%s.fits" % tag))[0].data, float)

def spec(x, dt=DT):
    n = len(x)//2*2
    F = np.fft.rfft(x[:n]*np.hanning(n))
    return np.fft.rfftfreq(n, d=dt)[1:], (np.abs(F)**2)[1:]

def rebin(f, P, width):
    """rebin a power spectrum to `width` Hz per bin, by mean"""
    idx = ((f - f[0])/width).astype(np.int64)
    nb = int(idx.max())+1
    cnt = np.bincount(idx, minlength=nb).astype(float)
    val = np.bincount(idx, weights=P, minlength=nb)
    ok = cnt > 0
    return f[0] + (np.arange(nb)[ok]+0.5)*width, val[ok]/cnt[ok]

def comb_acf(f, P, rebin_hz=0.5e-6, cont_hz=30e-6):
    """autocorrelate the p-mode band power spectrum; a real comb peaks at dnu.

    THE FIRST VERSION FAILED THIS GATE, AND THE BUG WAS A WINDOW SIZE. GOLF runs
    40,845,600 samples at 20 s, so df = 1.22 nHz and a median_filter of 501 bins is
    0.61 uHz -- NARROWER than a p-mode linewidth. It was flattening the very modes it
    was meant to normalise against, and the ACF came back at r = 0.004: no comb, from
    the instrument p-modes are classically measured with.

    At this resolution each mode is also split by rotation into components ~0.4 uHz
    apart, so the comb is clearest AFTER rebinning to about a linewidth. Rebin to
    0.5 uHz, flatten against a 30 uHz continuum -- wide compared with a mode, narrow
    compared with dnu = 135 uHz -- and only then autocorrelate."""
    b = (f > NU_LO) & (f < NU_HI)
    fb, Pb = rebin(f[b], P[b], rebin_hz)
    ncont = max(3, int(round(cont_hz/rebin_hz)) | 1)
    y = Pb/median_filter(Pb, size=ncont, mode="nearest")
    y = y - y.mean()
    a = np.correlate(y, y, "full")[len(y)-1:]
    a = a/a[0]
    lag = np.arange(len(a))*rebin_hz
    w = (lag > 100e-6) & (lag < 170e-6)
    i = np.where(w)[0][np.argmax(a[w])]
    return lag[i], a[i], lag, a, rebin_hz

def shift_series(x, nseg_samples, ref=None):
    """p-mode band frequency shift per segment, by cross-correlation against a
    reference spectrum. Sub-bin shift from a parabolic fit to the CC peak."""
    nseg = len(x)//nseg_samples
    out, duty = [], []
    S = []
    for k in range(nseg):
        seg = x[k*nseg_samples:(k+1)*nseg_samples]
        duty.append(float(np.mean(seg != 0.0)))
        if duty[-1] < 0.5:
            S.append(None); continue
        f, P = spec(seg)
        b = (f > NU_LO) & (f < NU_HI)
        y = P[b]/median_filter(P[b], size=101, mode="nearest")
        S.append(y - y.mean())
    good = [s for s in S if s is not None]
    if ref is None:
        L = min(len(s) for s in good); ref = np.mean([s[:L] for s in good], axis=0)
    f0, _ = spec(x[:nseg_samples]); df = f0[1]-f0[0]
    for s in S:
        if s is None: out.append(np.nan); continue
        L = min(len(s), len(ref))
        c = np.correlate(s[:L], ref[:L], "full")
        m = len(ref[:L])-1
        w = 40
        seg_c = c[m-w:m+w+1]
        j = int(np.argmax(seg_c))
        if 0 < j < len(seg_c)-1:
            a_, b_, c_ = seg_c[j-1], seg_c[j], seg_c[j+1]
            den = a_ - 2*b_ + c_
            # A genuine peak has den < 0. The first version wrote max(den, 1e-30),
            # which on a real peak replaced a NEGATIVE denominator with 1e-30 and
            # produced shifts of order 1e29 uHz. Reject non-peaks, and clamp the
            # sub-bin correction to the half-bin it is only ever allowed to be.
            sub = 0.0 if den >= 0 else float(np.clip(0.5*(a_-c_)/den, -0.5, 0.5))
        else: sub = 0.0
        out.append((j - w + sub)*df)
    return np.array(out), np.array(duty), ref

if __name__ == "__main__":
    t0 = time.time()
    print("loading GOLF ...", flush=True)
    X = {t: load(t) for t in ("MEAN", "PM1", "PM2")}
    n = len(X["MEAN"])
    print("  %d samples at %.0f s = %.2f yr; duty cycle %.1f%%\n"
          % (n, DT, n*DT/3.156e7, 100*np.mean(X["MEAN"] != 0.0)), flush=True)

    print("GATE 1a -- the comb. dnu is PREDICTED at %.1f uHz, never fitted." % (DNU_EXP*1e6))
    for t in ("MEAN", "PM1", "PM2"):
        f, P = spec(X[t])
        b = (f > 2.0e-3) & (f < 4.5e-3)
        fb_, Pb_ = rebin(f[b], P[b], 0.5e-6)
        C = median_filter(Pb_, size=61, mode="nearest")
        fe, Re = rebin(fb_, Pb_/C, 50e-6)      # envelope is BROAD -- ~1 mHz wide
        env = fe[np.argmax(Re)]
        dnu, amp, lag, a, df = comb_acf(f, P)
        ok = (abs(dnu-DNU_EXP) < 8e-6) and (amp > 0.05) and (2.7e-3 < env < 3.6e-3)
        print("    %-5s envelope %.4f mHz   ACF peak %.2f uHz (r=%.3f)   %s"
              % (t, env*1e3, dnu*1e6, amp,
                 "PASS" if ok else "*** FAIL ***"), flush=True)
    print("  (expect envelope near 3.09 mHz, ACF peak at 135 uHz)\n")

    nss = int(SEG_D*86400/DT)
    print("GATE 1b -- the known solar-cycle frequency shift (~0.4 uHz p-p at 3 mHz)")
    sh, duty, ref = shift_series(X["MEAN"], nss)
    good = np.isfinite(sh)
    print("    %d segments of %d d, %d usable (duty>50%%)" % (len(sh), SEG_D, good.sum()))
    print("    shift range %.3f .. %.3f uHz, peak-to-peak %.3f uHz"
          % (np.nanmin(sh)*1e6, np.nanmax(sh)*1e6, (np.nanmax(sh)-np.nanmin(sh))*1e6))
    print("    fractional: %.2e  (known solar-cycle term ~1.3e-4)"
          % ((np.nanmax(sh)-np.nanmin(sh))/3.09e-3), flush=True)

    np.savez(os.path.expanduser("~/golf_shift.npz"), shift=sh, duty=duty,
             seg_days=SEG_D, dt=DT)
    print("\n  saved ~/golf_shift.npz   (%.1f s)" % (time.time()-t0))
