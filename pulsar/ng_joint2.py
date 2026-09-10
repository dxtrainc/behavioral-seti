"""
The corrected surrogate: phase randomisation in a basis that respects sampling.

TWO FIXES WERE TRIED. The first added a red-noise Fourier basis to the joint fit,
on the theory that the null was conservative because NANOGrav residuals are
post-fit against a noise model the projection did not replicate. It changed
nothing: median p stayed at 0.948, KS D = 0.467. That diagnosis was wrong.

THE ACTUAL CAUSE. The surrogate used np.fft.rfft, which treats samples as evenly
spaced. Pulsar TOAs are not: 98% of gaps in B1855+09 are under an hour and the
largest is 477 days, a mean/median spacing ratio of 7e10. Randomising phases in
sample-index space and mapping back onto the real times manufactures large
ramp-like excursions, so surrogates beat the data and the bar was set too high.

THE FIX. Randomise phases in the Fourier basis evaluated AT THE REAL TOA TIMES.
Fit the red-noise coefficients, rotate each sine/cosine pair by a random angle --
which preserves that mode's amplitude exactly -- rebuild on the true epochs, and
add white noise at the measured TOA errors. The surrogate then has the pulsar's
own red spectrum and its own sampling, which is what the statistic sees.
"""
import os, sys, glob, json
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG
from ng_joint import basis, scan_joint, NMODE

def fit_modes(t, r, e, A, nmode=NMODE):
    """weighted least squares for the red-noise coefficients on real epochs"""
    w = 1.0/np.maximum(e, 1e-9)**2
    Aw = A*w[:, None]
    G = A.T @ Aw
    G += np.eye(G.shape[0])*1e-12*np.trace(G)/G.shape[0]
    c = np.linalg.pinv(G) @ (Aw.T @ r)
    nM = A.shape[1] - 2*nmode
    return c[:nM], c[nM:]

def make_surrogate(t, e, A, cF, rs, nmode=NMODE):
    """rotate each mode's phase, keeping its amplitude; rebuild on real epochs"""
    nM = A.shape[1] - 2*nmode
    F = A[:, nM:]
    c = cF.copy()
    for k in range(nmode):
        s, co = c[2*k], c[2*k+1]
        amp = np.hypot(s, co)
        th = rs.uniform(0, 2*np.pi)
        c[2*k], c[2*k+1] = amp*np.sin(th), amp*np.cos(th)
    return F @ c + rs.normal(0.0, e)

def run(name, t, r, e, M, NS=200):
    A = basis(t, M)
    grid, snr, rp, P, w = scan_joint(t, r, e, A)
    obs = float(np.max(np.abs(snr)))
    _, cF = fit_modes(t, r, e, A)
    rs = np.random.default_rng(11)
    nulls = np.empty(NS)
    for k in range(NS):
        sur = make_surrogate(t, e, A, cF, rs)
        _, s2, _, _, _ = scan_joint(t, sur, e, A, ngrid=60)
        nulls[k] = np.max(np.abs(s2))
    p = (1+(nulls >= obs).sum())/(1+NS)
    return obs, float(np.percentile(nulls, 95)), float(p)

def main():
    files = sorted(glob.glob(os.path.join(NG.DIR, "*.pkl")))
    print("surrogate: phase rotation in the Fourier basis, on real TOA epochs\n", flush=True)
    print("  %-14s %6s %8s %9s %8s" % ("pulsar", "nTOA", "max|S/N|", "null 95%", "p"), flush=True)
    out = []
    for f in files:
        try: name, t, r, e, M = NG.load(f)
        except Exception: continue
        if len(t) < 500 or (t.max()-t.min()) < 5*NG.YR: continue
        try: obs, n95, p = run(name, t, r, e, M)
        except Exception as ex:
            print("  %-14s failed: %s" % (name, ex), flush=True); continue
        print("  %-14s %6d %8.2f %9.2f %8.4f %s" %
              (name, len(t), obs, n95, p, "<-- excess" if p < 0.05 else ""), flush=True)
        out.append({"name": name, "ntoa": len(t), "max_snr": obs, "null_p95": n95,
                    "p": p, "rms_us": float(r.std()*1e6)})
    json.dump(out, open(os.path.expanduser("~/ng_joint2.json"), "w"), indent=1)
    pv = np.array([o["p"] for o in out])
    from scipy import stats
    D, kp = stats.kstest(pv, "uniform")
    print("\n=== null calibration ===", flush=True)
    print("  %d pulsars   median p = %.3f (uniform expects 0.500)" % (len(pv), np.median(pv)), flush=True)
    print("  KS vs uniform: D = %.3f, p = %.4f  -> %s"
          % (D, kp, "CALIBRATED" if kp > 0.05 else "still not uniform"), flush=True)
    print("  p<0.05: %d (expect %.1f)   Bonferroni %.5f   min p %.4f"
          % ((pv < 0.05).sum(), 0.05*len(pv), 0.05/len(pv), pv.min()), flush=True)
    print("\nsaved ~/ng_joint2.json", flush=True)

if __name__ == "__main__":
    main()
