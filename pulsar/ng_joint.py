"""
Joint fit of the step with the ephemeris AND a red-noise basis.

WHY THE FIRST VERSION'S NULL WAS CONSERVATIVE. Its p-values had median 0.948
against an expected 0.500 (KS D = 0.47, p < 1e-3): the real residuals carried
LESS ramp-like structure than their own phase-randomised surrogates. The cause is
not the projection, which is already equivalent to re-fitting the linearised
timing model. It is that NANOGrav residuals are post-fit against the timing model
*and a red-noise model*, while the projection removed only the former. Surrogates
therefore kept low-frequency structure that the real data had already had
absorbed, and the bar was set too high.

THE FIX. Fit everything together. The design matrix becomes

    [ M | F | T(t0) ]

with M the linearised timing model, F a Fourier basis of NMODE sine/cosine pairs
at k/T (the standard red-noise basis in pulsar timing), and T the ramp. Because
the same basis is applied to data and surrogate alike, the two become comparable,
and the amplitude estimate for T is the one that survives simultaneous
marginalisation over ephemeris and red noise -- which is what a real analysis
would report.

By Frisch-Waugh-Lovell the amplitude is identical to projecting T orthogonal to
[M|F] and correlating; the projection form is used because it is cheap inside a
scan over t0.
"""
import os, sys, glob, json
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG

NMODE = 30          # 30 Fourier pairs: the PTA convention

def basis(t, M, nmode=NMODE):
    T = t.max()-t.min()
    f = np.arange(1, nmode+1)/T
    ph = 2*np.pi*np.outer(t-t.min(), f)
    F = np.empty((len(t), 2*nmode))
    F[:, 0::2] = np.sin(ph); F[:, 1::2] = np.cos(ph)
    return np.hstack([M, F])

def scan_joint(t, r, e, A, ngrid=200, edge=0.1):
    w = 1.0/np.maximum(e, 1e-9)**2
    P = NG.projector(A, w)
    rp = P(r)
    lo, hi = t.min(), t.max(); span = hi-lo
    grid = np.linspace(lo+edge*span, hi-edge*span, ngrid)
    snr = np.zeros(ngrid)
    for i, t0 in enumerate(grid):
        Tp = P(np.where(t > t0, t-t0, 0.0))
        den = float(Tp @ (w*Tp))
        if den <= 0: continue
        snr[i] = float(Tp @ (w*rp))/np.sqrt(den)
    return grid, snr, rp, P, w

def run(name, t, r, e, M, NS=200):
    A = basis(t, M)
    grid, snr, rp, P, w = scan_joint(t, r, e, A)
    obs = float(np.max(np.abs(snr)))
    rs = np.random.default_rng(11)
    nulls = np.empty(NS)
    for k in range(NS):
        _, s2, _, _, _ = scan_joint(t, NG.surrogate(rp, rs), e, A, ngrid=60)
        nulls[k] = np.max(np.abs(s2))
    p = (1+(nulls >= obs).sum())/(1+NS)
    return obs, float(np.percentile(nulls, 95)), float(p), A, rp

def main():
    files = sorted(glob.glob(os.path.join(NG.DIR, "*.pkl")))
    print("joint fit: timing model + %d Fourier pairs + step\n" % NMODE, flush=True)
    print("  %-14s %6s %8s %9s %8s" % ("pulsar", "nTOA", "max|S/N|", "null 95%", "p"), flush=True)
    out = []
    for f in files:
        try: name, t, r, e, M = NG.load(f)
        except Exception: continue
        if len(t) < 500 or (t.max()-t.min()) < 5*NG.YR: continue
        try: obs, n95, p, A, rp = run(name, t, r, e, M)
        except Exception as ex:
            print("  %-14s failed: %s" % (name, ex), flush=True); continue
        print("  %-14s %6d %8.2f %9.2f %8.4f %s" %
              (name, len(t), obs, n95, p, "<-- excess" if p < 0.05 else ""), flush=True)
        out.append({"name": name, "ntoa": len(t), "max_snr": obs,
                    "null_p95": n95, "p": p, "rms_us": float(r.std()*1e6)})
    json.dump(out, open(os.path.expanduser("~/ng_joint.json"), "w"), indent=1)
    pv = np.array([o["p"] for o in out])
    from scipy import stats
    D, kp = stats.kstest(pv, "uniform")
    print("\n=== null calibration ===", flush=True)
    print("  %d pulsars   median p = %.3f (uniform expects 0.500)" % (len(pv), np.median(pv)), flush=True)
    print("  KS vs uniform: D = %.3f, p = %.4f  -> %s"
          % (D, kp, "CALIBRATED" if kp > 0.05 else "still not uniform"), flush=True)
    print("  p<0.05: %d (expect %.1f)   Bonferroni %.5f   min p %.4f"
          % ((pv < 0.05).sum(), 0.05*len(pv), 0.05/len(pv), pv.min()), flush=True)
    print("\nsaved ~/ng_joint.json", flush=True)

if __name__ == "__main__":
    main()
