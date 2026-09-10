"""
A SEARCH FOR TIMED SPIN-UP STEPS IN MILLISECOND PULSARS -- NANOGrav 15 yr.

The mechanism, from the companion paper's section 2.4: a sender delivers captured
prograde matter to a millisecond pulsar in a timed sequence, one small spin-up per
symbol. Nothing in the programme so far tests for this. The population and ladder
searches used ATNF catalogue values -- a single (P, Pdot) per pulsar, in which a
sequence of steps is invisible by construction -- and the six PTA searches asked
about cross-pulsar correlation and continuous waves, never about per-pulsar step
structure.

THE SIGNAL. A fractional frequency step dnu/nu at epoch t0 accumulates phase
linearly thereafter, so in timing residuals it is a RAMP switched on at t0:

    r(t) = -(dnu/nu) * (t - t0) * H(t - t0)

THE THING THAT MAKES IT HARD, and which decides whether any limit is honest: the
timing model has already fitted nu, nudot, position, proper motion, parallax and
binary parameters, so a ramp is PARTIALLY ABSORBED by that fit. Every template is
therefore projected into the space orthogonal to the design matrix M before it is
correlated with the residuals, exactly as a real fit would absorb it. Skipping
that step would overstate sensitivity by a large factor.

THE NULL. Phase randomisation of the residuals, which preserves the power
spectrum -- and so the red noise that dominates MSP timing -- while destroying a
phase-coherent feature. The companion paper warns that phase randomisation is a
NO-OP against power-spectrum statistics; this statistic is a time-domain matched
filter and is sensitive to phase, so the surrogate is matched to it here in the
sense of section 3.4. That distinction is the reason it is the right null.
"""
import os, sys, glob, types, pickle, json
import numpy as np

# enterprise pickles reference sksparse, which is not installed and is not needed
# to read TOAs and residuals; stub it so unpickling resolves.
_m = types.ModuleType("sksparse"); _c = types.ModuleType("sksparse.cholmod")
_c.cholesky = object; _c.CholmodError = Exception; _m.cholmod = _c
sys.modules.setdefault("sksparse", _m); sys.modules.setdefault("sksparse.cholmod", _c)

DIR = os.path.expanduser("~/ng15/ng15_psrs")
YR = 365.25*86400.0
rng = np.random.default_rng(7)

def load(path):
    o = pickle.load(open(path, "rb"))
    t = np.asarray(o.toas, float)
    r = np.asarray(o.residuals, float)
    e = np.asarray(o.toaerrs, float)
    M = np.asarray(o.Mmat, float)
    k = np.argsort(t)
    return o.name, t[k], r[k], e[k], M[k]

def projector(M, w):
    """weighted projection that removes the timing model's own span"""
    Wm = M*w[:, None]
    G = M.T @ Wm
    G += np.eye(G.shape[0])*1e-12*np.trace(G)/G.shape[0]
    Gi = np.linalg.pinv(G)
    def P(x):
        return x - M @ (Gi @ (Wm.T @ x))
    return P

def scan(t, r, e, M, ngrid=200, edge=0.1):
    """matched filter for a ramp switched on at t0, over a grid of t0"""
    w = 1.0/np.maximum(e, 1e-9)**2
    P = projector(M, w)
    rp = P(r)
    lo, hi = t.min(), t.max()
    span = hi-lo
    grid = np.linspace(lo+edge*span, hi-edge*span, ngrid)
    snr = np.zeros(ngrid)
    for i, t0 in enumerate(grid):
        T = np.where(t > t0, t-t0, 0.0)
        Tp = P(T)
        den = float(Tp @ (w*Tp))
        if den <= 0: continue
        num = float(Tp @ (w*rp))
        snr[i] = num/np.sqrt(den)
    return grid, snr, rp, P, w

def surrogate(r, rs):
    """phase randomisation: preserves the power spectrum, destroys phase"""
    F = np.fft.rfft(r)
    ph = rs.uniform(0, 2*np.pi, len(F)); ph[0] = 0
    if len(r) % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(F)*np.exp(1j*ph), n=len(r))

def inject(t, r, dnu_over_nu, t0):
    """a fractional frequency step at t0 shows in residuals as a ramp"""
    return r - dnu_over_nu*np.where(t > t0, t-t0, 0.0)

def sensitivity(name, t, r, e, M, depths=(1e-13, 3e-13, 1e-12, 3e-12, 1e-11), N=40):
    """95% recovery amplitude, measured the way section 4.2 of the main paper
    measures every limit: inject, run the identical pipeline, count."""
    _, _, rp, P, w = scan(t, r, e, M, ngrid=20)
    rs = np.random.default_rng(23)
    # threshold from this pulsar's own phase-randomised null
    nulls = []
    for _ in range(60):
        _, s2, _, _, _ = scan(t, surrogate(rp, rs), e, M, ngrid=40)
        nulls.append(np.max(np.abs(s2)))
    thr = float(np.percentile(nulls, 95))
    span = t.max()-t.min()
    out = []
    lim = None
    for d in depths:
        hits = 0
        for _ in range(N):
            base = surrogate(rp, rs)                       # red noise, no step
            t0 = t.min() + rs.uniform(0.25, 0.75)*span
            _, s2, _, _, _ = scan(t, inject(t, base, d, t0), e, M, ngrid=40)
            if np.max(np.abs(s2)) > thr: hits += 1
        rec = 100.0*hits/N
        out.append((d, rec))
        if lim is None and rec >= 95: lim = d
    return thr, out, lim

def main():
    files = sorted(glob.glob(os.path.join(DIR, "*.pkl")))
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    if only: files = [f for f in files if any(o in f for o in only)]
    print("pulsars: %d\n" % len(files), flush=True)
    print("  %-14s %6s %7s %9s %8s %8s %8s" %
          ("pulsar", "nTOA", "span/yr", "rms/us", "max|S/N|", "null 95%", "p"), flush=True)
    out = []
    for f in files:
        try:
            name, t, r, e, M = load(f)
        except Exception as ex:
            print("  %-14s load failed: %s" % (os.path.basename(f), ex), flush=True); continue
        if len(t) < 500 or (t.max()-t.min()) < 5*YR: continue
        grid, snr, rp, P, w = scan(t, r, e, M)
        obs = float(np.max(np.abs(snr)))
        NS = 200
        rs = np.random.default_rng(11)
        nulls = np.empty(NS)
        for k in range(NS):
            rsur = surrogate(rp, rs)
            _, s2, _, _, _ = scan(t, rsur, e, M, ngrid=60)
            nulls[k] = np.max(np.abs(s2))
        p = (1+(nulls >= obs).sum())/(1+NS)
        print("  %-14s %6d %7.2f %9.3f %8.2f %8.2f %8.4f %s" %
              (name, len(t), (t.max()-t.min())/YR, r.std()*1e6, obs,
               np.percentile(nulls, 95), p, "<-- excess" if p < 0.05 else ""), flush=True)
        out.append({"name": name, "ntoa": len(t), "span_yr": (t.max()-t.min())/YR,
                    "rms_us": r.std()*1e6, "max_snr": obs,
                    "null_p95": float(np.percentile(nulls, 95)), "p": float(p),
                    "best_epoch_mjd": float(grid[int(np.argmax(np.abs(snr)))]/86400.0)})
    json.dump(out, open(os.path.expanduser("~/ng_step.json"), "w"), indent=1)
    print("\nsaved ~/ng_step.json  (%d pulsars)" % len(out), flush=True)

    # sensitivity on the three best-timed pulsars in the set
    best = sorted(out, key=lambda d: d["rms_us"])[:3]
    print("\n=== sensitivity: 95% recovery in fractional frequency step ===", flush=True)
    for b in best:
        f = os.path.join(DIR, b["name"]+".pkl")
        if not os.path.exists(f): continue
        name, t, r, e, M = load(f)
        thr, curve, lim = sensitivity(name, t, r, e, M)
        print("  %-12s thr %6.1f   %s   95%% at %s" %
              (name, thr, " ".join("%.0e:%3.0f%%" % c for c in curve),
               ("%.1e" % lim) if lim else ">1e-11"), flush=True)
        b["thr"] = thr; b["curve"] = curve; b["lim95"] = lim
    json.dump(out, open(os.path.expanduser("~/ng_step.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
