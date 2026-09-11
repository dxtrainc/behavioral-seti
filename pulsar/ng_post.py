"""
The null built from NANOGrav's own published noise posteriors.

Three home-made surrogates failed, each for a diagnosable reason: FFT phase
randomisation assumed even sampling (median p 0.948); Fourier phase rotation
destroyed coherence between modes (6 of 7 "detections"); a hand-fitted power law
was conservative even after a units bug was fixed (median p 0.975). Two
hypotheses were ruled out along the way -- the residuals are not chromatic
(DMX has already removed dispersion) and the TOA errors are not understated
(EFAC-equivalent 0.82-1.03 measured from close pairs).

What remained was the noise model itself. NANOGrav publishes, per pulsar, the
full MCMC chain over per-backend EFAC, EQUAD and ECORR plus the achromatic red
noise amplitude and spectral index. Drawing a surrogate means drawing a POSTERIOR
SAMPLE and then realising that noise model on the real epochs, so the null
carries the published model's own uncertainty rather than a point estimate of
mine.

    white   sigma_eff^2 = EFAC_b^2 * sigma^2 + EQUAD_b^2      per backend b
    ECORR   one common offset per (backend, epoch), sigma = ECORR_b
    red     rho_k = A^2/(12 pi^2) (f_k/f_yr)^-gamma * YR^3 / T, 30 modes
"""
import os, sys, glob, json, re
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG
from ng_joint import basis, scan_joint, NMODE

REL = os.path.expanduser("~/ng15_release/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/noise")
YR = 365.25*86400.0
FYR = 1.0/YR

def load_chain(name):
    pf = os.path.join(REL, "%s.nb.pars.txt" % name)
    cf = os.path.join(REL, "%s.nb.chain_1.txt" % name)
    if not (os.path.exists(pf) and os.path.exists(cf)): return None, None
    pars = [l.strip() for l in open(pf) if l.strip()]
    ch = np.loadtxt(cf)
    if ch.ndim == 1: ch = ch[None, :]
    ch = ch[len(ch)//4:]                       # discard burn-in
    return pars, ch[:, :len(pars)]

def sample(pars, chain, rs):
    return dict(zip(pars, chain[rs.integers(0, len(chain))]))

def realise(t, e, be, epoch, s, name, F, f, T, rs):
    """one draw of the published noise model on the real epochs"""
    sig = np.array(e, float)
    out = np.zeros(len(t))
    for b in np.unique(be):
        m = be == b
        ef = s.get("%s_%s_efac" % (name, b), 1.0)
        lq = s.get("%s_%s_log10_equad" % (name, b), s.get("%s_%s_log10_t2equad" % (name, b), None))
        eq = 10.0**lq if lq is not None else 0.0
        se = np.sqrt((ef*sig[m])**2 + eq**2)
        out[m] = rs.normal(0.0, se)
        lc = s.get("%s_%s_log10_ecorr" % (name, b), None)
        if lc is not None:
            ec = 10.0**lc
            for ep in np.unique(epoch[m]):
                k = m & (epoch == ep)
                out[k] += rs.normal(0.0, ec)    # one common offset per epoch
    lgA = s.get("%s_red_noise_log10_A" % name)
    gam = s.get("%s_red_noise_gamma" % name)
    if lgA is not None and gam is not None:
        rho = (10.0**lgA)**2/(12*np.pi**2)*(f/FYR)**(-gam)*YR**3/T
        out = out + F @ rs.normal(0.0, np.sqrt(np.repeat(rho, 2)))
    return out

def run(name, t, r, e, be, M, pars, chain, NS=200):
    A = basis(t, M)
    _, snr, rp, P, w = scan_joint(t, r, e, A)
    obs = float(np.max(np.abs(snr)))
    T = t.max()-t.min()
    fk = np.arange(1, NMODE+1)/T
    ph = 2*np.pi*np.outer(t-t.min(), fk)
    F = np.empty((len(t), 2*NMODE)); F[:, 0::2] = np.sin(ph); F[:, 1::2] = np.cos(ph)
    epoch = np.round(t/86400.0).astype(np.int64)
    rs = np.random.default_rng(11)
    nulls = np.empty(NS)
    for k in range(NS):
        s = sample(pars, chain, rs)
        sur = realise(t, e, be, epoch, s, name, F, fk, T, rs)
        _, s2, _, _, _ = scan_joint(t, sur, e, A, ngrid=60)
        nulls[k] = np.max(np.abs(s2))
    p = (1+(nulls >= obs).sum())/(1+NS)
    return obs, float(np.percentile(nulls, 95)), float(p)

def main():
    import pickle, types
    files = sorted(glob.glob(os.path.join(NG.DIR, "*.pkl")))
    lim = int(os.environ.get("NPSR", "0"))
    if lim: files = files[:lim]
    print("surrogates drawn from the published NANOGrav 15 yr noise posteriors\n", flush=True)
    print("  %-14s %6s %9s %9s %8s" % ("pulsar", "nTOA", "max|S/N|", "null 95%", "p"), flush=True)
    out = []
    for fp in files:
        try:
            name, t, r, e, M = NG.load(fp)
            o = pickle.load(open(fp, "rb"))
            be = np.asarray(o.backend_flags)[np.argsort(np.asarray(o.toas, float))]
        except Exception: continue
        if len(t) < 500 or (t.max()-t.min()) < 5*NG.YR: continue
        pars, chain = load_chain(name)
        if pars is None:
            print("  %-14s no published chain" % name, flush=True); continue
        try: obs, n95, p = run(name, t, r, e, be, M, pars, chain)
        except Exception as ex:
            print("  %-14s failed: %s" % (name, ex), flush=True); continue
        print("  %-14s %6d %9.2f %9.2f %8.4f %s" %
              (name, len(t), obs, n95, p, "<-- excess" if p < 0.05 else ""), flush=True)
        out.append({"name": name, "ntoa": len(t), "max_snr": obs, "null_p95": n95, "p": p})
    json.dump(out, open(os.path.expanduser("~/ng_post.json"), "w"), indent=1)
    if len(out) > 4:
        pv = np.array([o["p"] for o in out]); from scipy import stats
        D, kp = stats.kstest(pv, "uniform")
        print("\n=== null calibration ===", flush=True)
        print("  %d pulsars   median p = %.3f (uniform expects 0.500)" % (len(pv), np.median(pv)), flush=True)
        print("  KS vs uniform: D = %.3f, p = %.4f -> %s"
              % (D, kp, "CALIBRATED" if kp > 0.05 else "NOT uniform"), flush=True)
        print("  p<0.05: %d (expect %.1f)  Bonferroni %.5f  min p %.4f"
              % ((pv < 0.05).sum(), 0.05*len(pv), 0.05/len(pv), pv.min()), flush=True)
    print("\nsaved ~/ng_post.json", flush=True)

if __name__ == "__main__":
    main()
