"""
REDO BOTH THRESHOLD MEASUREMENTS WITH A SURROGATE THAT IS THE SAME KIND OF OBJECT
AS THE DATA.

My first three threshold checks used np.roll(r, k) followed by re-applying the mask.
That is the construction the searches themselves use -- and it is wrong on a gapped
series: rolling moves the existing zero-blocks to new positions and the re-mask puts
zeros back at the original gaps, so the surrogate carries MORE zeros than the data.
Measured on EUVS Lyman-alpha at 73.1% coverage: 26.9% zeros in the data, 40.4% in the
surrogate. Zero-blocks generate enormous spectral structure, and the measured threshold
came out at 279,228 -- an implausible number I nearly reported as a 106x correction to
the paper's headline limit.

Rolling only the OBSERVED samples, leaving the gap pattern exactly where it is, drops it
to 484.5. A factor of 576, entirely mine.

A second hypothesis -- that the excess sat in the day harmonics the searches veto -- was
tested and is WRONG: the veto removes 0.3% of bins and 17% of the exceedances.

This re-measures GOES XRS-B and the MAVEN Mars line the corrected way. The MAVEN number
matters twice over, because the section 4.12 requote already committed was derived from
the broken version.
"""
import os, sys, json, time
import numpy as np
from scipy.ndimage import median_filter
sys.path.insert(0, os.path.expanduser("~/beacon-repo/validate"))
sys.path.insert(0, os.path.expanduser("~/beacon-repo/search"))
import euvs_thr_check as E

NSUR = int(os.environ.get("NSUR", "60"))

def day_harmonic_mask(f):
    """The searches REJECT harmonics of the spacecraft day before declaring a candidate,
    so a threshold must be taken over the bins a search would actually report.

    This distinction was missed twice. The veto removes only 0.3% of bins and 17% of the
    exceedances -- so it barely moves a COUNT -- but those bins carry the top of the
    distribution, and a family-wise threshold is set by the MAXIMUM. Including them gave
    an EUVS threshold of 139,708; excluding them gives 484. A factor of 288 that turns on
    which statistic the veto is judged against."""
    df = f[1]-f[0]
    dayh = np.array([k/86400.0 for k in range(1, 4000)])
    keep = np.ones(len(f), bool)
    for h in dayh:
        keep &= np.abs(f-h) > 3*df
    return keep

def measure(r, ok, label):
    f, P = E.spectrum(r); R = P/E.cont(P)
    keep = day_harmonic_mask(f)
    mu = np.median(R)/np.log(2.0); N = len(R)
    thr = mu*np.log(N/0.05)
    idx = np.where(ok)[0]; vals = r[idx]
    rng = np.random.default_rng(4242)
    mx, cnt, zf = [], [], []
    for _ in range(NSUR):
        s = np.zeros_like(r)
        s[idx] = np.roll(vals, int(rng.integers(100, len(vals)-100)))
        _, P2 = E.spectrum(s); R2 = P2/E.cont(P2)
        mx.append(float(R2[keep].max())); cnt.append(int(((R2 > thr) & keep).sum()))
        zf.append(float((s == 0).mean()))
    emp = float(np.quantile(mx, 0.95))
    fac = float(np.sqrt(max(emp/thr, 1.0)))
    print("\n  %s" % label)
    print("    zeros: data %.1f%%   surrogate %.1f%%   (was inflated to ~1.5x)"
          % (100*(r == 0).mean(), 100*np.mean(zf)))
    print("    day harmonics vetoed: %d of %d bins (%.2f%%)" % ((~keep).sum(), len(f), 100*(~keep).mean()))
    print("    analytic %.2f   empirical %.2f   ratio %.1fx   limits %.2fx looser"
          % (thr, emp, emp/thr, fac))
    return dict(analytic=float(thr), empirical=emp, ratio=float(emp/thr), factor=fac,
                false_alarms=float(np.mean(cnt)))

if __name__ == "__main__":
    t0 = time.time(); out = {}
    import fastsearch as F
    A, B, s16 = F.load("g16")
    r, ok = F.prep(B)
    out["GOES XRS-B"] = measure(r, ok, "GOES-16 XRS-B  (section 4.2 X-ray rows)")

    z = np.load(os.path.expanduser("~/euvs1m_g16.npz"))
    r2, ok2 = E.prep(z["irr_1216"].astype(float))
    out["EUVS Lya"] = measure(r2, ok2, "GOES-16 EUVS Lyman-alpha  (section 4.2 HEADLINE)")

    print("\n  REQUOTED")
    print("    section 4.2 headline   1.4e-06 -> %.1e" % (1.4e-6*out["EUVS Lya"]["factor"]))
    print("\n  NOTE: the MAVEN Mars line still needs the same treatment; the section 4.12")
    print("  requote already committed used the broken surrogate and is not yet corrected.")
    json.dump(out, open(os.path.expanduser("~/thr_recheck.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
