"""Is the residual structure chromatic? If so, an achromatic power-law red-noise
model is missing the dominant term and no surrogate built from it can calibrate."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG
from ng_joint import basis, scan_joint

for nm in ("B1855+09", "B1937+21", "J1909-3744"):
    f = os.path.join(NG.DIR, nm+".pkl")
    if not os.path.exists(f): continue
    name, t, r, e, M = NG.load(f)
    import pickle
    o = pickle.load(open(f, "rb"))
    nu = np.asarray(o.freqs, float)[np.argsort(np.asarray(o.toas, float))]
    A = basis(t, M)
    _, _, rp, P, w = scan_joint(t, r, e, A, ngrid=5)
    # split by observing frequency: chromatic noise makes the two halves differ
    med = np.median(nu)
    lo, hi = nu < med, nu >= med
    # correlation of residual with 1/nu^2, the dispersion signature
    x = 1.0/nu**2; x = x - x.mean()
    c = float(np.corrcoef(rp, x)[0, 1])
    print("%-12s nu %4.0f-%4.0f MHz  rms(low nu) %7.3f us  rms(high nu) %7.3f us  ratio %.2f   corr(resid, 1/nu^2) = %+.3f"
          % (name, nu.min(), nu.max(), rp[lo].std()*1e6, rp[hi].std()*1e6,
             rp[lo].std()/max(rp[hi].std(), 1e-30), c))
