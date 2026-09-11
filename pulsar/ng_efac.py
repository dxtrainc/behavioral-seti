"""Is the white noise larger than the quoted TOA errors?

NANOGrav fits a per-backend EFAC (error scale) and EQUAD (added variance),
because raw TOA uncertainties understate the true white noise. If a noise fit
uses the raw errors, the red-noise component absorbs the excess and surrogates
come out too structured -- which is the conservative direction, and is what the
p-values piling near 1.0 look like.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG
from ng_joint import basis, scan_joint

for nm in ("J1455-3330", "J1600-3053", "B1855+09", "J1909-3744"):
    f = os.path.join(NG.DIR, nm+".pkl")
    if not os.path.exists(f): continue
    name, t, r, e, M = NG.load(f)
    A = basis(t, M)
    _, _, rp, P, w = scan_joint(t, r, e, A, ngrid=5)
    # high-frequency residual power: differences of TOAs close in time are
    # dominated by white noise, so their scatter estimates the true white level
    o = np.argsort(t)
    dt = np.diff(t[o]); close = dt < 3600.0          # within an hour
    if close.sum() < 50:
        print("%-12s too few close pairs" % name); continue
    d = np.diff(rp[o])[close]
    white_meas = d.std()/np.sqrt(2)
    white_quoted = np.sqrt(np.mean(e**2))
    print("%-12s quoted white %7.3f us   measured from close pairs %7.3f us   EFAC-equivalent %.2f"
          % (name, white_quoted*1e6, white_meas*1e6, white_meas/white_quoted))
