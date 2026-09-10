import sys, os
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG
from ng_joint import basis, scan_joint, NMODE
from ng_joint2 import fit_modes, make_surrogate

name, t, r, e, M = NG.load(os.path.join(NG.DIR, "B1855+09.pkl"))
A = basis(t, M)
g, s, rp, P, w = scan_joint(t, r, e, A)
_, cF = fit_modes(t, r, e, A)
rs = np.random.default_rng(1)
sur = [make_surrogate(t, e, A, cF, rs) for _ in range(20)]
sur_p = [P(x) for x in sur]
print("real projected residual rms   : %.4e s" % rp.std())
print("surrogate rms, pre-projection : %.4e s" % np.mean([x.std() for x in sur]))
print("surrogate rms, post-projection: %.4e s" % np.mean([x.std() for x in sur_p]))
print("ratio real / surrogate(post)  : %.2f" % (rp.std()/np.mean([x.std() for x in sur_p])))
print()
print("mean TOA error (white floor)  : %.4e s" % e.mean())
print("-> the 30-mode + white model reproduces %.0f%% of the residual power"
      % (100*np.mean([x.std() for x in sur_p])/rp.std()))
