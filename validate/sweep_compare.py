"""
COMPARE THE RESEEDED SWEEPS AGAINST THE COMMITTED ONES.

The seeds changed from Python's per-process hash() to a stable BLAKE2b, so the surrogate
draws are different. The p-values SHOULD differ in detail and AGREE in distribution:
they are two realisations of the same null, not two different nulls.

Three outcomes and only one is a problem:

  distributions agree, same verdict      expected. The conclusions hold and become
                                         reproducible for the first time.
  distributions agree, a marginal test
  crosses threshold in one and not the   also expected, and worth reporting. With 1,074
  other                                  and 10,996 tests at a family-wise 0.05, the
                                         cases near threshold are exactly where a single
                                         surrogate realisation is not decisive.
  DISTRIBUTIONS DISAGREE                 the seed change did something other than redraw
                                         surrogates, and neither set should be trusted
                                         until that is understood.

A Kolmogorov-Smirnov test on the two p-value sets answers the third directly.
"""
import os, json, sys
import numpy as np
from scipy.stats import ks_2samp

R = os.path.expanduser("~/beacon-repo/results")

def load(fn):
    p = os.path.join(R, fn)
    if not os.path.exists(p): return None
    return json.load(open(p))

for old_fn, new_fn, label, thr in (("sweepT_pairs_zeus.json", "sweepT_pairs_reseed.json",
                                    "PAIR SWEEP", 0.05/1074),
                                   ("sweepT_triples.json", "sweepT_triples_reseed.json",
                                    "TRIPLE SWEEP", 0.05/10996)):
    o, n = load(old_fn), load(new_fn)
    if o is None or n is None:
        print("\n%s: missing (%s / %s)" % (label, o is not None, n is not None)); continue
    po = np.array([r["p"] for r in o["results"]])
    pn = np.array([r["p"] for r in n["results"]])
    print("\n%s" % label)
    print("  tests            %d committed, %d reseeded" % (len(po), len(pn)))
    print("  median p         %.4f  ->  %.4f" % (np.median(po), np.median(pn)))
    print("  fraction < 0.05  %.4f  ->  %.4f" % ((po < 0.05).mean(), (pn < 0.05).mean()))
    print("  min p            %.5f  ->  %.5f" % (po.min(), pn.min()))
    print("  survivors at Bonferroni %.2e:  %d  ->  %d"
          % (thr, int((po < thr).sum()), int((pn < thr).sum())))
    if len(po) == len(pn):
        ks = ks_2samp(po, pn)
        print("  two-sample KS    D = %.4f, p = %.4f  -> %s"
              % (ks.statistic, ks.pvalue,
                 "same distribution" if ks.pvalue > 0.05 else "*** DISTRIBUTIONS DIFFER ***"))
    # which tests moved across threshold
    if len(po) == len(pn):
        ко = po < thr; кn = pn < thr
        moved = np.where(ко != кn)[0]
        if len(moved):
            print("  tests crossing threshold in one run only: %d" % len(moved))
            for i in moved[:6]:
                r = n["results"][i]
                print("     %-40s %-12s  %.5f -> %.5f"
                      % (r["combo"][:40], r["form"], po[i], pn[i]))
        else:
            print("  no test crosses the threshold in one run and not the other")
