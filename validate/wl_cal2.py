"""CALIBRATION OF THE EXACT SWEEP PIPELINE -- gaps, forms, statistic and all.

The first calibration tested the surrogate construction on a low-gap quadruple.
This one tests what will actually be run: the six-way common dataset with 1,064
breaks, gap-spanning increments dropped, all three four-body forms.

Two questions, and the second is the one that has caught every bad null in this
programme:
  1. Is p uniform on data built to contain NO four-way structure?
  2. Does the gap masking actually work -- i.e. if the breaks are left IN, does the
     null visibly fail? A control that cannot fail is not a control.
"""
import numpy as np, os, sys, time
from scipy.stats import kstest
from multiprocessing import Pool
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "search"))
import wl_sweep as W

NSUR=int(os.environ.get("NSUR","300"))
MASK=os.environ.get("MASK","1")=="1"

def stat_nomask(v):
    d=np.diff(v); s=d.std()
    if s<=0: return np.nan
    from scipy.stats import kurtosis
    return float(kurtosis(d,fisher=True)-4.0*np.median(np.abs(d))/s)

def all_stats(Z,mask):
    f=W.stat if mask else stat_nomask
    return np.array([f(W.FORMS[fn](*[Z[i] for i in q])) for q,fn in W.TESTS])

def trial(seed):
    rng=np.random.default_rng(seed)
    Y=W.common_phase(W.X,rng)                 # NO four-way structure by construction
    obs=all_stats(Y,MASK)
    S=np.array([all_stats(W.common_phase(Y,rng),MASK) for _ in range(NSUR)])
    return np.array([(1.0+(S[:,t]>=obs[t]).sum())/(1.0+NSUR) for t in range(len(W.TESTS))])

if __name__=="__main__":
    NT=int(sys.argv[1]) if len(sys.argv)>1 else 200
    print("gap masking: %s   %d trials x %d surrogates x %d tests"
          %("ON" if MASK else "OFF (deliberate control)",NT,NSUR,len(W.TESTS)))
    t0=time.time()
    with Pool(80) as p: P=np.array(p.map(trial,range(5000,5000+NT),chunksize=1))
    print("%.1f s\n"%(time.time()-t0))
    print("  %-46s %-10s %8s %9s %9s"%("quadruple","form","mean p","KS p","frac<.05"))
    bad=0
    for t,(q,fn) in enumerate(W.TESTS):
        ps=P[:,t]; ps=ps[np.isfinite(ps)]
        if len(ps)<10: continue
        ks=kstest(ps,"uniform").pvalue
        if ks<0.05/len(W.TESTS): bad+=1
        print("  %-46s %-10s %8.4f %9.2e %9.4f%s"
              %(" / ".join(W.WL[i] for i in q),fn,ps.mean(),ks,(ps<0.05).mean(),
                "   <-- NOT UNIFORM" if ks<0.05/len(W.TESTS) else ""))
    allp=P[np.isfinite(P)]
    print("\n  POOLED: mean p %.4f (expect 0.5), frac<0.05 = %.4f (expect 0.05), KS p = %.2e"
          %(allp.mean(),(allp<0.05).mean(),kstest(allp,"uniform").pvalue))
    print("  tests failing uniformity at Bonferroni: %d of %d"%(bad,len(W.TESTS)))
