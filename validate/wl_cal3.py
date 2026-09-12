"""GAP CONTROL THAT CAN ACTUALLY FAIL.

The previous control was a no-op: its synthetic data was phase-randomised on the
COMPRESSED index axis, so it was smooth in index space and had no gap-jumps. Turning
the masking off cannot hurt a series with nothing to mask.

The real asymmetry is this. Observed data lives on a CALENDAR with 1,064 breaks; we
compress it to an index axis, which silently places a 38-day drift next to a 1-day
drift. Surrogates are generated on the compressed axis and are smooth there. So the
observed statistic sees jumps its own null never does.

To reproduce that, synthetic null data must be built on the FULL grid and THEN
subsampled to the overlap days -- giving genuine gap-jumps with no four-way
structure. If the masking matters, the unmasked arm must now over-fire.
"""
import numpy as np, os, sys, time
from scipy.stats import kstest, norm
from multiprocessing import Pool
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "search"))
import wl_sweep as W

FULL=np.arange(len(W.OV))
def fill(ch):
    v=W.M[W.ix[ch]].astype(float); m=np.isfinite(v)
    return np.interp(FULL,FULL[m],v[m])
XF=np.vstack([norm.ppf((np.argsort(np.argsort(fill(c)))+0.5)/len(FULL)) for c in W.WL])
NF=XF.shape[1]
SUB=np.where(W.OV)[0]

def full_phase(rng):
    F=np.fft.rfft(XF,axis=1)
    ph=rng.uniform(0,2*np.pi,F.shape[1]); ph[0]=0.0
    if NF%2==0: ph[-1]=0.0
    return np.fft.irfft(F*np.exp(1j*ph)[None,:],n=NF,axis=1)

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
    Y=full_phase(rng)[:,SUB]                  # gap-jumps present, no 4-way structure
    obs=all_stats(Y,MASK)
    S=np.array([all_stats(W.common_phase(Y,rng),MASK) for _ in range(NSUR)])
    return np.array([(1.0+(S[:,t]>=obs[t]).sum())/(1.0+NSUR) for t in range(len(W.TESTS))])

if __name__=="__main__":
    NT=int(sys.argv[1]) if len(sys.argv)>1 else 200
    print("GAP-REALISTIC control. masking: %s  %d trials x %d sur x %d tests"
          %("ON" if MASK else "OFF",NT,NSUR,len(W.TESTS)))
    print("increment size at breaks vs non-breaks, real data:")
    for i,c in enumerate(W.WL):
        d=np.abs(np.diff(W.X[i]))
        print("   %-12s break %.4f   non-break %.4f   ratio %.2f"
              %(c,d[~W.KEEP].mean(),d[W.KEEP].mean(),d[~W.KEEP].mean()/d[W.KEEP].mean()))
    t0=time.time()
    with Pool(80) as p: P=np.array(p.map(trial,range(7000,7000+NT),chunksize=1))
    print("\n%.1f s"%(time.time()-t0))
    bad=sum(1 for t in range(len(W.TESTS))
            if kstest(P[:,t][np.isfinite(P[:,t])],"uniform").pvalue<0.05/len(W.TESTS))
    allp=P[np.isfinite(P)]
    print("  POOLED mean p %.4f (expect 0.5)   frac<0.05 = %.4f (expect 0.05)   ratio %.2fx"
          %(allp.mean(),(allp<0.05).mean(),(allp<0.05).mean()/0.05))
    print("  tests failing uniformity at Bonferroni: %d of %d"%(bad,len(W.TESTS)))
