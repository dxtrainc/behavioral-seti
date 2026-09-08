"""Does log-2nd-diff manufacture survivors on correlated red noise with no signal?

log(a) - 2log(b) + log(c) CANCELS common-mode variation when the three series are
aligned and does not cancel when two of them are circularly shifted. The statistic
could therefore differ between observed and null for reasons of structure alone.
This measures the false-positive rate directly, on synthetic channels built to have
the red spectra and common driver of the real ones and NO embedded signal.
"""
import sys, time, numpy as np
from multiprocessing import Pool
from scipy import stats
from scipy.ndimage import median_filter

NDAY   = int(sys.argv[1]) if len(sys.argv)>1 else 5000
NTRIP  = int(sys.argv[2]) if len(sys.argv)>2 else 240
NSH    = int(sys.argv[3]) if len(sys.argv)>3 else 10000
POOL   = int(sys.argv[4]) if len(sys.argv)>4 else 8
ALPHA  = 0.05/10996          # the real triple sweep's Bonferroni threshold
W      = 181

def runmed(y,w): return median_filter(y,size=w,mode="nearest")
def clean(x):
    d=np.diff(x); s=1.4826*np.median(np.abs(d-np.median(d)))
    if s<=0: return None
    d=np.clip(d,-3.0*s,3.0*s); y=np.concatenate([[x[0]],x[0]+np.cumsum(d)])
    sd=y.std(); return (y-y.mean())/sd if sd>0 else None
def stat(x):
    z=clean(x)
    if z is None: return np.nan
    d=np.diff(z); s=d.std()
    if s<=0: return np.nan
    return float(stats.kurtosis(d,fisher=True)-4.0*np.median(np.abs(d))/s)
def prep(x):
    r=x/np.median(x); d=r-runmed(r,W); sd=d.std()
    return r,(d/sd if sd>0 else d)
def rollp(p,s): return (np.roll(p[0],s),np.roll(p[1],s))
def l2d(pa,pb,pc):
    la=np.log(np.maximum(pa[0],1e-12)); lb=np.log(np.maximum(pb[0],1e-12)); lc=np.log(np.maximum(pc[0],1e-12))
    return la-2.0*lb+lc
def rtp(pa,pb,pc): return pa[1]*pb[1]*pc[1]

def red(rng,n,alpha=0.92):
    """AR(1) red noise, positive, log-normal-ish -- like a solar activity proxy."""
    e=rng.standard_normal(n); x=np.zeros(n)
    for i in range(1,n): x[i]=alpha*x[i-1]+e[i]
    return np.exp(0.35*(x-x.mean())/max(x.std(),1e-9))+0.05

def one(seed):
    rng=np.random.default_rng(seed)
    n=NDAY
    drv=red(rng,n,0.97)                       # shared solar-activity driver
    w=rng.uniform(0.4,0.9,size=3)             # how strongly each channel follows it
    ch=[drv**w[k]*red(rng,n,0.90)**(1-w[k]) for k in range(3)]
    pa,pb,pc=[prep(c) for c in ch]
    out={}
    for nm,fn in (("log-2nd-diff",l2d),("resid-triple",rtp)):
        obs=stat(fn(pa,pb,pc))
        if not np.isfinite(obs): continue
        null=[]
        for s1,s2 in zip(rng.integers(1,n,size=NSH),rng.integers(1,n,size=NSH)):
            v=stat(fn(pa,rollp(pb,int(s1)),rollp(pc,int(s2))))
            if np.isfinite(v): null.append(v)
        null=np.asarray(null); ne=int((null>=obs).sum())
        out[nm]=dict(p=(1.0+ne)/(1.0+len(null)),obs=obs,
                     null_med=float(np.median(null)),ne=ne,
                     z=float((obs-np.mean(null))/max(np.std(null),1e-12)))
    return out

if __name__=="__main__":
    print("synthetic triples: %d, %d days, %d shifts, alpha=%.3e, pool=%d"
          %(NTRIP,NDAY,NSH,ALPHA,POOL),flush=True)
    t0=time.time()
    with Pool(POOL) as p: R=[r for r in p.map(one,range(1000,1000+NTRIP)) if r]
    print("done in %.1f min\n"%((time.time()-t0)/60),flush=True)
    for nm in ("log-2nd-diff","resid-triple"):
        v=[r[nm] for r in R if nm in r]
        if not v: continue
        ps=np.array([x["p"] for x in v]); zs=np.array([x["z"] for x in v])
        floor=1.0/(NSH+1)
        print("%-14s n=%d" % (nm,len(v)))
        print("   p at the %.2e floor : %d (%.1f%%)"%(floor,(ps<=floor+1e-15).sum(),100*(ps<=floor+1e-15).mean()))
        print("   p < 0.05            : %d (%.1f%%)   expected 5.0%%"%((ps<0.05).sum(),100*(ps<0.05).mean()))
        print("   p < 1e-3            : %d (%.1f%%)   expected 0.1%%"%((ps<1e-3).sum(),100*(ps<1e-3).mean()))
        print("   median z of obs vs null : %+.3f"%np.median(zs))
        print()
    print("VERDICT: a sound form gives ~5%% below 0.05 and ~0%% at the floor.")
