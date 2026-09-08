"""Which threshold did the searches actually use, and what recovery does it give?

The original sets mu = median(R)/ln2 -- an EMPIRICAL scale read off the observed
power distribution -- and alpha = 0.05. The injection test used mu = 1 (Exp(1))
and alpha = 0.01. If mu_emp differs from 1, my claim that "the analytic threshold
is optimistic" mischaracterises the method: the empirical scale already absorbs
the non-exponentiality.
"""
import numpy as np, os
from scipy.ndimage import median_filter
def _p(n):
    for c in (n, os.path.join(os.path.expanduser("~"),n)):
        if os.path.exists(c): return c
    return n
def prep(x,wmin=60):
    ok=np.isfinite(x)&(x>0); y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med)); y=np.clip(y,med-5*s,med+5*s)
    a=np.nan_to_num(y); m=ok.astype(float)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    trend=np.where(cm[hi]-cm[lo]>5,(ca[hi]-ca[lo])/np.maximum(cm[hi]-cm[lo],1),med)
    return np.where(ok,y-trend,0.0)

for tag,f,col,band,quoted in (("Ly-alpha 2min","euvs1m_g16.npz","irr_1216",(1/180.,1/90.),4.0e-7),
                              ("X-ray 2min","xrs1m_g16.npz",None,(1/180.,1/90.),2.6e-4),
                              ("X-ray 1 day","xrs1m_g16.npz",None,(1/129600.,1/57600.),3.1e-2)):
    z=np.load(_p(f)); k=col or [c for c in z.files if "b" in c.lower()][0]
    r=prep(np.asarray(z[k],float)); n=len(r)//2*2; w=np.hanning(n)
    X=np.fft.rfft(r[:n]*w); fr=np.fft.rfftfreq(n,d=60.0); P=np.abs(X)**2
    sel=(fr>=band[0])&(fr<=band[1]); ns=int(sel.sum())
    if ns<20: print("%-14s band too narrow (%d bins)"%(tag,ns)); continue
    C=np.exp(median_filter(np.log(np.maximum(P[sel],1e-300)),size=min(2001,ns//4*2+1),mode="nearest"))
    Pn=P[sel]/C
    mu=float(np.median(Pn)/np.log(2.0))
    thr_orig=mu*np.log(ns/0.05); thr_mine=np.log(ns/0.01)
    eps=lambda t: 2.0*np.sqrt(t*np.median(C))/w.sum()
    print("%-14s bins %8d  mu_emp %.3f  thr_orig %6.2f  thr_mine %6.2f"%(tag,ns,mu,thr_orig,thr_mine))
    print("   threshold-crossing eps:  original convention %.3e   quoted %.1e   ratio %.2f"
          %(eps(thr_orig),quoted,eps(thr_orig)/quoted))
    # recovery at the quoted amplitude, under BOTH thresholds
    rng=np.random.default_rng(3); Xs=X[sel]
    for lbl,t in (("original",thr_orig),("mine",thr_mine)):
        hit=0
        for _ in range(400):
            i=int(rng.integers(5,ns-5)); ph=rng.uniform(0,2*np.pi)
            if abs(Xs[i]+quoted*w.sum()/2*np.exp(1j*ph))**2/C[i]>t: hit+=1
        print("     recovery at the quoted %.1e under the %-8s threshold: %5.1f%%"%(quoted,lbl,100*hit/400))
