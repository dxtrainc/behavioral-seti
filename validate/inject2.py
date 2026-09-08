"""Injection-recovery curves, v2.

v1 had two errors: the injected sinusoid was given A*sum(w)/4 in its bin instead
of A*sum(w)/2, understating every amplitude by a factor of two; and the detection
threshold came from a sample quantile, which reports the largest value present
rather than a false-alarm rate. Both are fixed. The threshold is now analytic --
for continuum-normalised power that is exponential under the null, a family-wise
false-alarm probability alpha over n_bins needs z > ln(n_bins/alpha) -- and the
exponential assumption is CHECKED against the data rather than assumed.
"""
import numpy as np, os, json
from scipy.ndimage import median_filter
D=os.path.expanduser("~")

def prep(x, wmin=60):
    ok=np.isfinite(x)&(x>0)
    y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-5*s,med+5*s)
    a=np.nan_to_num(y); m=ok.astype(float)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    num=ca[hi]-ca[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>5,num/np.maximum(den,1),med)
    return np.where(ok,y-trend,0.0)

def run(x, dt, band, name, amps, ntrial=600, alpha=0.01, seed=7):
    r=prep(np.asarray(x,float)); n=len(r)//2*2
    w=np.hanning(n); rw=r[:n]*w; wsum=w.sum()
    X=np.fft.rfft(rw); f=np.fft.rfftfreq(n,d=dt); P=np.abs(X)**2
    sel=(f>=band[0])&(f<=band[1]); ns=int(sel.sum())
    C=np.exp(median_filter(np.log(np.maximum(P[sel],1e-300)),size=2001,mode="nearest"))
    Pn=P[sel]/C
    thr=np.log(ns/alpha)                       # analytic FWER threshold, Exp(1) null
    # is the null actually exponential? mean should be ~1 and the tail should match
    mean=float(Pn.mean()); frac=float((Pn>thr).mean()); expect=float(np.exp(-thr))
    rng=np.random.default_rng(seed)
    Cs=C; out=[]
    for a in amps:
        add=a*wsum/2.0                          # correct on-bin amplitude
        hit=0
        for _ in range(ntrial):
            k=int(rng.integers(5,ns-5)); ph=rng.uniform(0,2*np.pi)
            Xk=X[sel][k]+add*np.exp(1j*ph)
            if (abs(Xk)**2)/Cs[k]>thr: hit+=1
        out.append(hit/ntrial)
    return dict(name=name,amps=[float(v) for v in amps],rec=out,thr=float(thr),
                nbins=ns,null_mean=mean,tail_obs=frac,tail_exp=expect)

JOBS=[]
z=np.load(D+"/euvs1m_g16.npz")
JOBS.append(("Lyman-alpha, 2 min",z["irr_1216"],60.0,(1/180.,1/90.),
             [1e-7,2e-7,4e-7,6e-7,1e-6,2e-6]))
z2=np.load(D+"/xrs1m_g16.npz"); k2=[k for k in z2.files if "b" in k.lower()] or list(z2.files)
JOBS.append(("soft X-ray, 2 min",z2[k2[0]],60.0,(1/180.,1/90.),
             [5e-5,1e-4,2e-4,2.6e-4,4e-4,8e-4]))
z3=np.load(D+"/nm_OULU_60s.npz")
JOBS.append(("cosmic ray, 5 min - 3 h",z3["v"],60.0,(1/10800.,1/300.),
             [1e-5,3e-5,6.3e-5,1e-4,2e-4,4e-4]))

OUT=[]
for name,arr,dt,band,amps in JOBS:
    c=run(arr,dt,band,name,amps); OUT.append(c)
    print("%-24s bins %8d  thr %.1f   null mean %.3f (want 1.0)   tail obs %.2e vs exp %.2e"
          %(name,c["nbins"],c["thr"],c["null_mean"],c["tail_obs"],c["tail_exp"]),flush=True)
    for a,r in zip(c["amps"],c["rec"]): print("      %9.2e -> %5.1f%%"%(a,100*r),flush=True)
json.dump(OUT,open(D+"/injection.json","w"))
print("saved ~/injection.json")
