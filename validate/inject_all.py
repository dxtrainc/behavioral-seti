"""Injection curves for every periodogram-band row of section 4.2, on ONE
convention: the threshold the searches themselves use, mu = median(R)/ln2 with
alpha = 0.05, not an assumed Exp(1). Reports 90 / 95 / 99% recovery amplitudes.
"""
import numpy as np, os, json
from scipy.ndimage import median_filter
def _p(n):
    for c in (n, os.path.join("results",os.path.basename(n)),
              os.path.join(os.path.expanduser("~"),os.path.basename(n))):
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

def curve(x,band,label,quoted,ntrial=500,alpha=0.05,seed=5):
    r=prep(np.asarray(x,float)); n=len(r)//2*2; w=np.hanning(n); ws=w.sum()
    X=np.fft.rfft(r[:n]*w); f=np.fft.rfftfreq(n,d=60.0); P=np.abs(X)**2
    sel=(f>=band[0])&(f<=band[1]); ns=int(sel.sum())
    if ns<50: return dict(label=label,error="band has only %d bins"%ns)
    wid=min(2001, max(51, ns//4*2+1))
    C=np.exp(median_filter(np.log(np.maximum(P[sel],1e-300)),size=wid,mode="nearest"))
    Pn=P[sel]/C
    mu=float(np.median(Pn)/np.log(2.0))            # the searches' own scale
    thr=mu*np.log(ns/alpha)
    Xs=X[sel]; rng=np.random.default_rng(seed)
    lo,hi=quoted/6.0, quoted*30.0
    amps=np.exp(np.linspace(np.log(lo),np.log(hi),14))
    rec=[]
    for a in amps:
        add=a*ws/2.0; hit=0
        for _ in range(ntrial):
            k=int(rng.integers(3,ns-3)); ph=rng.uniform(0,2*np.pi)
            if abs(Xs[k]+add*np.exp(1j*ph))**2/C[k]>thr: hit+=1
        rec.append(hit/ntrial)
    def at(pct):
        for i in range(1,len(rec)):
            if rec[i]>=pct and rec[i-1]<pct:
                t=(pct-rec[i-1])/max(rec[i]-rec[i-1],1e-9)
                return float(np.exp(np.log(amps[i-1])+t*(np.log(amps[i])-np.log(amps[i-1]))))
        return float(amps[-1]) if rec[-1]<pct else float(amps[0])
    return dict(label=label,quoted=quoted,mu=mu,thr=float(thr),nbins=ns,
                amps=[float(v) for v in amps],rec=rec,
                a90=at(0.90),a95=at(0.95),a99=at(0.99))

E16=np.load(_p("euvs1m_g16.npz")); X16=np.load(_p("xrs1m_g16.npz"))
xk=[k for k in X16.files if "b" in k.lower()][0]
OU=np.load(_p("nm_OULU_60s.npz"))["v"]
JOBS=[("Lyman-alpha, 2 min",       E16["irr_1216"],(1/180.,1/90.),      4.0e-7),
      ("Lyman-alpha, 5 min",       E16["irr_1216"],(1/420.,1/240.),     2.7e-6),
      ("soft X-ray, 2 min",        X16[xk],        (1/180.,1/90.),      2.6e-4),
      ("soft X-ray, 1 day",        X16[xk],        (1/129600.,1/57600.),3.1e-2),
      ("cosmic ray, 5 min - 3 h",  OU,             (1/10800.,1/300.),   6.3e-5),
      ("cosmic ray, 15 min - 3 h", OU,             (1/10800.,1/900.),   4.2e-5)]
OUT=[]
print("%-26s %8s %6s %10s %10s %10s %10s"%("row","bins","mu","quoted","90%","95%","99%"))
for lab,arr,band,q in JOBS:
    c=curve(arr,band,lab,q); OUT.append(c)
    if "error" in c: print("%-26s %s"%(lab,c["error"])); continue
    print("%-26s %8d %6.3f %10.2e %10.2e %10.2e %10.2e"
          %(lab,c["nbins"],c["mu"],q,c["a90"],c["a95"],c["a99"]),flush=True)
json.dump(OUT,open(_p("injection_all.json") if os.path.isdir("results") else os.path.expanduser("~/injection_all.json"),"w"))
print("\nsaved injection_all.json")
