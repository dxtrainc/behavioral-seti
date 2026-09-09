"""Injection for two of the five non-periodogram rows.

  TSI 30-day occultation box : matched box filter on the daily TSI composite
  sidereal fold              : fixed-frequency amplitude at 366.2422 cyc/yr,
                               with the anti-sidereal line as the noise scale

Both detectors are unambiguous enough to implement from the definition. The
other three rows (achromatic, multi-scale, self-keyed) have their own scripts
and should be injected by adding a flag to those, not by re-implementation --
re-implementing a detector is what produced the 127 uHz mode-comb error.
"""
import numpy as np, os, json
from scipy.ndimage import median_filter, uniform_filter1d
def _p(n):
    for c in (n, os.path.join(os.path.expanduser("~"),os.path.basename(n))): 
        if os.path.exists(c): return c
    return n

# ---------------- TSI: 30-day box dip -------------------------------------
d=np.loadtxt(_p("tsi.csv"),delimiter=",",skiprows=1)
jd,v=d[:,0],d[:,1]
ok=np.isfinite(v)&(v>1000)
v=v[ok]; n=len(v)
r=v/np.median(v)
r=r-median_filter(r,size=181,mode="nearest")
W=30
def boxstat(x):
    b=uniform_filter1d(x,size=W,mode="nearest")
    s=1.4826*np.median(np.abs(b-np.median(b)))
    return (np.median(b)-b)/max(s,1e-12)          # depth in robust sigma
base=boxstat(r); thr=float(np.quantile(base,1.0-0.05/len(base)))
print("TSI: %d daily points, box %d d, threshold %.2f sigma (FWER 0.05)"%(n,W,thr))
rng=np.random.default_rng(4)
amps=np.array([1e-5,2e-5,4e-5,6.2e-5,1e-4,2e-4,4e-4])
rec=[]
for a in amps:
    hit=0
    for _ in range(300):
        k=int(rng.integers(200,n-200-W))
        y=r.copy(); y[k:k+W]-=a
        if boxstat(y)[k+W//2]>thr: hit+=1
    rec.append(hit/300.); print("   depth %8.1e -> %5.1f%%"%(a,100*rec[-1]),flush=True)
def at(p,A,R):
    for i in range(1,len(R)):
        if R[i]>=p and R[i-1]<p:
            t=(p-R[i-1])/max(R[i]-R[i-1],1e-9)
            return float(np.exp(np.log(A[i-1])+t*(np.log(A[i])-np.log(A[i-1]))))
    return None
tsi95=at(.95,amps,rec)
print("   quoted 62 ppm (6.2e-5);  95%% recovery %s\n"%("%.1e"%tsi95 if tsi95 else ">4e-4"))

# ---------------- sidereal fold -------------------------------------------
z=np.load(_p("nm_OULU_60s.npz")); w=z["v"].astype(float)
m=np.isfinite(w); w=np.where(m,w,np.nanmedian(w[m]))
w=w/np.median(w)-1.0
w=w-median_filter(w,size=1441,mode="nearest")
step=5; t=np.arange(0,len(w),step)*60.0; x=w[::step]
yr=365.2422*86400.0
def amp(f):
    ph=2*np.pi*f*t
    return 2.0*np.hypot(np.dot(x,np.cos(ph)),np.dot(x,np.sin(ph)))/len(x)
noise=np.array([amp(c/yr) for c in np.linspace(358,363,26)])
sd=float(np.std(noise)); thr_s=5.0*sd
print("sidereal: %d samples used, noise sd %.3e, threshold 5 sd = %.3e"%(len(x),sd,thr_s))
fs=366.2422/yr
rec2=[]; amps2=np.array([5e-5,1e-4,1.8e-4,3e-4,6e-4,1e-3])
for a in amps2:
    hit=0
    for _ in range(60):
        ph0=rng.uniform(0,2*np.pi)
        xi=x+a*np.cos(2*np.pi*fs*t+ph0)
        ph=2*np.pi*fs*t
        A=2.0*np.hypot(np.dot(xi,np.cos(ph)),np.dot(xi,np.sin(ph)))/len(xi)
        if A>thr_s: hit+=1
    rec2.append(hit/60.); print("   amp %8.1e -> %5.1f%%"%(a,100*rec2[-1]),flush=True)
sid95=at(.95,amps2,rec2)
print("   quoted 1.8e-4;  95%% recovery %s"%("%.1e"%sid95 if sid95 else "below the grid"))
json.dump(dict(tsi=dict(amps=amps.tolist(),rec=rec,a95=tsi95,quoted=6.2e-5),
               sidereal=dict(amps=amps2.tolist(),rec=rec2,a95=sid95,quoted=1.8e-4)),
          open(os.path.expanduser("~/injection_two.json"),"w"))
print("\nsaved injection_two.json")
