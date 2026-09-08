"""The 1-day X-ray row, on the DAILY series the original search used.

Injecting into the 1-minute series in a 1-day band saturated at the lowest
amplitude tried, which suggests the quoted 3.1e-2 came from the daily background
product (13,289 points) rather than from 1-minute data. This redoes it there.
"""
import numpy as np, os
from scipy.ndimage import median_filter
def _p(n):
    for c in (n, os.path.join(os.path.expanduser("~"),os.path.basename(n))):
        if os.path.exists(c): return c
    return n
rows=[l.split(",") for l in open(_p("goes_xray_bg.csv")).read().strip().split("\n")[1:]]
v=np.array([float(r[1]) for r in rows]); ok=np.isfinite(v)&(v>1e-10)&(v<1e-4)
y=np.log(v[ok]); print("daily X-ray background: %d points"%len(y))
y=y-median_filter(y,size=181,mode="nearest")
n=len(y)//2*2; w=np.hanning(n); ws=w.sum()
X=np.fft.rfft(y[:n]*w); f=np.fft.rfftfreq(n,d=86400.0); P=np.abs(X)**2
sel=f>0
C=np.exp(median_filter(np.log(np.maximum(P[sel],1e-300)),size=201,mode="nearest"))
Pn=P[sel]/C; ns=int(sel.sum())
mu=float(np.median(Pn)/np.log(2.0)); thr=mu*np.log(ns/0.05)
print("bins %d  mu %.3f  thr %.2f"%(ns,mu,thr))
Xs=X[sel]; rng=np.random.default_rng(9)
amps=np.exp(np.linspace(np.log(1e-3),np.log(2e-1),16))
rec=[]
for a in amps:
    add=a*ws/2.0; hit=0
    for _ in range(600):
        k=int(rng.integers(2,ns-2)); ph=rng.uniform(0,2*np.pi)
        if abs(Xs[k]+add*np.exp(1j*ph))**2/C[k]>thr: hit+=1
    rec.append(hit/600.)
    print("   %9.3e -> %5.1f%%"%(a,100*rec[-1]),flush=True)
def at(p):
    for i in range(1,len(rec)):
        if rec[i]>=p and rec[i-1]<p:
            t=(p-rec[i-1])/max(rec[i]-rec[i-1],1e-9)
            return float(np.exp(np.log(amps[i-1])+t*(np.log(amps[i])-np.log(amps[i-1]))))
    return None
print("\nquoted 3.10e-02   90%% %.2e   95%% %.2e   99%% %.2e"%(at(.90) or -1,at(.95) or -1,at(.99) or -1))
a99=at(.99)
if a99: print("cadence ratio, injection-verified 99%%: 1 day / 2 min = %.2e / 6.87e-04 = %.0f"%(a99,a99/6.87e-4))
