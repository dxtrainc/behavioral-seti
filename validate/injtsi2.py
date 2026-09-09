"""TSI 30-day occultation box, on QUIET EPOCHS as the original search selected.

The first attempt ran the box filter over all 16,801 days and recovered nothing
at six times the quoted depth. The original restricts to quiet days -- 4,919 of
them -- chosen by sunspot number, and that selection is the whole sensitivity of
the search: solar activity, not instrument noise, sets the floor.
"""
import numpy as np, os, json
from scipy.ndimage import median_filter, uniform_filter1d
def _p(n):
    for c in (n, os.path.join(os.path.expanduser("~"),os.path.basename(n))):
        if os.path.exists(c): return c
    return n
d=np.loadtxt(_p("tsi.csv"),delimiter=",",skiprows=1)
jd_t,v=d[:,0]+2309100.5,d[:,1]
ok=np.isfinite(v)&(v>1000); jd_t,v=jd_t[ok],v[ok]
sn=np.genfromtxt(_p("sunspot.csv"),delimiter=";")
m=sn[:,5]>=0
jd_s=(sn[m,3]-2000.0)*365.25+2451545.0; ss=sn[m,5]
S=np.interp(jd_t,jd_s,ss)
for pct in (30,):
    cut=np.percentile(S,pct); sel=S<=cut
    r=v[sel]/np.median(v[sel]); r=r-median_filter(r,size=181,mode="nearest")
    n=len(r); W=30
    def boxstat(x):
        b=uniform_filter1d(x,size=W,mode="nearest")
        s=1.4826*np.median(np.abs(b-np.median(b)))
        return (np.median(b)-b)/max(s,1e-12)
    thr=float(np.quantile(boxstat(r),1.0-0.05/n))
    print("quiet %d%% of days (SN <= %.0f): %d days, threshold %.2f sigma"%(pct,cut,n,thr),flush=True)
    rng=np.random.default_rng(11)
    amps=np.array([2e-5,4e-5,6.2e-5,1e-4,1.5e-4,2.5e-4])
    rec=[]
    for a in amps:
        hit=0
        for _ in range(250):
            k=int(rng.integers(100,n-100-W))
            y=r.copy(); y[k:k+W]-=a
            if boxstat(y)[k+W//2]>thr: hit+=1
        rec.append(hit/250.); print("   depth %8.1e -> %5.1f%%"%(a,100*rec[-1]),flush=True)
    def at(p):
        for i in range(1,len(rec)):
            if rec[i]>=p and rec[i-1]<p:
                t=(p-rec[i-1])/max(rec[i]-rec[i-1],1e-9)
                return float(np.exp(np.log(amps[i-1])+t*(np.log(amps[i])-np.log(amps[i-1]))))
        return None
    a95=at(.95)
    print("   quoted 6.2e-5;  95%% recovery %s"%("%.1e"%a95 if a95 else "above the grid"),flush=True)
    json.dump(dict(amps=amps.tolist(),rec=rec,a95=a95,quoted=6.2e-5,ndays=n),
              open(os.path.expanduser("~/injection_tsi.json"),"w"))
