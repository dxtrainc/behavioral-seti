"""Back-extrapolated initial-period search, v2 -- with a detector that can fail.

v1's statistic sat at 0.954 for data AND surrogates: autocorrelating the pairwise
log-ratio histogram measured its smooth envelope, not comb structure. Same failure
as the p-mode search against a global median. Fixed by dividing out a LOCAL
continuum before autocorrelating, and gated on an injection that must be recovered.
"""
import os as _os
def _p(name):
    """Resolve a data or result path: as given, then ./results, then $HOME.
    These scripts were written to run from a home directory; this lets the
    repository be cloned anywhere without editing them."""
    for c in (name, _os.path.join("results", _os.path.basename(name)),
              _os.path.join(_os.path.dirname(__file__), "..", "results", _os.path.basename(name)),
              _os.path.join(_os.path.expanduser("~"), _os.path.basename(name))):
        if _os.path.exists(c): return c
    return _os.path.expanduser(name)

import json, os, sys, numpy as np
from scipy.ndimage import median_filter
yr=3.155693e7; Gyr=1e9*yr
d=json.load(open(_p("psrcat.json")))
NB=200; RANGE=(0.0,1.2); CONT=25; LAGLO,LAGHI=3,80

def sel(pmax,pdmax=None,pmin=0.0,field=True):
    o=[(x["P"],x["Pd"]) for x in d
       if pmin<=x["P"]<pmax and x["Pd"] and x["Pd"]>0
       and (pdmax is None or x["Pd"]<pdmax)
       and (not field or "GC" not in (x["assoc"] or ""))]
    a=np.array(o); return a[:,0].copy(), a[:,1].copy()

def stat(P0):
    if len(P0)<25: return np.nan
    L=np.log10(P0); n=len(L)
    D=np.abs(L[:,None]-L[None,:])[np.triu_indices(n,1)]
    h,_=np.histogram(D,bins=NB,range=RANGE); h=h.astype(float)
    cont=median_filter(h,size=CONT,mode="nearest")     # LOCAL continuum
    r=h/np.maximum(cont,0.5)-1.0
    if r.std()<=0: return np.nan
    ac=np.correlate(r,r,"full")[len(r)-1:]
    if ac[0]<=0: return np.nan
    return float(np.max(ac[LAGLO:LAGHI]/ac[0]))

def scan(P,Pd,tg):
    o=np.full(len(tg),np.nan)
    for i,t in enumerate(tg):
        v=P*P-2.0*P*Pd*t; m=v>0
        if m.sum()>=25: o[i]=stat(np.sqrt(v[m]))
    return o

def surrogate(P,Pd,rng):
    lp=np.log10(P); lpd=np.log10(Pd)
    c=np.polyfit(lp,lpd,1); pred=np.polyval(c,lp)
    return 10**(pred+rng.permutation(lpd-pred))

def inject(P,Pd,t_true,frac,ratio,rng):
    """snap a fraction of the TRUE initial periods onto a geometric ladder,
    then evolve them forward so 'today' is self-consistent."""
    K=P*Pd                                   # under n=3, K = P*Pdot is CONSTANT
    v=P*P-2.0*K*t_true; m=v>0
    P0=np.sqrt(v[m]); Km=K[m]
    k=rng.choice(len(P0),size=max(3,int(frac*len(P0))),replace=False)
    base=np.median(P0)
    step=np.round(np.log(P0[k]/base)/np.log(ratio))
    P0[k]=base*ratio**step
    Pnew=np.sqrt(P0**2+2.0*Km*t_true)        # same K -- self-consistent
    return Pnew, Km/Pnew                     # Pdot implied by the constant K

if __name__=="__main__":
    NT=int(sys.argv[1]); NS=int(sys.argv[2])
    tg=np.linspace(0.0,13.8,NT)*Gyr
    rng=np.random.default_rng(20260908)
    P,Pd=sel(0.03,1e-20)
    print("detector check on %d field MSPs (Pdot<1e-20)"%len(P),flush=True)
    print("  baseline statistic on the real set at t=0: %.4f"%stat(P),flush=True)

    print("\nINJECTION -- the detector must recover a planted ladder:",flush=True)
    t_true=5.0*Gyr
    for frac,ratio in ((1.0,1.15),(0.7,1.15),(0.5,1.15),(0.3,1.15),(0.7,1.30),(1.0,1.05)):
        Pi,Pdi=inject(P,Pd,t_true,frac,ratio,rng)
        S=scan(Pi,Pdi,tg); k=int(np.nanargmax(S))
        null=np.array([np.nanmax(scan(Pi,surrogate(Pi,Pdi,rng),tg)) for _ in range(40)])
        z=(S[k]-null.mean())/max(null.std(),1e-12)
        print("   frac=%3.0f%% ratio=%4.2f -> peak %.4f at t=%4.2f Gyr (true 5.00)  z=%+6.2f  %s"
              %(100*frac,ratio,S[k],tg[k]/Gyr,z,"RECOVERED" if z>4 and abs(tg[k]/Gyr-5.0)<1.5 else "missed"),flush=True)

    print("\nREAL DATA:",flush=True)
    for lbl,a in (("field MSPs Pdot<1e-20",dict(pmax=0.03,pdmax=1e-20)),
                  ("field MSPs all",       dict(pmax=0.03)),
                  ("CONTROL normal pulsars",dict(pmax=10.0,pmin=0.1))):
        Px,Pdx=sel(**a)
        S=scan(Px,Pdx,tg)
        if np.all(np.isnan(S)): print("   %-24s no usable t"%lbl); continue
        k=int(np.nanargmax(S))
        null=np.array([np.nanmax(scan(Px,surrogate(Px,Pdx,rng),tg)) for _ in range(NS)])
        p=(1.0+(null>=S[k]).sum())/(1.0+len(null))
        z=(S[k]-null.mean())/max(null.std(),1e-12)
        print("   %-24s N=%4d peak %.4f at t=%5.2f Gyr  null %.4f+-%.4f  z=%+5.2f  p=%.4f"
              %(lbl,len(Px),S[k],tg[k]/Gyr,null.mean(),null.std(),z,p),flush=True)
