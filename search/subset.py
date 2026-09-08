"""Pulsar SUBSET search: does any small GROUP of millisecond pulsars have
back-extrapolated initial periods lying on a geometric ladder at a common epoch?

The population scan (backex2) needed ~40% of the sample on the ladder before it
fired, so it cannot address the hypothesis as actually proposed -- "several
long-lived pulsars in a group". Enumerating subsets is hopeless (C(137,5)=3.5e8
per epoch), so instead, for each epoch t and each candidate ratio rho, count how
many pulsars fall within eps of a rung of the ladder at the BEST phase. That is a
Hough transform over (rho, phase): O(N) per ratio rather than O(C(N,k)).

STATISTIC  M(t) = max over rho and phase of the rung occupancy.
NULL       the identical maximisation on surrogates with Pdot resampled from the
           local P-Pdot relation. The surrogate performs the SAME search, so the
           enormous look-elsewhere factor over (t, rho, phase) is paid by both
           sides and needs no analytic correction.
GATE       an injection carrying only k of N pulsars must be recovered, at the
           right epoch, before any null is believed.
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
from multiprocessing import Pool
yr=3.155693e7; Gyr=1e9*yr
d=json.load(open(_p("psrcat.json")))

EPS=0.002                      # tolerance in log10(P0): ~0.5% in period
NB=400                         # phase bins
# TWO ratio grids. The broad one is agnostic and pays an enormous
# look-elsewhere price. The restricted one follows sec 2.2: a sender who wants to
# be found must choose a key whose dictionary is cheap to exhaust, so the ladder
# ratio should be RECOGNISABLE -- a small-integer ratio or a familiar constant.
LR_BROAD=np.log10(np.linspace(1.02,2.00,220))
_R=[]
for q in range(1,9):
    for pn in range(q+1,3*q+1):
        r=pn/q
        if 1.05<=r<=3.0: _R.append(r)
_R += [2**0.5, 2**(1/3.), 3**0.5, (1+5**0.5)/2, np.pi/2, np.e/2]
LR_SMALL=np.log10(np.array(sorted(set(round(x,6) for x in _R))))
LR=LR_SMALL
if os.environ.get("ONE_RATIO"):
    LR=np.log10(np.array([float(os.environ["ONE_RATIO"])]))

def sel(pmax,pdmax=None,pmin=0.0,field=True,taumin=None):
    o=[(x["P"],x["Pd"]) for x in d
       if pmin<=x["P"]<pmax and x["Pd"] and x["Pd"]>0
       and (pdmax is None or x["Pd"]<pdmax)
       and (not field or "GC" not in (x["assoc"] or ""))
       and (taumin is None or x["P"]/(2.0*x["Pd"]) > taumin)]
    a=np.array(o); return a[:,0].copy(), a[:,1].copy()

def tmax_of(P,Pd):
    """a common epoch cannot exceed the YOUNGEST characteristic age present"""
    return float(np.min(P/(2.0*Pd)))

def best_at_ratio(L,lr,W=0.02,FLOOR=5e-5):
    """Excess rung occupancy in sigma, computed exactly.

    Binning was the previous failure: with an ABSOLUTE tolerance the window
    spans 2*EPS/lr of the phase circle, so the chance expectation collapses at
    large ratios and the maximisation is driven by Poisson spikes rather than by
    real ladders. Here the tolerance is a fixed FRACTION W of the rung spacing,
    so the expectation N*W is identical for every ratio and the statistic is
    comparable across the grid. The count is an exact sliding window on the
    sorted phases -- no histogram, no bin-edge artefacts."""
    if W*lr < FLOOR: return -1e9          # absolute tolerance below measurement precision
    n=len(L)
    frac=np.sort(np.mod(L/lr,1.0))
    ext=np.concatenate([frac,frac+1.0])
    j=np.searchsorted(ext,frac+W,side="right")
    cnt=int((j-np.arange(n)).max())
    exp=n*W
    return float((cnt-exp)/np.sqrt(exp))

def M_of_t(P,Pd,t):
    v=P*P-2.0*P*Pd*t
    m=v>0
    if m.sum()<20: return -1e9
    L=np.log10(np.sqrt(v[m]))
    return max(best_at_ratio(L,lr) for lr in LR)

def scan(P,Pd,tg): return np.array([M_of_t(P,Pd,t) for t in tg],dtype=float)

def surrogate(P,Pd,rng):
    lp=np.log10(P); lpd=np.log10(Pd)
    c=np.polyfit(lp,lpd,1); pred=np.polyval(c,lp)
    return 10**(pred+rng.permutation(lpd-pred))

_G={}
def _one(seed):
    rng=np.random.default_rng(seed)
    return float(scan(_G["P"],surrogate(_G["P"],_G["Pd"],rng),_G["tg"]).max())

def run(lbl,P,Pd,tg,ns,pool):
    S=scan(P,Pd,tg); k=int(np.argmax(S)); obs=float(S[k])
    _G.update(P=P,Pd=Pd,tg=tg)
    with Pool(pool) as pl: null=np.array(pl.map(_one,range(9000,9000+ns)))
    p=(1.0+(null>=obs).sum())/(1.0+len(null))
    z=(obs-null.mean())/max(null.std(),1e-12)
    print("   %-30s N=%4d  best excess %6.2f sd at t=%5.2f Gyr   null %5.2f+-%4.2f  z=%+5.2f  p=%.4f"
          %(lbl,len(P),obs,tg[k]/Gyr,null.mean(),null.std(),z,p),flush=True)
    return obs,tg[k],null

def inject(P,Pd,t_true,k,ratio,rng):
    K=P*Pd; v=P*P-2.0*K*t_true; m=v>0
    P0=np.sqrt(v[m]); Km=K[m]
    idx=rng.choice(len(P0),size=k,replace=False)
    base=np.median(P0)
    step=np.round(np.log(P0[idx]/base)/np.log(ratio))
    P0[idx]=base*ratio**step
    Pn=np.sqrt(P0**2+2.0*Km*t_true)
    return Pn, Km/Pn

if __name__=="__main__":
    NT=int(sys.argv[1]); NS=int(sys.argv[2]); POOL=int(sys.argv[3])
    rng=np.random.default_rng(7717)
    TAU=float(sys.argv[4])*Gyr if len(sys.argv)>4 else 6.0*Gyr
    P,Pd=sel(0.03,None,taumin=TAU)
    tmax=tmax_of(P,Pd)
    tg=np.linspace(0.0,0.98*tmax,NT)
    print("subset ladder search")
    print("  sample: field MSPs with characteristic age > %.1f Gyr  -> N=%d"%(TAU/Gyr,len(P)))
    print("  epoch scan bounded by the youngest characteristic age: 0 - %.2f Gyr"%(tmax/Gyr))
    print("  %d ratios, eps=%.3f dex, %d epochs; N is constant across t"%(len(LR),EPS,NT),flush=True)
    INJ_RATIO=float(os.environ.get("INJ_RATIO","1.5"))
    t_true=0.5*tmax
    print("  ratio grid: %d recognisable ratios (small-integer + constants)"%len(LR),flush=True)
    print("\nINJECTION at t=%.2f Gyr, ladder ratio %.4f -- how large must the group be?"%(t_true/Gyr,INJ_RATIO),flush=True)
    for k in (40,30,24,18,14,10,6):
        Pi,Pdi=inject(P,Pd,t_true,k,INJ_RATIO,rng)
        S=scan(Pi,Pdi,tg); j=int(np.argmax(S))
        ok=abs(tg[j]-t_true)/Gyr<0.4
        print("   k=%2d  best excess %6.2f sd at t=%5.2f Gyr (true %.2f)   %s"
              %(k,S[j],tg[j]/Gyr,t_true/Gyr,"RECOVERED" if ok else "missed"),flush=True)
    print("\nREAL DATA:",flush=True)
    run("field MSPs, tau > %.0f Gyr"%(TAU/Gyr),P,Pd,tg,NS,POOL)
    Pn,Pdn=sel(10.0,None,pmin=0.1,taumin=TAU)
    if len(Pn)>30:
        tn=np.linspace(0.0,0.98*tmax_of(Pn,Pdn),NT)
        run("CONTROL normal pulsars",Pn,Pdn,tn,NS,POOL)
