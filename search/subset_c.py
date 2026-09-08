"""Subset ladder search on Shklovskii+Galactic-CORRECTED Pdot."""
import json, os, sys, math, numpy as np
from multiprocessing import Pool
yr=3.155693e7; Gyr=1e9*yr
R=json.load(open(os.path.expanduser("~/psrcat_shk.json")))
NB=400; EPS=0.002
_r=[]
for q in range(1,9):
    for pn in range(q+1,3*q+1):
        v=pn/q
        if 1.05<=v<=3.0: _r.append(v)
_r+=[2**0.5,2**(1/3.),3**0.5,(1+5**0.5)/2,math.pi/2,math.e/2]
LR=np.log10(np.array(sorted(set(round(v,6) for v in _r))))
if os.environ.get("ONE_RATIO"): LR=np.log10(np.array([float(os.environ["ONE_RATIO"])]))

def sel(corrected, taumin):
    o=[]
    for x in R:
        if not (x["P"]<0.03): continue
        pd=x["intr"] if corrected else x["Pd"]
        if pd is None or pd<=0: continue
        if x["P"]/(2.0*pd) <= taumin: continue
        o.append((x["P"],pd))
    a=np.array(o); return a[:,0].copy(), a[:,1].copy()

def best_at_ratio(L,lr,W=0.02,FLOOR=5e-5):
    if W*lr<FLOOR: return -1e9
    n=len(L); frac=np.sort(np.mod(L/lr,1.0))
    ext=np.concatenate([frac,frac+1.0])
    j=np.searchsorted(ext,frac+W,side="right")
    cnt=int((j-np.arange(n)).max()); exp=n*W
    return float((cnt-exp)/np.sqrt(exp))
def M_of_t(P,Pd,t):
    v=P*P-2.0*P*Pd*t; m=v>0
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
def inject(P,Pd,t_true,k,ratio,rng):
    K=P*Pd; v=P*P-2.0*K*t_true; m=v>0
    P0=np.sqrt(v[m]); Km=K[m]
    idx=rng.choice(len(P0),size=min(k,len(P0)),replace=False)
    base=np.median(P0); step=np.round(np.log(P0[idx]/base)/np.log(ratio))
    P0[idx]=base*ratio**step
    Pn=np.sqrt(P0**2+2.0*Km*t_true)
    return Pn, Km/Pn

if __name__=="__main__":
    NT=int(sys.argv[1]); NS=int(sys.argv[2]); POOL=int(sys.argv[3]); TAU=float(sys.argv[4])*Gyr
    rng=np.random.default_rng(4242)
    for corrected in (False,True):
        P,Pd=sel(corrected,TAU)
        tmax=float(np.min(P/(2.0*Pd))); tg=np.linspace(0.0,0.98*tmax,NT)
        lbl="CORRECTED Pdot" if corrected else "observed Pdot"
        print("\n=== %s ===  N=%d, epoch scan 0-%.2f Gyr"%(lbl,len(P),tmax/Gyr),flush=True)
        t_true=0.5*tmax
        for k in (20,14,10,6):
            Pi,Pdi=inject(P,Pd,t_true,k,1.5,rng)
            S=scan(Pi,Pdi,tg); j=int(np.argmax(S))
            print("   inject k=%2d -> %6.2f sd at %5.2f Gyr (true %.2f) %s"
                  %(k,S[j],tg[j]/Gyr,t_true/Gyr,"OK" if abs(tg[j]-t_true)/Gyr<0.4 else "missed"),flush=True)
        S=scan(P,Pd,tg); j=int(np.argmax(S)); obs=float(S[j])
        _G.update(P=P,Pd=Pd,tg=tg)
        with Pool(POOL) as pl: null=np.array(pl.map(_one,range(500,500+NS)))
        z=(obs-null.mean())/max(null.std(),1e-12)
        p=(1.0+(null>=obs).sum())/(1.0+len(null))
        print("   REAL: best %.2f sd at %5.2f Gyr   null %.2f+-%.2f   z=%+.2f  p=%.4f"
              %(obs,tg[j]/Gyr,null.mean(),null.std(),z,p),flush=True)
