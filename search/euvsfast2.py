"""
THE FAST SEARCH ON EUVS, VERSION 2 -- after version 1 found the spacecraft's day.

Version 1 returned 145 distinct candidates, 50 of them in more than one line, and
they PASSED the GOES-17 twin gate. They are still not astrophysics. Every one is an
exact integer multiple of 1/86400 Hz:

    6.1111111e-3 Hz = 528/day = 22/hour      6.0995e-3 = 527/day
    8.0555556e-3 Hz = 696/day = 29/hour      8.0671e-3 = 697/day

A comb at harmonics 500-700 of one cycle per day is the Fourier series of a SHARP
FEATURE THAT RECURS AT THE SAME TIME EVERY DAY -- an eclipse entry, a calibration
event, a recurring gap. Nothing else produces a comb that far up.

TWO GATES FAILED HERE, AND THE REASON IS THE SAME:

  gate 2, the twin.  GOES-16 and GOES-17 carry the same INSTRUMENT DESIGN and fly
  the same geostationary orbit with the same calibration schedule. A twin separates
  a per-UNIT artifact from the sky. It cannot separate a per-DESIGN artifact, which
  appears identically in both. The Oulu/Kiel gate worked because the stations are
  different hardware built decades apart; this one does not.

  gate 3, multi-line.  25.6, 28.4 and 30.4 nm are all EUVS-A, one sensor. A sensor
  artifact appears in all three "independent" lines at once. Those three channels
  returned 27,000-35,000 bins above threshold against 10-33 for the others, which is
  the tell.

THE FIX removes the whole comb in one operation rather than notching 700 harmonics:
FOLD AT ONE DAY, SUBTRACT THE MEAN DAILY PROFILE. Every k/86400 component dies at
once, including the ones too weak to have been listed. What survives is anything not
locked to the spacecraft's day -- which is the only thing that could be solar.
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


import os, sys
import numpy as np
from scipy.ndimage import median_filter

D=os.path.expanduser("~")
CH=["irr_256","irr_284","irr_304","irr_1175","irr_1216","irr_1335","irr_1405","MgII_EXIS"]
NICE={"irr_256":"25.6nm","irr_284":"28.4nm","irr_304":"30.4nm","irr_1175":"117.5nm",
      "irr_1216":"Ly-alpha","irr_1335":"133.5nm","irr_1405":"140.5nm","MgII_EXIS":"MgII"}

def load(tok):
    z=np.load(_p("euvs1m_%s.npz")%tok)
    return {c:z[c] for c in CH}

def daily_profile(r, ok, day=1440):
    """subtract the mean profile over one day. Kills every harmonic of 1/day at
    once -- the comb, its leakage, and the components below the listing threshold."""
    n=(len(r)//day)*day
    v=r[:n].reshape(-1,day); m=ok[:n].reshape(-1,day)
    prof=np.where(m.sum(axis=0)>10, (v*m).sum(axis=0)/np.maximum(m.sum(axis=0),1), 0.0)
    out=r.copy(); out[:n]=(v-prof).ravel(); out[~ok]=0.0
    return out, prof

def prep(x, wmin=1440):
    ok=np.isfinite(x)&(x>0)
    y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-5*s,med+5*s)
    a=np.nan_to_num(y); m=ok.astype(np.float64)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    num=ca[hi]-ca[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>10,num/np.maximum(den,1),med)
    return np.where(ok,y-trend,0.0), ok

def notch(r, ok, periods, step=60.0):
    t=np.arange(len(r),dtype=np.float64)*step
    A=np.vstack(sum(([np.cos(2*np.pi*t/P),np.sin(2*np.pi*t/P)] for P in periods),[])).T
    coef,*_=np.linalg.lstsq(A[ok],r[ok],rcond=None)
    out=r-A@coef; out[~ok]=0.0
    return out, coef

def spec(r,step=60.0):
    n=len(r)//2*2; w=np.blackman(n)
    P=np.abs(np.fft.rfft(r[:n]*w))**2
    f=np.fft.rfftfreq(n,d=step)
    return f[1:],P[1:],w
def cont(P,w=801):
    return np.exp(median_filter(np.log(np.maximum(P,1e-300)),size=w,mode="nearest"))

LINES=[86400.0/k for k in (1,2,3,4,5,6)]
G16=load("g16")
print("EUVS 1-MINUTE FAST SEARCH v2 -- mean daily profile removed\n")
res={}
for c in CH:
    r,ok=prep(G16[c])
    r,prof=daily_profile(r,ok)
    if ok.mean()<0.3:
        print("  %-10s only %.1f%% present -- skipped"%(NICE[c],100*ok.mean())); continue
    rn,coef=notch(r,ok,LINES)
    f,P,w=spec(rn); C=cont(P); R=P/C
    N=len(R); thr=(np.median(R)/np.log(2.0))*np.log(N/0.05)
    idx=np.where(R>thr)[0]
    res[c]=(f,R,thr,C,w,idx)
    d24=np.hypot(coef[0],coef[1])
    print("  %-10s %5.1f%% present  %d bins  thr %5.1f  |  24h line %.3e (%.4f%%)  |  above thr: %d"
          %(NICE[c],100*ok.mean(),N,thr,d24,100*d24,len(idx)),flush=True)

print("\nGATE 3 -- do any candidates appear in MORE THAN ONE LINE?")
allc={}
for c,(f,R,thr,C,w,idx) in res.items():
    for i in idx[np.argsort(-R[idx])][:40]:
        allc.setdefault(round(f[i],11),[]).append((NICE[c],R[i]))
multi={k:v for k,v in allc.items() if len(v)>1}
print("   distinct candidate frequencies: %d ; appearing in >1 line: %d"%(len(allc),len(multi)))
for fr,hits in sorted(multi.items(),key=lambda kv:-max(h[1] for h in kv[1]))[:12]:
    print("   %.7e Hz  period %9.4f h  ->  %s"
          %(fr,1.0/fr/3600," ".join("%s:%.0f"%h for h in hits)))
if not multi: print("   none -- every candidate is confined to a single detector channel")

print("\nGATE 2 -- GOES-17")
try:
    G17=load("g17")
    for fr,hits in sorted(multi.items(),key=lambda kv:-max(h[1] for h in kv[1]))[:10]:
        line=[c for c in CH if NICE[c]==hits[0][0]][0]
        r2,ok2=prep(G17[line]); r2,_=daily_profile(r2,ok2); rn2,_=notch(r2,ok2,LINES)
        f2,P2,_=spec(rn2); R2=P2/cont(P2)
        N2=len(R2); thr2=(np.median(R2)/np.log(2.0))*np.log(N2/0.05)
        j=np.argmin(np.abs(f2-fr)); m=R2[max(0,j-2):j+3].max()
        print("   %.7e Hz  %-10s g16 R=%7.1f  g17 R=%7.1f  %s"
              %(fr,hits[0][0],hits[0][1],m,"BOTH" if m>thr2 else "not present"))
    if not multi: print("   nothing reached gate 2")
except FileNotFoundError:
    print("   GOES-17 not available yet")

print("\nTHE 3-5 MINUTE BAND (reported, not gated)")
for c in ("irr_304","irr_1216","MgII_EXIS"):
    if c not in res: continue
    f,R,thr,C,w,idx=res[c]
    for lo_,hi_ in ((1/360.,1/240.),(1/300.,1/180.)):
        m=(f>lo_)&(f<hi_)
        print("   %-10s %4.1f-%4.1f min: max R = %6.1f  (threshold %.1f)"
              %(NICE[c],1/hi_/60,1/lo_/60,R[m].max(),thr))

print("\nSENSITIVITY, Lyman-alpha")
if "irr_1216" in res:
    f,R,thr,C,w,idx=res["irr_1216"]
    for pmin in (2,5,15,60,360,1440):
        i=np.argmin(np.abs(f-1.0/(pmin*60)))
        print("   period %6d min   eps = %.3e  (%.5f%%)"%(pmin,2*np.sqrt(thr*C[i])/np.sum(w),
                                                          100*2*np.sqrt(thr*C[i])/np.sum(w)))
