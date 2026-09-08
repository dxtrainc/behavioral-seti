"""
NEUTRON MONITOR FAST SEARCH, VERSION 3 -- THE CLEAN ERA.

Version 2 reached 4.2e-5 across sixty-two years and was still not a null: sixty bins
stayed above threshold, every one at an exact harmonic of one hour. That comb is the
archive's own time structure -- data reported on hour boundaries, and eras of
different native resolution spliced onto a single grid. Subtracting each fundamental
did not clear it, because the comb's amplitude changes across the era boundaries and
one sinusoid per line cannot follow that.

So this version drops the spliced eras entirely. Oulu and Kiel from 2000, where both
run at a uniform one-minute cadence throughout -- 26 years, 13.7 million samples,
one sampling grid. Shorter baseline, but the band now reaches a 2-minute period
instead of 10, and there is no splice for a comb to live on.

THE HOUR HARMONICS ARE DELIBERATELY NOT SUBTRACTED HERE. If they were an artifact of
splicing they should be absent from a uniform grid, and removing them by hand would
hide the one diagnostic that settles it. Only the solar family is notched.

A NOTE ON WHAT COHERENT INTEGRATION SELECTS FOR. Version 2 failed to recover the
sidereal anisotropy at 23.9345 h despite nominally reaching its amplitude. The reason
is physical and it matters more than the gate did: the sidereal signal's amplitude and
phase track the 22-year heliospheric polarity cycle, so it is not phase-stable across
decades, and coherent integration destroys exactly that. Deep coherent integration is
therefore a FILTER -- it rejects natural modulations whose phase drifts and keeps
signals whose phase holds. For a search after an engineered beacon that is not a
limitation, it is the point.
"""

import os
import numpy as np
from scipy.ndimage import median_filter

D=os.path.expanduser("~")
def load(st):
    z=np.load(D+"/nm_%s_60s.npz"%st); return z["v"],int(z["y0"]),int(z["step"])

def prep(v,step,wsec=86400):
    ok=np.isfinite(v)&(v>0)
    y=np.full(len(v),np.nan); y[ok]=np.log(v[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-6*s,med+6*s)
    w=int(wsec//step); a=np.nan_to_num(y); m=ok.astype(np.float64)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=w//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    num=ca[hi]-ca[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>10,num/np.maximum(den,1),med)
    return np.where(ok,y-trend,0.0), ok

def notch(r, ok, step, periods):
    """least-squares removal of known deterministic lines, in the time domain.
    Only the covered samples enter the fit, so gaps cannot fabricate a component."""
    t=np.arange(len(r),dtype=np.float64)*step
    cols=[]
    for P in periods:
        cols += [np.cos(2*np.pi*t/P), np.sin(2*np.pi*t/P)]
    A=np.vstack(cols).T
    coef,*_=np.linalg.lstsq(A[ok],r[ok],rcond=None)
    out=r-A@coef
    out[~ok]=0.0
    return out, coef

SOLAR=[86400.0/k for k in (1,2,3,4,5,6)]          # diurnal and harmonics
HOURLY=[3600.0/k for k in (1,2,3,4,5,6)]          # the archive's own hour grid
LINES=SOLAR                                       # hour comb NOT removed: the point is whether it is even there

def spec(r,step,win):
    n=len(r)//2*2
    P=np.abs(np.fft.rfft(r[:n]*win[:n]))**2
    f=np.fft.rfftfreq(n,d=float(step))
    return f[1:],P[1:]
def cont(P,w=801):
    return np.exp(median_filter(np.log(np.maximum(P,1e-300)),size=w,mode="nearest"))

V,Y0,STEP=load("OULU"); r,ok=prep(V,STEP)
n=len(r)//2*2
win=np.blackman(n)                                 # -58 dB sidelobes, vs -31 for Hanning
f0,P0=spec(r,STEP,win); R0=P0/cont(P0)
rn,coef=notch(r,ok,STEP,LINES)
f,P=spec(rn,STEP,win); C=cont(P); R=P/C
N=len(R); mu=np.median(R)/np.log(2.0); thr=mu*np.log(N/0.05)

print("OULU (clean era, 1-min)  %.1f yr, %.1f%% present, %d bins, threshold R > %.1f"%(len(V)*STEP/3.156e7,100*ok.mean(),N,thr))
print("\nLINES REMOVED (amplitude in fractional units)")
for i,P_ in enumerate(LINES):
    amp=np.hypot(coef[2*i],coef[2*i+1])
    print("   %8.4f h   %.3e  (%.4f%%)"%(P_/3600,amp,100*amp))

print("\nGATE 1 -- with the solar line subtracted, is the SIDEREAL anisotropy visible?")
for name,per in (("solar diurnal    24.0000 h",86400.0),("sidereal diurnal 23.9345 h",86164.09)):
    ff=1.0/per; j=np.argmin(np.abs(f-ff))
    print("   %-28s before: R = %8.1f    after notch: R = %8.1f"%(name,R0[j],R[j]))
j=np.argmin(np.abs(f-1.0/86164.09)); k=slice(j-30,j+31)
jm=np.argmax(R[k])+j-30
print("   strongest bin within +/-30 of sidereal: %+d bins, R = %.1f"%(jm-j,R[jm]))

print("\nBLIND SEARCH ON THE RESIDUAL")
idx=np.where(R>thr)[0]
cand=[]
for i in idx[np.argsort(-R[idx])][:20]:
    cand.append((f[i],R[i],1.0/f[i]))
    print("   f = %.7e Hz   period %11.5f h   R = %8.1f"%(f[i],1.0/f[i]/3600,R[i]))
print("   bins above threshold: %d of %d  (expected by chance: 0.05)"%(len(idx),N))
if not len(idx): print("   nothing above threshold")

print("\nGATE 2 -- Kiel, same threshold, +/-2 bins")
V2,_,S2=load("KIEL"); r2,ok2=prep(V2,S2)
n2=len(r2)//2*2; win2=np.blackman(n2)
rn2,_=notch(r2,ok2,S2,LINES)
f2,P2=spec(rn2,S2,win2); C2=cont(P2); R2=P2/C2
N2=len(R2); thr2=(np.median(R2)/np.log(2.0))*np.log(N2/0.05)
print("   KIEL %.1f%% present, threshold R > %.1f"%(100*ok2.mean(),thr2))
for fr,Rv,per in cand[:15]:
    j=np.argmin(np.abs(f2-fr)); m=R2[max(0,j-2):j+3].max()
    print("   %.7e Hz (%10.5f h)  OULU R=%8.1f  KIEL R=%8.1f  %s"
          %(fr,per/3600,Rv,m,"BOTH" if m>thr2 else "not present"))
if not cand: print("   no candidates to test")

print("\nSENSITIVITY (Nyquist edge excluded -- the continuum estimate is unreliable there)")
for pmin in (2,5,15,30,60,180,360,1440):
    i=np.argmin(np.abs(f-1.0/(pmin*60)))
    eps=2.0*np.sqrt(thr*C[i])/np.sum(win)
    print("   period %6d min   eps = %.3e  (%.5f%%)"%(pmin,eps,100*eps))
