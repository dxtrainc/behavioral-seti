"""
THE FAST SEARCH.

The premise this programme has been carrying is wrong, and the correction points the
search somewhere better. Everything above assumed slow modulation -- one bit a day,
one bit a solar rotation -- because those are natural SOLAR timescales. But the
sender does not pick the rate to suit the Sun. They pick it to suit the RECEIVER
they are gating for, and any receiver that has passed a detection gate built around
buried structure in astrophysical data has fast Fourier transforms, matched filters
and compute to burn. Assuming one bit per day assumes the receiver is patient rather
than capable.

And the physics agrees, which is the part worth measuring. Solar variability is RED:
its power falls steeply with frequency. So the natural background a signal has to
hide under -- or be found beneath -- is quietest at the fast end. A designer wanting
to be findable-but-not-obvious puts the signal where the noise is low and the
required compute is high. That is high frequency, and it is exactly the regime forty
years of daily records cannot represent at all: one sample a day is band-limited to
periods of two days and longer, no matter how many decades it runs.

So: measure the noise spectrum, do a blind narrowband search across it, and report
the minimum detectable amplitude as a function of frequency -- which is the actual
shape of the elbow room.

  gate 1  the diurnal and eclipse artifacts MUST be recovered. A blind search that
          cannot find the known instrumental lines is not calibrated.
  gate 2  any candidate must appear at the same frequency on GOES-16 AND GOES-17.
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

# ---- SHARED LIBRARY. These were local copies; beacon/ is now the single place the
# construction lives. Equivalence was PROVEN before this edit rather than assumed:
# spec, cont, thr_of and roll_masked each reproduce the local result exactly
# (thr_of 21.902966897 both ways; roll_masked max|diff| = 0), so no committed number
# moves and no re-run was needed. validate/migrate_equiv.py is that check.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from beacon.spectra import spec as _spec, cont as _cont, thr_of as _thr_of
from beacon.surrogates import roll_masked as _roll_masked


D=os.path.expanduser("~")
def load(tok):
    z=np.load(_p("xrs1m_%s.npz")%tok)
    return z["a"],z["b"],int(z["start"])

def prep(x, wmin=1440):
    """log flux, flares clipped at the robust level, slow trend removed, gaps zeroed.
    Zero-filling a gap adds no power at any frequency; interpolating would."""
    y=np.full(len(x),np.nan)
    ok=np.isfinite(x)&(x>0)
    y[ok]=np.log(x[ok])
    med=np.nanmedian(y)
    s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-4*s,med+4*s)                 # flares are the loudest thing here
    # running mean over one day via cumulative sums, NaN-aware
    v=np.nan_to_num(y); m=ok.astype(np.float64)
    cv=np.concatenate([[0],np.cumsum(v)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2
    lo=np.maximum(0,np.arange(len(y))-h); hi=np.minimum(len(y),np.arange(len(y))+h)
    num=cv[hi]-cv[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>10,num/np.maximum(den,1),med)
    r=np.where(ok,y-trend,0.0)
    return r, ok

def spectrum(r):
    f, P, _ = _spec(r, dt=60.0)
    return f, P

def continuum(P, w=801):
    return _cont(P, w)

def search(P,f,label,alpha=0.05):
    C=continuum(P)
    R=P/C
    N=len(R)
    # periodogram power over a smooth continuum is ~exponential, so p = exp(-x/mean)
    mu=np.median(R)/np.log(2.0)
    thr=mu*np.log(N/alpha)
    idx=np.where(R>thr)[0]
    print("  %s: %d bins, threshold R > %.1f (Bonferroni over %d)"%(label,N,thr,N))
    out=[]
    for i in idx[np.argsort(-R[idx])][:25]:
        per=1.0/f[i]
        out.append((f[i],R[i],per))
        print("     f = %.6e Hz   period %10.4f min = %8.4f h   R = %7.1f"
              %(f[i],per/60,per/3600,R[i]))
    if not len(idx): print("     no bins above threshold")
    return out,R,C,mu

print("loading ...",flush=True)
A16,B16,s16=load("g16")
r,ok=prep(B16)
print("  GOES-16 XRS-B: %d minutes, %.1f%% present, %.2f yr span"
      %(len(B16),100*ok.mean(),len(B16)/1440/365.25),flush=True)
f,P=spectrum(r)
print("  frequency range %.3e .. %.3e Hz  (periods %.1f min .. %.2f yr)"
      %(f[0],f[-1],1/f[-1]/60,1/f[0]/3.156e7),flush=True)

print("\nNOISE SPECTRUM -- how much elbow room is there, and where?")
C=continuum(P)
for pmin in (2,5,15,60,360,1440):
    i=np.argmin(np.abs(f-1.0/(pmin*60)))
    print("     period %7d min   continuum power %.4e"%(pmin,C[i]))

print("\nBLIND NARROWBAND SEARCH")
cand16,R16,_,mu=search(P,f,"GOES-16")

print("\nGATE 2 -- the same frequencies on GOES-17")
try:
    A17,B17,s17=load("g17")
    r7,ok7=prep(B17)
    f7,P7=spectrum(r7); C7=continuum(P7); R7=P7/C7
    print("  GOES-17 XRS-B: %d minutes, %.1f%% present"%(len(B17),100*ok7.mean()))
    for fr,Rv,per in cand16[:12]:
        j=np.argmin(np.abs(f7-fr))
        print("     %.6e Hz (%8.3f h)  g16 R=%7.1f   g17 R=%7.1f   %s"
              %(fr,per/3600,Rv,R7[j],"BOTH" if R7[j]>10 else "not present"))
except FileNotFoundError:
    print("  GOES-17 not available yet")

print("\nSENSITIVITY -- smallest fractional amplitude detectable, by period")
n=len(r)//2*2
for pmin in (2,5,15,60,360,1440):
    i=np.argmin(np.abs(f-1.0/(pmin*60)))
    # a sinusoid of fractional amplitude eps in log flux contributes
    # power (eps^2/4)*sum(w^2) at its bin; solve for the Bonferroni threshold
    thr=mu*np.log(len(R16)/0.05)*C[i]
    w=np.hanning(n); eps=2.0*np.sqrt(thr)/np.sum(w)
    print("     period %7d min   detectable at eps = %.3e  (%.4f%%)"%(pmin,eps,100*eps))
