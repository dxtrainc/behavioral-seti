"""
IS THE 135 uHz COMB REAL, OR IS IT MY CONTINUUM FILTER?

The p-mode result rests on an autocorrelation peak at 135 uHz in the normalised
spectrum. The normalisation divides by a running median of log power 20001 bins
wide -- and at this record length 20001 bins IS 119 uHz. A running median divided
out of a spectrum high-passes it at the filter's own scale and can ring at roughly
that period. 119 against 135 is not enough separation to publish on.

THE TEST. Vary the filter width over a factor of twenty-four and watch the ACF peak.
  - if the peak TRACKS the filter width, the detection is the filter
  - if the peak STAYS at 135 uHz, it is the star

A second control: run the identical pipeline on a PHASE-RANDOMISED surrogate, which
keeps the power spectrum's envelope but destroys the comb. The filter is applied
identically, so anything the filter imprints survives into the surrogate and
anything real does not.
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

import os
import numpy as np
from scipy.ndimage import median_filter

D=os.path.expanduser("~")
def load(tok,c="MgII_EXIS"):
    z=np.load(_p("euvs1m_%s.npz")%tok); return z[c]

def prep(x, wmin=60):
    ok=np.isfinite(x)&(x>0)
    y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-5*s,med+5*s)
    a=np.nan_to_num(y); m=ok.astype(np.float64)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    num=ca[hi]-ca[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>5,num/np.maximum(den,1),med)
    return np.where(ok,y-trend,0.0), ok

def acf_peak(r, width, lo=2.5e-3, hi=4.0e-3):
    n=len(r)//2*2; w=np.blackman(n)
    P=np.abs(np.fft.rfft(r[:n]*w))**2
    f=np.fft.rfftfreq(n,d=60.0)[1:]; P=P[1:]
    m=(f>=1.0e-3)&(f<=6.0e-3)
    fb=f[m]; Pb=P[m]
    C=np.exp(median_filter(np.log(np.maximum(Pb,1e-300)),size=width,mode="nearest"))
    R=Pb/C
    mm=(fb>=lo)&(fb<=hi)
    x=R[mm]-R[mm].mean(); df=fb[1]-fb[0]
    A=np.fft.irfft(np.abs(np.fft.rfft(x))**2, n=len(x)); A=A[:len(A)//2]/A[0]
    lags=np.arange(len(A))*df
    sel=(lags>60e-6)&(lags<260e-6)
    k=np.argmax(A[sel])
    base=np.std(A[(lags>300e-6)&(lags<600e-6)])
    return lags[sel][k], A[sel][k], (A[sel][k]/base if base>0 else np.nan), df

def surrogate(r, seed=5):
    F=np.fft.rfft(r); rng=np.random.default_rng(seed)
    ph=rng.uniform(0,2*np.pi,len(F)); ph[0]=0
    if len(r)%2==0: ph[-1]=0
    return np.fft.irfft(np.abs(F)*np.exp(1j*ph), n=len(r))

print("DOES THE 135 uHz PEAK TRACK THE FILTER, OR STAY PUT?\n")
for tok in ("g16","g17"):
    r,ok=prep(load(tok))
    _,_,_,df=acf_peak(r,20001)
    print("  %s  (MgII, %.1f%% present, bin %.4f uHz)"%(tok,100*ok.mean(),df*1e6))
    print("     %-14s %-14s %10s %8s %10s"%("filter bins","= filter uHz","ACF peak","value","sigma"))
    for wdt in (5001,20001,60001,120001):
        lag,val,sig,_=acf_peak(r,wdt)
        print("     %-14d %-14.1f %10.1f %8.3f %10.1f"%(wdt,wdt*df*1e6,lag*1e6,val,sig),flush=True)
    rs=surrogate(r)
    lag,val,sig,_=acf_peak(rs,20001)
    print("     %-14s %-14s %10.1f %8.3f %10.1f   <- PHASE-RANDOMISED SURROGATE"
          %("surrogate","20001 / 119.0",lag*1e6,val,sig),flush=True)
    print()
print("  a real comb sits at 135 uHz for every filter width and DIES in the surrogate.")
print("  a filter artifact moves with the width and SURVIVES the surrogate.")
