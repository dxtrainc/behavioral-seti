"""p-mode comb curves. Pipeline identical to pmode_cutoff.py except that the
20001-wide running median is evaluated on a stride and interpolated: the filter
output varies only on the scale of its own window, so sampling it every 400 bins
oversamples its smoothness by 50x while turning an O(n*w)=1e10 operation into
seconds. Validated by whether it reproduces the published 135 uHz peak."""
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

import os, json, time, numpy as np
D=os.path.expanduser("~")

def logmedian_fast(y, size=20001, stride=400):
    n=len(y); h=size//2
    idx=np.arange(0,n,stride)
    if idx[-1]!=n-1: idx=np.append(idx,n-1)
    vals=np.empty(len(idx))
    for j,i in enumerate(idx):
        a=max(0,i-h); b=min(n,i+h+1)
        vals[j]=np.median(y[a:b])
    return np.interp(np.arange(n), idx, vals)

def load(tok,c="MgII_EXIS"):
    z=np.load(_p("euvs1m_%s.npz")%tok); return z[c]

def prep(x, wmin=60):                      # verbatim from pmode_cutoff.py
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

def band_acf_curve(r, lo, hi):
    n=len(r)//2*2; w=np.blackman(n)
    P=np.abs(np.fft.rfft(r[:n]*w))**2
    f=np.fft.rfftfreq(n,d=60.0)[1:]; P=P[1:]
    m=(f>=lo-1e-3)&(f<=hi+1e-3)
    fb=f[m]; Pb=P[m]
    C=np.exp(logmedian_fast(np.log(np.maximum(Pb,1e-300))))
    R=Pb/C
    mm=(fb>=lo)&(fb<=hi)
    x=R[mm]-R[mm].mean(); df=fb[1]-fb[0]
    A=np.fft.irfft(np.abs(np.fft.rfft(x))**2, n=len(x)); A=A[:len(A)//2]/A[0]
    lags=np.arange(len(A))*df
    sel=(lags>60e-6)&(lags<260e-6)
    k=int(np.argmax(A[sel]))
    base=float(np.std(A[(lags>300e-6)&(lags<600e-6)]))
    return (lags[sel]*1e6, A[sel], float(lags[sel][k]*1e6), float(A[sel][k]),
            float(A[sel][k]/base if base>0 else float("nan")))

OUT={}; t0=time.time()
for tok in ("g16","g17"):
    r,_=prep(load(tok))
    for tag,lo,hi in (("mode",2.5e-3,4.0e-3),("cutoff",6.0e-3,8.0e-3)):
        L,A,pk,val,sig=band_acf_curve(r,lo,hi)
        s=max(1,len(L)//440)
        OUT["%s_%s"%(tok,tag)]=dict(lag=[round(float(v),2) for v in L[::s]],
                                    acf=[round(float(v),4) for v in A[::s]],
                                    peak=round(pk,1), val=round(val,4), sigma=round(sig,1))
        print("%-4s %-7s peak %7.1f uHz  A=%.3f  %.1f sigma   [%.0fs]"
              %(tok,tag,pk,val,sig,time.time()-t0),flush=True)
json.dump(OUT,open(_p("figdata_comb.json"),"w"))
print("saved ~/figdata_comb.json")
