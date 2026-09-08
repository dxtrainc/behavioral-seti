"""
THE CONTROL THAT ACTUALLY DISCRIMINATES -- because the last one could not.

The filter-width test passed cleanly: the 135 uHz peak sits at 135.1 uHz on GOES-16
across filter widths from 30 to 714 uHz, and at 135.0 on GOES-17 across 52 to 620.
It does not track the filter.

The phase-randomised surrogate ALSO showed the peak, at 135.3 uHz and the same
significance -- and that is a broken control, not a failed detection. Phase
randomisation preserves |FFT|^2 of the time series EXACTLY. The comb is a feature of
the power spectrum. The statistic is computed from the power spectrum. So the
surrogate preserves the very thing it was meant to destroy, and could never have
discriminated anything here. It was the right null for the multi-scale search, where
the statistic depended on fold phase; it is a null-op for this one. Sixth control in
this programme that had to be thrown away.

THE PHYSICS SUPPLIES A BETTER ONE. Solar p-modes are acoustic waves TRAPPED inside
the Sun, and trapping fails above the acoustic cutoff at ~5.3 mHz. Above the cutoff
waves escape, there are no standing modes, and there is no comb. So:

    2.5 - 4.0 mHz   below cutoff, mode band      -> comb MUST be present
    6.0 - 8.0 mHz   above cutoff, no modes       -> comb MUST be absent

Identical pipeline, identical filter, same record, same instrument -- the only thing
that changes is whether the Sun can trap a wave there. If 135 uHz appears in both
bands the method is manufacturing it. If it appears only below the cutoff, it is the
star, and the control is one no instrument artifact can imitate.
"""
import os
import numpy as np
from scipy.ndimage import median_filter

D=os.path.expanduser("~")
def load(tok,c="MgII_EXIS"):
    z=np.load(D+"/euvs1m_%s.npz"%tok); return z[c]

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

def band_acf(r, lo, hi, width=20001):
    n=len(r)//2*2; w=np.blackman(n)
    P=np.abs(np.fft.rfft(r[:n]*w))**2
    f=np.fft.rfftfreq(n,d=60.0)[1:]; P=P[1:]
    m=(f>=lo-1e-3)&(f<=hi+1e-3)
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
    # value of the ACF exactly at 135 uHz, whatever the peak is
    j=np.argmin(np.abs(lags-135e-6))
    return lags[sel][k], A[sel][k], (A[sel][k]/base if base>0 else np.nan), A[j]

BANDS=[("2.5-4.0 mHz  BELOW cutoff (modes)",2.5e-3,4.0e-3),
       ("6.0-8.0 mHz  ABOVE cutoff (none) ",6.0e-3,8.0e-3),
       ("1.0-2.0 mHz  below mode band     ",1.0e-3,2.0e-3)]
print("ACOUSTIC-CUTOFF CONTROL -- MgII, identical pipeline in each band\n")
print("  %-6s %-36s %10s %8s %8s %10s"%("sat","band","ACF peak","value","sigma","A(135uHz)"))
for tok in ("g16","g17"):
    r,ok=prep(load(tok))
    for name,lo,hi in BANDS:
        lag,val,sig,a135=band_acf(r,lo,hi)
        print("  %-6s %-36s %10.1f %8.3f %8.1f %10.3f"%(tok,name,lag*1e6,val,sig,a135),flush=True)
    print()
print("  p-modes are trapped only below ~5.3 mHz. A comb in the 6-8 mHz band would")
print("  mean the method makes combs; its absence there and presence below is physics.")
