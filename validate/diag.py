"""Why does the reproduction not return the tabulated analytic amplitude?

eps = 2*sqrt(mu*ln(N/alpha)*C)/sum(w). Four inputs can differ between the
original and the reproduction: mu, N, alpha, and -- the suspect -- C, which the
original evaluates AT THE TARGET BIN while the reproduction used the median over
the search band. Solar power is red, so those are not the same number.
"""
import numpy as np, os
from scipy.ndimage import median_filter
def _p(n):
    for c in (n,os.path.join(os.path.expanduser("~"),os.path.basename(n))):
        if os.path.exists(c): return c
    return n
def prep(x,wmin=60):
    ok=np.isfinite(x)&(x>0); y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    m0=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-m0)); y=np.clip(y,m0-5*s,m0+5*s)
    a=np.nan_to_num(y); mk=ok.astype(float)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(mk)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    return np.where(ok,y-np.where(cm[hi]-cm[lo]>5,(ca[hi]-ca[lo])/np.maximum(cm[hi]-cm[lo],1),m0),0.0)

z=np.load(_p("euvs1m_g16.npz")); r=prep(np.asarray(z["irr_1216"],float))
n=len(r)//2*2; w=np.hanning(n); ws=w.sum()
X=np.fft.rfft(r[:n]*w); f=np.fft.rfftfreq(n,d=60.0)[1:]; P=(np.abs(X)**2)[1:]
C_all=np.exp(median_filter(np.log(np.maximum(P,1e-300)),size=20001,mode="nearest"))
R=P/C_all
mu_full=float(np.median(R)/np.log(2.0)); N_full=len(R)
band=(f>=1/180.)&(f<=1/90.)
mu_band=float(np.median(R[band])/np.log(2.0)); N_band=int(band.sum())
i2=int(np.argmin(np.abs(f-1/120.)))                 # the 2-minute bin exactly
def eps(mu,N,alpha,C): return 2.0*np.sqrt(mu*np.log(N/alpha)*C)/ws
print("tabulated analytic value: 4.0e-07\n")
print("%-52s %10s"%("variant","eps"))
for lbl,mu,N,al,C in (
  ("full spectrum, mu empirical, a=0.05, C at 2-min bin",  mu_full,N_full,0.05,C_all[i2]),
  ("  same but C = median over the 90-180 s band",         mu_full,N_full,0.05,float(np.median(C_all[band]))),
  ("  same but N = band bins only",                        mu_full,N_band,0.05,C_all[i2]),
  ("  same but alpha = 0.01",                              mu_full,N_full,0.01,C_all[i2]),
  ("  same but mu = 1 (Exp(1) assumed)",                   1.0,    N_full,0.05,C_all[i2]),
  ("band-only mu and N, C at bin  [the reproduction]",     mu_band,N_band,0.05,float(np.median(C_all[band]))),
):
    print("%-52s %10.3e"%(lbl,eps(mu,N,al,C)))
print("\nC at the 2-minute bin      %.4e"%C_all[i2])
print("C median over 90-180 s     %.4e   ratio %.2f (amplitude %.2f)"
      %(np.median(C_all[band]),np.median(C_all[band])/C_all[i2],np.sqrt(np.median(C_all[band])/C_all[i2])))
print("mu full %.3f   mu band %.3f   N full %d   N band %d"%(mu_full,mu_band,N_full,N_band))
