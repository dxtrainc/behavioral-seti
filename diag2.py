"""Why the reproduction misses the tabulated analytic amplitude.

eps = 2*sqrt(mu*ln(N/alpha)*C)/sum(w). The suspect is C: fastsearch.py evaluates
it AT THE TARGET BIN, the reproduction used the median over the search band, and
solar power is red so those differ. Continuum computed with the strided median
that reproduced the mode comb exactly, since the exact filter does not finish.
"""
import numpy as np, os
def _p(n):
    for c in (n,os.path.join(os.path.expanduser("~"),os.path.basename(n))):
        if os.path.exists(c): return c
    return n
def logmed_fast(y,size=20001,stride=400):
    n=len(y); h=size//2; idx=np.arange(0,n,stride)
    if idx[-1]!=n-1: idx=np.append(idx,n-1)
    v=np.array([np.median(y[max(0,i-h):min(n,i+h+1)]) for i in idx])
    return np.interp(np.arange(n),idx,v)
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
C=np.exp(logmed_fast(np.log(np.maximum(P,1e-300))))
R=P/C
mu_f=float(np.median(R)/np.log(2.0)); N_f=len(R)
band=(f>=1/180.)&(f<=1/90.)
mu_b=float(np.median(R[band])/np.log(2.0)); N_b=int(band.sum())
i2=int(np.argmin(np.abs(f-1/120.)))
eps=lambda mu,N,al,c: 2.0*np.sqrt(mu*np.log(N/al)*c)/ws
print("tabulated analytic value for Ly-alpha at 2 min: 4.0e-07\n")
print("%-54s %10s"%("variant","eps"))
for lbl,mu,N,al,c in (
 ("full spectrum, mu=median/ln2, a=0.05, C at the 2-min bin", mu_f,N_f,0.05,C[i2]),
 ("  C = median over the 90-180 s band instead",              mu_f,N_f,0.05,float(np.median(C[band]))),
 ("  N = band bins only",                                     mu_f,N_b,0.05,C[i2]),
 ("  alpha = 0.01",                                           mu_f,N_f,0.01,C[i2]),
 ("  mu = 1",                                                 1.0, N_f,0.05,C[i2]),
 ("the earlier reproduction (band mu, band N, band-median C)",mu_b,N_b,0.05,float(np.median(C[band]))),
): print("%-54s %10.3e"%(lbl,eps(mu,N,al,c)))
print("\nC at 2-min bin %.4e | C median over band %.4e | ratio %.2f -> amplitude x%.2f"
      %(C[i2],np.median(C[band]),np.median(C[band])/C[i2],np.sqrt(np.median(C[band])/C[i2])))
print("mu full %.3f  mu band %.3f  N full %d  N band %d"%(mu_f,mu_b,N_f,N_b))
