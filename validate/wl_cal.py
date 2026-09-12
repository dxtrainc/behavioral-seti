"""CALIBRATION. A control must be able to fail.

The synthetic null is self-consistent: Y = common_phase(X) has the REAL pairwise
structure of the data (preserved exactly) and, by construction, no four-way phase
alignment. Run the full test on Y using fresh common_phase surrogates. If the null
is matched, p is uniform. If it is not, this is where it shows -- the same place the
Gamma(5) null failed by 3-5x and the circular shift failed by 9x.
"""
import numpy as np, os, sys, time
from scipy.stats import norm, kstest
from multiprocessing import Pool
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "search"))
from wl_null2 import normal_score, common_phase

M=np.load(os.path.expanduser("~/wl_M.npy"))
names=open(os.path.expanduser("~/wl_names.txt")).read().split("\n")
ix={n:i for i,n in enumerate(names)}
Q=["TSI","F10.7","cosmic ray","X-ray bg"]
ov=np.all(np.isfinite(M[[ix[n] for n in Q]]),axis=0)
X=np.vstack([normal_score(M[ix[n]][ov]) for n in Q])
n=X.shape[1]

def stat(v):
    d=np.diff(v); s=d.std()
    if s<=0: return np.nan
    m=np.median(np.abs(d))
    from scipy import stats as st
    return float(st.kurtosis(d,fisher=True)-4.0*m/s)

FORMS={"3rd-diff":  lambda z: z[0]-3*z[1]+3*z[2]-z[3],
       "cross":     lambda z: z[0]-z[1]-z[2]+z[3],
       "quad-prod": lambda z: z[0]*z[1]*z[2]*z[3]}

NSUR=400
def trial(seed):
    rng=np.random.default_rng(seed)
    Y=common_phase(X,rng)                     # data with NO four-way structure
    out={}
    for fn,f in FORMS.items():
        obs=stat(f(Y))
        null=np.array([stat(f(common_phase(Y,rng))) for _ in range(NSUR)])
        null=null[np.isfinite(null)]
        out[fn]=(1.0+int((null>=obs).sum()))/(1.0+len(null))
    return out

if __name__=="__main__":
    NT=int(sys.argv[1]) if len(sys.argv)>1 else 300
    t0=time.time()
    with Pool(80) as p: res=p.map(trial,range(9000,9000+NT))
    print("calibration: %d trials x %d surrogates, %.1f s\n"%(NT,NSUR,time.time()-t0))
    print("  %-11s %8s %8s %8s %8s"%("form","mean p","KS p","frac<.05","expected"))
    for fn in FORMS:
        ps=np.array([r[fn] for r in res]); ps=ps[np.isfinite(ps)]
        ks=kstest(ps,"uniform").pvalue
        flag="" if ks>0.05 else "   <-- NOT UNIFORM"
        print("  %-11s %8.4f %8.2e %8.4f %8.4f%s"%(fn,ps.mean(),ks,(ps<0.05).mean(),0.05,flag))
