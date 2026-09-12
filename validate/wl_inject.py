"""INJECTION. A limit is not a limit until something injected comes back out.

A STRUCTURAL POINT FOUND WHILE DESIGNING THIS, which belongs in the paper.
Suppose a sender adds c_k*s(t) to channel k. A linear four-body form picks up
(c_a - 3c_b + 3c_c - c_d)s; a linear pair form picks up (c_i - c_j)s. For the signal
to be invisible to EVERY pair, all c_k must be equal -- and then the four-body form
picks up (1-3+3-1)s = 0. So a LINEAR four-body form cannot carry a signal that pairs
cannot also see. Only the multiplicative form can: a four-way phase alignment lives
in the trispectrum and leaves every pairwise cross-correlation untouched. The linear
forms are therefore higher-order differences, not independent carriers, and the
quadruple space is genuinely new only through quad-prod.

So two injections are run:
  A. one-channel quantised modulation -- gives a sensitivity in fractional channel
     units, comparable with the rest of the programme. Visible to pairs too; stated.
  B. four-way-exclusive sign alignment, using the 8 even-parity sign patterns on
     four channels. Every pairwise product averages to zero over the injection
     times while the four-way product is always +1. This is the signal ONLY a
     four-body form can see, and it is what the quadruple space is for.

The injected period is 307 d -- incommensurate with the day, the year and the
27-d solar rotation, after the earlier lesson that a 600 s injection sat on the
144th harmonic of the daily fold and was annihilated by it.
"""
import numpy as np, os, sys, time, json
from multiprocessing import Pool
import sys as _s, os as _o
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "search"))
import wl_sweep as W

NSUR=int(os.environ.get("NSUR","400"))
QUAD=(0,1,2,4)                      # TSI / F10.7 / cosmic ray / X-ray bg
P=307.0; LEV=6                      # quantised sinusoid, 6 levels
t=np.arange(W.N,dtype=float)
SIG=np.round(LEV*np.sin(2*np.pi*t/P))/LEV

EVEN=np.array([s for s in np.ndindex(2,2,2,2)],dtype=float)*2-1
EVEN=np.array([s for s in EVEN if np.prod(s)>0])          # 8 even-parity patterns

def inject(X,mode,eps,rng):
    Y=X.copy()
    if mode=="A":
        Y[QUAD[0]]=Y[QUAD[0]]+eps*SIG
    else:
        pick=EVEN[rng.integers(0,len(EVEN),W.N)]
        for r,k in enumerate(QUAD):
            Y[k]=Y[k]+eps*pick[:,r]
    return Y

def one(arg):
    mode,eps,seed=arg
    rng=np.random.default_rng(seed)
    Y=inject(W.X,mode,eps,rng)
    res={}
    for fn in W.FORMS:
        obs=W.stat(W.FORMS[fn](*[Y[i] for i in QUAD]))
        null=np.array([W.stat(W.FORMS[fn](*[Z[i] for i in QUAD]))
                       for Z in (W.common_phase(Y,rng) for _ in range(NSUR))])
        null=null[np.isfinite(null)]
        res[fn]=(1.0+(null>=obs).sum())/(1.0+len(null))
    return mode,eps,res

if __name__=="__main__":
    NREP=int(sys.argv[1]) if len(sys.argv)>1 else 24
    EPS=[0.0,0.02,0.05,0.10,0.20,0.35,0.50,0.75,1.00]
    jobs=[(m,e,1000+7919*i) for m in ("A","B") for e in EPS for i in range(NREP)]
    print("quadruple: %s"%" / ".join(W.WL[i] for i in QUAD))
    print("%d jobs (%d amplitudes x %d reps x 2 modes), %d surrogates each\n"%(len(jobs),len(EPS),NREP,NSUR))
    t0=time.time()
    with Pool(80) as p: out=p.map(one,jobs,chunksize=1)
    print("%.1f s\n"%(time.time()-t0))
    ALPHA=0.05/len(W.TESTS)
    for mode,lbl in (("A","one-channel quantised modulation"),("B","four-way-exclusive alignment")):
        print("  === mode %s: %s ==="%(mode,lbl))
        print("  %8s %12s %12s %12s"%("eps","3rd-diff","cross","quad-prod"))
        for e in EPS:
            rows=[r for m,ee,r in out if m==mode and ee==e]
            cells=[]
            for fn in W.FORMS:
                ps=np.array([r[fn] for r in rows])
                cells.append("%.2f (%.3f)"%((ps<ALPHA).mean(),np.median(ps)))
            print("  %8.2f %12s %12s %12s"%(e,*cells))
        print()
    json.dump([{"mode":m,"eps":e,"p":r} for m,e,r in out],
              open(os.path.expanduser("~/wl_inject.json"),"w"),indent=1)
    print("recovery fraction (median p) at Bonferroni alpha=%.2e"%ALPHA)
