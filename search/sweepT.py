"""
THE COMBINATION SWEEP.

WHY THIS EXISTS. The twenty-four searches in this programme each picked a channel
opportunistically -- whichever archive was open when the idea arrived -- and asked
whether something was in it. That samples a space without covering it. This sweep
poses the question the other way round: enumerate the parameters, enumerate the
combinations, and test them all against one trials budget.

THE REGULARISER, WHICH IS WHAT MAKES IT TRACTABLE. A sender choosing something to
modulate cannot know our units, our calibration, our zero points, our cadence, or
which proxy we happen to use for a quantity. Whatever they embed must therefore
survive all of that -- which admits only DIMENSIONLESS constructions: ratios, log
ratios, normalised correlations, phase differences. Anything carrying units is
unusable to them. That takes 2^22 subsets down to a few hundred tests. It is not a
convenience; it is what the sender is forced into.

THE NULL. For a combination of two channels the wrong null is independence -- both
are solar-activity driven and genuinely related. The right one is a CIRCULAR SHIFT
of one channel against the other: it preserves each channel's own distribution,
spectrum, red noise and non-Gaussianity exactly, and destroys only the specific
alignment being tested. Same device as the sky scramble and the shifted key.

THE DETECTOR is the quantisation statistic validated earlier in this programme --
MAD-based clipping (std-based fails because the std is set by the outliers being
clipped), then excess kurtosis of increments minus a scaled median absolute
increment. It is the one instrument here with a measured ROC.
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

import os, json, itertools
from multiprocessing import Pool
import numpy as np
from scipy import stats

rng = np.random.default_rng(4242)
D = os.path.expanduser("~")

from scipy.ndimage import median_filter
def runmed(y, w):
    # O(n log w) instead of O(n*w). The pure-python version made the null loop
    # 260 billion operations -- days of runtime for a search meant to take minutes.
    return median_filter(y, size=w, mode="nearest")

def clean(x):
    d=np.diff(x)
    s=1.4826*np.median(np.abs(d-np.median(d)))
    if s<=0: return None
    d=np.clip(d,-3.0*s,3.0*s)
    y=np.concatenate([[x[0]],x[0]+np.cumsum(d)])
    sd=y.std()
    return (y-y.mean())/sd if sd>0 else None

def stat(x):
    z=clean(x)
    if z is None: return np.nan
    d=np.diff(z); s=d.std()
    if s<=0: return np.nan
    return float(stats.kurtosis(d,fisher=True)-4.0*np.median(np.abs(d))/s)

# ---------- load every channel onto one daily JD grid ----------
def daily(jd, val, grid, minpts=1):
    idx=np.digitize(jd,grid)-1
    out=np.full(len(grid)-1,np.nan)
    for k in range(len(grid)-1):
        m=idx==k
        if m.sum()>=minpts: out[k]=np.median(val[m])
    return out

LO,HI = 2444240.5, 2461041.5     # 1980-01-01 .. 2025-12-31
grid = np.arange(np.floor(LO), np.ceil(HI)+1, 1.0)
CH = {}

d=np.loadtxt(_p("tsi.csv"),delimiter=",",skiprows=1)
ok=np.isfinite(d[:,1])&(d[:,1]>1000); CH["TSI"]=daily(d[ok,0]+2309100.5,d[ok,1],grid)

d=np.loadtxt(_p("f107.csv"),delimiter=",",skiprows=1)
ok=np.isfinite(d[:,1])&(d[:,1]>0); CH["F10.7"]=daily(d[ok,0],d[ok,1],grid)

sn=np.genfromtxt(_p("sunspot.csv"),delimiter=";")
m=sn[:,5]>=0; CH["sunspot"]=daily((sn[m,3]-2000.0)*365.25+2451545.0,sn[m,5],grid)

mg=np.loadtxt(_p("mgii.csv"),delimiter=",",skiprows=1)
ok=np.isfinite(mg[:,1])&(mg[:,1]>0); CH["MgII"]=daily(mg[ok,0],mg[ok,1],grid)

# neutron monitor: "YYYY-MM-DD HH:MM:SS; value"
import datetime as dt
T,V=[],[]
for line in open(_p("nm_oulu.txt"),errors="ignore"):
    if not line[:4].isdigit(): continue
    try:
        a,b=line.split(";")
        y,mo,dd=int(a[:4]),int(a[5:7]),int(a[8:10])
        v=float(b)
        T.append((dt.date(y,mo,dd)-dt.date(2000,1,1)).days+2451544.5); V.append(v)
    except Exception: pass
T=np.array(T); V=np.array(V); ok=np.isfinite(V)&(V>0)
CH["cosmic ray"]=daily(T[ok],V[ok],grid)

# OMNI: |B| col 8, density 23, speed 24
import glob
# OMNI2 carries 55 columns. Three of them are DIMENSIONLESS BY CONSTRUCTION --
# alpha/proton ratio, plasma beta, Alfven Mach number -- which is exactly the class
# a sender is restricted to, since they are invariant to our units and calibration.
# Those are natural candidates rather than ratios we impose.
# ALGEBRAIC INDEPENDENCE. A first version included plasma beta (col 36) and the
# Alfven Mach number (col 37) because they are dimensionless -- and they produced
# 16 Bonferroni "survivors", every one of them a pair involving those two. They are
# not measurements: beta = nkT/(B^2/2mu0) and M_A = v/(B/sqrt(mu0 rho)) are
# FUNCTIONS of density, temperature, speed and |B|, all already in the set. Testing
# |B| against beta tests |B| against a formula containing |B|. Dropped, along with
# flow pressure and the vxB electric field for the same reason, and the geomagnetic
# indices, which are a downstream RESPONSE to the solar wind rather than an
# independent observable. Dimensionless is necessary but not sufficient: the
# quantity must also be independently measured.
OM = {"IMF |B|":(8,0,100), "proton temp":(22,1e3,1e7), "wind density":(23,0,200),
      "wind speed":(24,100,3000), "alpha/proton":(27,1e-3,0.3)}
acc={k:[] for k in OM}; JD=[]
for f in sorted(glob.glob(D+"/omni/omni2_*.dat")):
    a=np.genfromtxt(f)
    if a.ndim!=2 or a.shape[1]<45: continue
    yr,doy,hr=a[:,0],a[:,1],a[:,2]
    jd=367*yr-np.floor(7*(yr+np.floor(11/12))/4)+np.floor(275/9)+1721013.5+doy-1+hr/24.0
    JD.append(jd)
    for k,(c,lo_,hi_) in OM.items():
        v=a[:,c].copy()
        v[(v<lo_)|(v>hi_)]=np.nan          # OMNI fill values fall outside physical range
        acc[k].append(v)
jd=np.concatenate(JD)
for k in OM:
    v=np.concatenate(acc[k]); m=np.isfinite(v)
    CH[k]=daily(jd[m],v[m],grid,6)

# ---------------------------------------------------------------- NOAA / GOES
# Three channels the earlier sweep did not have, all from the open NCEI archive.
# Every one is measured at Earth orbit, so each carries a 3.4% peak-to-peak
# annual modulation from the eccentricity of our own orbit -- nothing to do with
# the Sun. Left in, that ONE-YEAR periodicity is common to all three and would
# couple them to each other and to any other uncorrected channel. Correct it.
def au2(jdv):
    """(r/1AU)^2 -- multiply an at-Earth flux by this to put it at 1 AU"""
    g=np.radians(357.529+0.98560028*(jdv-2451545.0))
    r=1.00014-0.01671*np.cos(g)-0.00014*np.cos(2*g)
    return r*r

def csvcol(path, col, lo_, hi_):
    """date,... -> (jd, value); dates are ISO so parse them directly"""
    J,V=[],[]
    for line in open(path):
        f=line.strip().split(",")
        if len(f)<=col or not f[0][:4].isdigit(): continue
        try:
            y,mo,dd=int(f[0][:4]),int(f[0][5:7]),int(f[0][8:10]); v=float(f[col])
        except ValueError: continue
        if not (lo_<v<hi_): continue
        J.append((dt.date(y,mo,dd)-dt.date(2000,1,1)).days+2451544.5); V.append(v)
    return np.array(J), np.array(V)

# X-ray background 1-8 A. 1983-2019 reduced from the 1-minute archive as the daily
# 10th percentile (the standard definition -- immune to flares); 2017-2025 from the
# GOES-16 bkd1d product, which is NOAA's own version of the same quantity.
jx1,vx1 = csvcol(_p("goes_xray_bg.csv"),1,1e-10,1e-4)
jx2,vx2 = csvcol(D+"/goesr/g16_xray_bg.csv",1,1e-10,1e-4)
SPLIT   = 2457791.5                                    # 2017-02-07, GOES-16 start
keep1   = jx1 < SPLIT
jx=np.concatenate([jx1[keep1],jx2]); vx=np.concatenate([vx1[keep1],vx2])
CH["X-ray bg"]=daily(jx,vx*au2(jx),grid)

# Integral proton flux >10 MeV. EPS p3_flux_ic to 2009, EPEAD ZPGT10 E/W after --
# NOAA designed the channels to correspond, so the join is by construction.
jp,vp = csvcol(_p("goes_proton10.csv"),1,1e-4,1e5)
CH["proton >10MeV"]=daily(jp,vp*au2(jp),grid)

# Lyman-alpha 121.6 nm. GOES 13/14/15 EUVE daily 2006-2016, GOES-16 EUVS 2019-2025.
# Sparse, with a three-year hole at the instrument change -- which the pair overlap
# handles and the handover check quantifies.
JL,VL=[],[]
for src in sorted(glob.glob(D+"/goesr/G1?_EUVE_daily_????_2016_v4.txt")):
    c=6                                    # Irrad_ly, the 1 nm Lyman-alpha band
    for line in open(src,errors="ignore"):
        if line.startswith(";") or not line[:4].isdigit(): continue
        f=line.split()
        if len(f)<=7: continue
        try: v=float(f[c]); au=float(f[7])
        except ValueError: continue
        if not (1e-4<v<1e-1): continue
        y,mo,dd=int(line[:4]),int(line[5:7]),int(line[8:10])
        JL.append((dt.date(y,mo,dd)-dt.date(2000,1,1)).days+2451544.5); VL.append(v/au)
jl2,vl2 = csvcol(D+"/goesr/g16_euvs.csv",4,1e-4,1e-1)   # irr_1216
JL=np.concatenate([np.array(JL),jl2]); VL=np.concatenate([np.array(VL),vl2*au2(jl2)])
o=np.argsort(JL)
CH["Lyman-alpha"]=daily(JL[o],VL[o],grid)



# ---------- ADDITIONAL CHANNELS (acquired by ~/pull_chan.py) ----------
# Each is already one value per day, sorted and unique, so the O(n) map below
# replaces daily() -- which loops 16,802 bins per channel and would dominate
# the runtime for seventeen more channels.
def daily_pre(jdv, val, grid):
    out=np.full(len(grid)-1,np.nan)
    idx=(np.round(jdv)-grid[0]).astype(int)
    m=(idx>=0)&(idx<len(out))
    out[idx[m]]=val[m]
    return out

NEW=["spot_area","spot_asym","wso_meanfield",
     "cak_emdx","cak_k2vk3","cak_k3","cak_delk1",
     "cme_rate","cme_speed",
     "swics_o7to6","swics_c6to5","swics_avqfe","swics_avqsi",
     "swics_fetoo","swics_hetoo","swics_avqo","swics_ctoo"]
for _nm in NEW:
    _f=os.path.expanduser("~/chan/%s.npz"%_nm)
    if not os.path.exists(_f):
        print("  MISSING %s"%_nm,flush=True); continue
    _z=np.load(_f)
    CH[_nm]=daily_pre(_z["jd"],_z["v"],grid)

names=list(CH.keys())
M=np.vstack([CH[n] for n in names])
# Requiring every channel simultaneously collapses to the intersection -- 463 days
# out of 16,071, because one short record truncates all eight. Each PAIR is tested
# on its own overlap instead, which is what the question actually needs.
good=np.all(np.isfinite(M),axis=0)
print("channels: %s"%", ".join(names),flush=True)
for n in names:
    print("  %-14s finite %6d / %d"%(n,np.isfinite(CH[n]).sum(),len(grid)-1),flush=True)
print("days with ALL channels: %d  (using PER-PAIR overlap instead)"%good.sum(),flush=True)
N=len(names)

# ---------- dimensionless combinations ----------
W=181
def prep(x):
    """normalise and detrend ONCE. A circular shift of a series has a circularly
    shifted running median, so the detrend never needs recomputing inside the
    null loop -- shift the prepared arrays instead."""
    r=x/np.median(x)
    d=r-runmed(r,W)
    sd=d.std()
    return r, (d/sd if sd>0 else d)

def forms_pre(ra,da,rb,db):
    out={}
    out["ratio"]=ra/rb
    out["log-ratio"]=np.log(np.maximum(ra,1e-12))-np.log(np.maximum(rb,1e-12))
    out["resid-product"]=da*db
    out["resid-ratio"]=da/np.where(np.abs(db)<0.1,np.nan,db)
    return {k:v for k,v in out.items() if np.all(np.isfinite(v))}

import time

# ================== TAIL-FIT ENGINE =================================
# Counting exceedances floors p at 1/(M+1), so M must scale with the test
# count to stay under Bonferroni -- cost N^2. Fitting a generalised Pareto
# to the upper tail of the shift null gives p below that floor from a FIXED
# surrogate budget, restoring linear cost. Knijnenburg et al. 2009,
# Bioinformatics 25, i161.
import argparse
from scipy.stats import genpareto, kstest

# ---- DETERMINISTIC SEEDS. Python hash() is salted per process, so the seed below was
# different in every run and no p-value produced through it could be regenerated:
# three consecutive processes gave 1689220225, 56146563 and 1311111642 for the same
# key. Not a correctness fault -- an arbitrary seed is still a valid seed -- but a
# reproducibility one, and for a paper whose claim is that a null is worth the
# fraction of a space it excludes, an unreproducible null is worth less than it
# looks. beacon.seeds.stable_seed is BLAKE2b over a canonical repr: fixed across
# processes, interpreter versions and platforms.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from beacon.seeds import stable_seed as _stable_seed



def _outpath(name):
    """Bare filename -> written to ~; anything with a directory component or an
    absolute path -> used as given. The old form joined "~/" unconditionally."""
    name = os.path.expanduser(name)
    if os.path.isabs(name) or os.path.dirname(name):
        return name
    return os.path.join(os.path.expanduser("~"), name)


ap=argparse.ArgumentParser()
ap.add_argument("--mode",default="pairs",choices=["pairs","triples"])
ap.add_argument("--shifts",type=int,default=10000)
ap.add_argument("--pool",type=int,default=80)
ap.add_argument("--shard",default="0/1")      # i/n -- take every nth job
ap.add_argument("--minov",type=int,default=2000)
ap.add_argument("--out",default="sweepT.json")
A=ap.parse_args()
NSH=A.shifts
SI,SN=(int(x) for x in A.shard.split("/"))
MINOV=A.minov

def pval(obs,null):
    """Empirical p where the exceedance count resolves it; GPD tail
    extrapolation where it does not. Returns (p, method, gof)."""
    null=np.asarray(null,float); Mn=len(null)
    ne=int((null>=obs).sum())
    emp=(1.0+ne)/(1.0+Mn)
    if ne>=10: return emp,"emp",float("nan")
    for frac in (0.10,0.15,0.20,0.30):
        Nt=int(Mn*frac)
        if Nt<100: continue
        u=float(np.partition(null,Mn-Nt)[Mn-Nt])
        exc=null[null>u]-u
        if len(exc)<100 or obs<=u: continue
        try: c,_,sc=genpareto.fit(exc,floc=0)
        except Exception: continue
        if not np.isfinite(sc) or sc<=0: continue
        try: gof=float(kstest(exc,"genpareto",args=(c,0,sc)).pvalue)
        except Exception: gof=0.0
        if gof<0.05: continue          # bad fit -- widen the tail and retry
        p=(len(exc)/Mn)*float(genpareto.sf(obs-u,c,loc=0,scale=sc))
        if np.isfinite(p) and p>0: return p,"gpd",gof
    return emp,"emp-floor",float("nan")

def rollp(p,s): return (np.roll(p[0],s), np.roll(p[1],s))

def forms_tri(pa,pb,pc):
    """Dimensionless three-body forms. resid-triple is the third-order
    analogue of the pair product and does NOT factor into pair statistics;
    log-2nd-diff is a curvature in log space, with pb as the middle term."""
    ra,da=pa; rb,db=pb; rc,dc=pc
    out={}
    out["resid-triple"]=da*db*dc
    la=np.log(np.maximum(ra,1e-12)); lb=np.log(np.maximum(rb,1e-12)); lc=np.log(np.maximum(rc,1e-12))
    out["log-2nd-diff"]=la-2.0*lb+lc
    return {k:v for k,v in out.items() if np.all(np.isfinite(v))}

def one_pair(arg):
    i,j,fname=arg
    ov=np.isfinite(M[i])&np.isfinite(M[j])
    if ov.sum()<MINOV: return None
    pa=prep(M[i][ov]); pb=prep(M[j][ov]); n=int(ov.sum())
    f=forms_pre(pa[0],pa[1],pb[0],pb[1])
    if fname not in f: return None
    obs=stat(f[fname])
    if not np.isfinite(obs): return None
    rng=np.random.default_rng(_stable_seed(i,j,fname))
    null=[]
    for s in rng.integers(1,n,size=NSH):
        pb2=rollp(pb,int(s))
        fb=forms_pre(pa[0],pa[1],pb2[0],pb2[1])
        if fname not in fb: continue
        v=stat(fb[fname])
        if np.isfinite(v): null.append(v)
    if len(null)<NSH//4: return None
    p,meth,gof=pval(obs,null)
    return dict(combo="%s / %s"%(names[i],names[j]),form=fname,obs=float(obs),
                p=p,method=meth,gof=gof,nnull=len(null),ndays=n)

def one_tri(arg):
    i,j,k,fname=arg
    ov=np.isfinite(M[i])&np.isfinite(M[j])&np.isfinite(M[k])
    if ov.sum()<MINOV: return None
    pa=prep(M[i][ov]); pb=prep(M[j][ov]); pc=prep(M[k][ov]); n=int(ov.sum())
    f=forms_tri(pa,pb,pc)
    if fname not in f: return None
    obs=stat(f[fname])
    if not np.isfinite(obs): return None
    rng=np.random.default_rng(_stable_seed(i,j,k,fname))
    null=[]
    for s1,s2 in zip(rng.integers(1,n,size=NSH),rng.integers(1,n,size=NSH)):
        fb=forms_tri(pa,rollp(pb,int(s1)),rollp(pc,int(s2)))
        if fname not in fb: continue
        v=stat(fb[fname])
        if np.isfinite(v): null.append(v)
    if len(null)<NSH//4: return None
    p,meth,gof=pval(obs,null)
    return dict(combo="%s / %s / %s"%(names[i],names[j],names[k]),form=fname,
                obs=float(obs),p=p,method=meth,gof=gof,nnull=len(null),ndays=n)

# ---- build the job list -------------------------------------------------
jobs=[]
if A.mode=="pairs":
    for i,j in itertools.combinations(range(N),2):
        ov=np.isfinite(M[i])&np.isfinite(M[j])
        if ov.sum()<MINOV: continue
        pa=prep(M[i][ov]); pb=prep(M[j][ov])
        for fn in forms_pre(pa[0],pa[1],pb[0],pb[1]): jobs.append((i,j,fn))
    worker=one_pair
else:
    drop=0
    for i,j,k in itertools.combinations(range(N),3):
        ov=np.isfinite(M[i])&np.isfinite(M[j])&np.isfinite(M[k])
        if ov.sum()<MINOV: drop+=1; continue
        jobs.append((i,j,k,"resid-triple"))
        # log-2nd-diff is not symmetric: each channel gets a turn as the middle
        jobs.append((i,j,k,"log-2nd-diff"))
        jobs.append((j,i,k,"log-2nd-diff"))
        jobs.append((i,k,j,"log-2nd-diff"))
    print("  triples dropped for overlap < %d d: %d of %d"%(MINOV,drop,N*(N-1)*(N-2)//6),flush=True)
    worker=one_tri

jobs=[x for n_,x in enumerate(jobs) if n_%SN==SI]
print("\nmode=%s shard=%d/%d  %d jobs x %d shifts = %.1fM evaluations, pool=%d"
      %(A.mode,SI,SN,len(jobs),NSH,len(jobs)*NSH/1e6,A.pool),flush=True)

t0=time.time()
with Pool(A.pool) as pool:
    R=[r for r in pool.map(worker,jobs,chunksize=1) if r]
el=time.time()-t0
print("completed %d tests in %.1f min  (%.0f evaluations/s)"
      %(len(R),el/60,len(R)*NSH/max(el,1)),flush=True)

json.dump(dict(mode=A.mode,shard=A.shard,shifts=NSH,ntests=len(R),
               elapsed_s=el,results=R),open(_outpath(A.out),"w"))
print("saved ~/%s"%A.out,flush=True)

R.sort(key=lambda r:r["p"])
meth={}
for r in R: meth[r["method"]]=meth.get(r["method"],0)+1
print("  p by method: %s"%meth,flush=True)
print("\n  %-46s %-14s %9s %10s %6s"%("combination","form","stat","p","meth"))
for r in R[:15]:
    print("  %-46s %-14s %9.4f %10.3e %6s"%(r["combo"][:46],r["form"],r["obs"],r["p"],r["method"]),flush=True)
