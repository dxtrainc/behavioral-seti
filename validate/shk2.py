"""Shklovskii + Galactic-acceleration correction. psrcat stores no GL/GB, so
galactic coordinates are computed from RAJ/DECJ (J2000 -> galactic rotation)."""
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

import tarfile, os, json, math, numpy as np
R0=8.178; TH0=236.0e3; c=2.99792458e8
kpc=3.0856775814913673e19; yr=3.155693e7
mas_yr=(1e-3/3600.0)*(math.pi/180.0)/yr
AG=math.radians(192.85948); DG=math.radians(27.12825); LNCP=math.radians(122.93192)

def hms(s, is_ra):
    p=s.split(":")
    try: v=[float(x) for x in p]
    except ValueError: return None
    while len(v)<3: v.append(0.0)
    sign=-1.0 if s.strip().startswith("-") else 1.0
    a=abs(v[0])+v[1]/60.0+v[2]/3600.0
    a*=sign
    return math.radians(a*15.0) if is_ra else math.radians(a)

def galactic(ra,dec):
    sb=math.sin(dec)*math.sin(DG)+math.cos(dec)*math.cos(DG)*math.cos(ra-AG)
    b=math.asin(max(-1,min(1,sb)))
    y=math.cos(dec)*math.sin(ra-AG)
    x=math.sin(dec)*math.cos(DG)-math.cos(dec)*math.sin(DG)*math.cos(ra-AG)
    l=LNCP-math.atan2(y,x)
    return (math.degrees(l)%360.0), math.degrees(b)

with tarfile.open(_p("psrcat_pkg.tar.gz")) as t:
    for m in t.getmembers():
        if m.name.endswith("psrcat.db"):
            db=t.extractfile(m).read().decode("utf8","replace"); break
recs=[]; cur={}
for line in db.split("\n"):
    if line.startswith("@"):
        if cur: recs.append(cur); 
        cur={}; continue
    p=line.split()
    if len(p)>=2: cur.setdefault(p[0],p[1])
if cur: recs.append(cur)
def f(r,k):
    try: return float(r[k])
    except Exception: return None

def a_z(zkpc):
    z=abs(zkpc)
    return (2.27*z + 3.68*(1.0-math.exp(-4.31*z)))*1e-11

rows=[]; have_pm=0
for r in recs:
    P=f(r,"P0"); Pd=f(r,"P1"); F0=f(r,"F0"); F1=f(r,"F1")
    if P is None and F0: P=1.0/F0
    if Pd is None and F0 and F1 is not None: Pd=-F1/(F0*F0)
    if P is None or Pd is None: continue
    pmra=f(r,"PMRA"); pmdec=f(r,"PMDEC")
    if pmra is None or pmdec is None: continue
    have_pm+=1
    px=f(r,"PX"); d=None; dsrc="dm"
    if px and px>0: d=1.0/px; dsrc="px"
    else:
        for k in ("DIST_A","DIST_DM","DIST_DM1"):
            if f(r,k): d=f(r,k); break
    if not d or d<=0: continue
    if "RAJ" not in r or "DECJ" not in r: continue
    ra=hms(r["RAJ"],True); dec=hms(r["DECJ"],False)
    if ra is None or dec is None: continue
    gl,gb=galactic(ra,dec)
    mu=math.hypot(pmra,pmdec)*mas_yr
    shk=P*mu*mu*(d*kpc)/c
    l=math.radians(gl); b=math.radians(gb)
    beta=(d*math.cos(b))/R0-math.cos(l)
    dr=-(TH0**2/(c*R0*kpc))*(math.cos(l)+beta/(math.sin(l)**2+beta**2))*P
    z=d*math.sin(b)
    vz=-a_z(z)*abs(math.sin(b))/c*P
    intr=Pd-shk-dr-vz
    rows.append(dict(name=r.get("PSRJ","?"),P=P,Pd=Pd,shk=shk,gal=dr+vz,intr=intr,
                     d=d,dsrc=dsrc,gl=gl,gb=gb,
                     frac=((shk+dr+vz)/Pd if Pd>0 else float("nan"))))

msp=[r for r in rows if r["P"]<0.03 and r["Pd"]>0]
print("with proper motion:            %d"%have_pm)
print("fully corrigible (PM+dist+pos): %d   of which millisecond: %d"%(len(rows),len(msp)))
print("  distance from parallax %d, from DM %d"%(sum(1 for r in msp if r["dsrc"]=="px"),
                                                 sum(1 for r in msp if r["dsrc"]=="dm")))
fr=np.array([r["frac"] for r in msp if np.isfinite(r["frac"])])
print("\nkinematic terms as a fraction of OBSERVED Pdot (millisecond pulsars):")
for q in (10,25,50,75,90): print("   %2dth pct %+.3f"%(q,np.percentile(fr,q)))
print("   >100%% (intrinsic would be negative): %d of %d"%(int((fr>1).sum()),len(fr)))
lo_obs=[r for r in msp if r["Pd"]<1e-20]
lo_int=[r for r in msp if 0<r["intr"]<1e-20]
print("\nthe population the hypothesis needs:")
print("   Pdot_obs  < 1e-20 : %d"%len(lo_obs))
print("   Pdot_intr < 1e-20 : %d"%len(lo_int))
if lo_obs:
    print("   median kinematic fraction among Pdot_obs<1e-20: %+.2f"
          %np.median([r["frac"] for r in lo_obs if np.isfinite(r["frac"])]))
tau=[r["P"]/(2*r["intr"])/(1e9*yr) for r in msp if r["intr"]>0]
print("\n   corrected characteristic ages: %d with tau>6 Gyr (was 155 uncorrected)"
      %sum(1 for t_ in tau if t_>6))
json.dump(rows,open(_p("psrcat_shk.json"),"w"))
print("\nsaved ~/psrcat_shk.json (%d rows)"%len(rows))
