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

import tarfile, os, json, math
tgz=_p("psrcat_pkg.tar.gz")
with tarfile.open(tgz) as t:
    for m in t.getmembers():
        if m.name.endswith("psrcat.db"):
            db=t.extractfile(m).read().decode("utf8","replace"); break
recs=[]; cur={}
for line in db.split("\n"):
    if line.startswith("@"):
        if cur: recs.append(cur)
        cur={}; continue
    p=line.split()
    if len(p)>=2: cur.setdefault(p[0],p[1])
if cur: recs.append(cur)
def f(r,k):
    try: return float(r[k])
    except Exception: return None
out=[]
for r in recs:
    P=f(r,"P0"); Pd=f(r,"P1")
    F0=f(r,"F0"); F1=f(r,"F1")
    if P is None and F0: P=1.0/F0                      # period from frequency
    if Pd is None and F0 and F1 is not None: Pd=-F1/(F0*F0)   # Pdot = -F1/F0^2
    if P is None or Pd is None: continue
    dist=None
    for k in ("DIST","DIST_DM","DIST_A","DIST_AMN","DIST_DM1"):
        if f(r,k): dist=f(r,k); break
    out.append(dict(name=r.get("PSRJ") or r.get("PSRB") or "?", P=P, Pd=Pd,
                    pmra=f(r,"PMRA"), pmdec=f(r,"PMDEC"), px=f(r,"PX"), dist=dist,
                    assoc=r.get("ASSOC",""), binary=r.get("BINARY","")))
json.dump(out,open(_p("psrcat.json"),"w"))
msp=[x for x in out if x["P"]<0.03 and x["Pd"]>0]
field=[x for x in msp if "GC" not in (x["assoc"] or "")]
lo=[x for x in field if x["Pd"]<1e-20]
lo2=[x for x in field if x["Pd"]<3e-21]
print("records                     %5d"%len(recs))
print("with usable P and Pdot      %5d"%len(out))
print("MSPs  P<30 ms, Pdot>0       %5d"%len(msp))
print("  field (non-globular)      %5d"%len(field))
print("  field, Pdot < 1e-20       %5d   <-- retains most of P0 over 10 Gyr"%len(lo))
print("  field, Pdot < 3e-21       %5d"%len(lo2))
print("  with proper motion        %5d"%sum(1 for x in field if x["pmra"] is not None))
print("  with parallax or distance %5d"%sum(1 for x in field if x["px"] or x["dist"]))
print("saved ~/psrcat.json")
