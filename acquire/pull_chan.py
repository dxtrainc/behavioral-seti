#!/usr/bin/env python3
"""Acquire additional dimensionless-eligible solar/heliospheric channels.
Every channel is saved as a 2-col npz: jd (float days) + value. Reports coverage."""
import urllib.request, ssl, os, sys, re, io
import numpy as np
ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
D=os.path.expanduser("~/chan"); os.makedirs(D,exist_ok=True)
UA={"User-Agent":"Mozilla/5.0"}
def fetch(u,fn=None,timeout=300):
    fn=fn or os.path.join(D,u.split("/")[-1].split("?")[0])
    if os.path.exists(fn) and os.path.getsize(fn)>1000: return fn
    r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=timeout,context=ctx)
    open(fn,"wb").write(r.read()); return fn
def jd(y,m,d):
    y=np.asarray(y,float); m=np.asarray(m,float); d=np.asarray(d,float)
    a=np.floor((14-m)/12); yy=y+4800-a; mm=m+12*a-3
    return d+np.floor((153*mm+2)/5)+365*yy+np.floor(yy/4)-np.floor(yy/100)+np.floor(yy/400)-32045
OUT={}
def save(name,j,v,note=""):
    j=np.asarray(j,float); v=np.asarray(v,float)
    ok=np.isfinite(j)&np.isfinite(v); j,v=j[ok],v[ok]
    if len(j)==0: print("  %-16s EMPTY"%name); return
    o=np.argsort(j); j,v=j[o],v[o]
    # collapse duplicate days by mean
    uj,inv=np.unique(j,return_inverse=True)
    uv=np.bincount(inv,weights=v)/np.bincount(inv)
    np.savez_compressed(os.path.join(D,name+".npz"),jd=uj,v=uv)
    yr=(uj[-1]-uj[0])/365.25
    OUT[name]=(len(uj),yr)
    print("  %-16s %6d days  %5.1f yr  %s"%(name,len(uj),yr,note))

# ---------------- 1. ACE SWICS charge states & abundance ratios -------------
print("[1] ACE SWICS (ssv4 1998-2011 daily)")
try:
    from pyhdf.HDF import HDF, HC
    from pyhdf.VS import VS
    fn=fetch("https://izw1.caltech.edu/ACE/ASC/DATA/level2/ssv4/ssv4_data_1day.hdf",
             os.path.join(D,"ssv4_1day.hdf"))
    h=HDF(fn,HC.READ); vs=h.vstart(); ref=vs.next(-1); vd=vs.attach(ref)
    nrec,_,fields,_,_=vd.inquire()
    arr=np.array(vd.read(nrec),dtype=object); vd.detach()
    col={f:i for i,f in enumerate(fields)}
    def C(f): return np.array([float(r[col[f]]) for r in arr])
    J=jd(C("year"),1,1)+C("day")-1
    for f,q in [("O7to6","qf_O7to6"),("C6to5","qf_C6to5"),("avqFe","qf_avqFe"),
                ("FetoO","qf_FetoO"),("HetoO","qf_HetoO"),("avqO","qf_avqO"),
                ("CtoO","qf_CtoO"),("avqSi","qf_avqSi")]:
        v=C(f); qq=C(q)
        v[(qq!=0)|(v<=-1e30)|(v<=0)]=np.nan
        save("swics_"+f.lower(),J,v,"qf==0 only")
    vs.end(); h.close()
except Exception as e: print("  FAILED:",repr(e)[:180])

# ---------------- 2. sunspot area, total and hemispheric asymmetry ----------
print("[2] sunspot area (MSFC/RGO+USAF)")
try:
    fn=fetch("https://solarscience.msfc.nasa.gov/greenwch/daily_area.txt")
    d=np.genfromtxt(fn,skip_header=1)
    J=jd(d[:,0],d[:,1],d[:,2]); tot,n,s=d[:,3],d[:,4],d[:,5]
    tot=np.where(tot<0,np.nan,tot); n=np.where(n<0,np.nan,n); s=np.where(s<0,np.nan,s)
    save("spot_area",J,tot,"total corrected area")
    asym=(n-s)/np.where((n+s)>0,(n+s),np.nan)
    save("spot_asym",J,asym,"(N-S)/(N+S) -- dimensionless by construction")
except Exception as e: print("  FAILED:",repr(e)[:180])

# ---------------- 3. WSO mean solar magnetic field --------------------------
print("[3] WSO mean solar magnetic field")
try:
    fn=fetch("http://wso.stanford.edu/meanfld/MF_timeseries.txt")
    J=[];V=[]
    for line in open(fn,errors="replace"):
        m=re.match(r"\s*(\d{4}):(\d{2}):(\d{2})_\d+h\s+([-\d.]+)",line)
        if m:
            J.append(jd(int(m.group(1)),int(m.group(2)),int(m.group(3)))); V.append(float(m.group(4)))
    save("wso_meanfield",J,V,"daily mean line-of-sight field")
except Exception as e: print("  FAILED:",repr(e)[:180])

# ---------------- 4. LASCO CME rate and mean speed --------------------------
print("[4] SOHO/LASCO CME catalogue")
try:
    fn=fetch("https://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/text_ver/univ_all.txt")
    days={}; spd={}
    for line in open(fn,errors="replace"):
        m=re.match(r"\s*(\d{4})/(\d{2})/(\d{2})\s+(\d{2}):(\d{2}):(\d{2})\s+\S+\s+(\S+)\s+(\S+)",line)
        if not m: continue
        k=float(jd(int(m.group(1)),int(m.group(2)),int(m.group(3))))
        days[k]=days.get(k,0)+1
        try:
            v=float(m.group(8))
            if v>0: spd.setdefault(k,[]).append(v)
        except ValueError: pass
    ks=sorted(days)
    save("cme_rate",ks,[days[k] for k in ks],"CMEs per day")
    ks2=sorted(spd)
    save("cme_speed",ks2,[np.mean(spd[k]) for k in ks2],"daily mean linear speed")
except Exception as e: print("  FAILED:",repr(e)[:180])

# ---------------- 5. Ca II K (Sac Peak, via LISIRD) -------------------------
print("[5] Ca II K indices (LISIRD)")
try:
    fn=fetch("https://lasp.colorado.edu/lisird/latis/dap/cak.csv",os.path.join(D,"cak.csv"))
    txt=open(fn,errors="replace").read().strip().split("\n")
    hdr=[h.strip() for h in txt[0].split(",")]
    rows=[r.split(",") for r in txt[1:] if r.count(",")>=len(hdr)-1]
    def cc(name):
        i=hdr.index(name); return np.array([float(r[i]) for r in rows])
    mm=np.array([int(r[0].split()[0]) for r in rows]); dd=np.array([int(r[0].split()[1]) for r in rows])
    yy=np.array([int(r[0].split()[2]) for r in rows])
    J=jd(yy,mm,dd)
    for f in ("emdx","k2vk3","k3","delk1"):
        if f in hdr: save("cak_"+f,J,cc(f),"Sac Peak Ca II K")
except Exception as e: print("  FAILED:",repr(e)[:180])

print("\n=== SUMMARY: %d channels written to ~/chan ==="%len(OUT))
for k,(n,y) in sorted(OUT.items(),key=lambda x:-x[1][1]):
    print("  %-18s %6d days  %5.1f yr"%(k,n,y))
