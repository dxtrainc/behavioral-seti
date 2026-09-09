"""Fetch Wind/MFI h0 over the ACE span and reduce to daily medians of |B|.

Files are streamed, reduced and DELETED -- the analysis needs one number per
day, not 7.6 GB of 3-second vectors. Daily median, >100 samples required, is
exactly the reduction ace_key.py applies to ACE, so the two series are
comparable by construction.
"""
import os, io, sys, json, datetime as dt
import numpy as np, requests, cdflib
from concurrent.futures import ThreadPoolExecutor

BASE="https://cdaweb.gsfc.nasa.gov/pub/data/wind/mfi/mfi_h0"
TMP="/home/dxtra/windtmp"; os.makedirs(TMP, exist_ok=True)
S=requests.Session()

def listing(year):
    r=S.get("%s/%d/"%(BASE,year), timeout=60); r.raise_for_status()
    import re
    return sorted(set(re.findall(r"wi_h0_mfi_(\d{8})_v(\d+)\.cdf", r.text)))

def one(args):
    ymd, ver = args
    url="%s/%s/wi_h0_mfi_%s_v%s.cdf"%(BASE, ymd[:4], ymd, ver)
    fp=os.path.join(TMP, "w_%s.cdf"%ymd)
    try:
        r=S.get(url, timeout=180)
        if r.status_code!=200: return ymd, None, 0
        open(fp,"wb").write(r.content)
        c=cdflib.CDF(fp); i=c.cdf_info()
        zv=list(getattr(i,"zVariables",[]) or [])+list(getattr(i,"rVariables",[]) or [])
        if "BF1" not in zv or "Epoch" not in zv: return ymd, None, 0
        b=np.asarray(c.varget("BF1"), float).ravel()
        m=np.isfinite(b)&(b>0)&(b<1e3)          # Wind flags fill as -1e31
        if m.sum()<100: return ymd, None, int(m.sum())
        return ymd, float(np.median(b[m])), int(m.sum())
    except Exception as e:
        return ymd, None, -1
    finally:
        try: os.remove(fp)
        except Exception: pass

def main():
    todo=[]
    for y in range(2012, 2020):
        try:
            L=listing(y); todo += [(d,v) for d,v in L]
            print("%d: %d files"%(y,len(L)), flush=True)
        except Exception as e:
            print("%d: listing failed %s"%(y,e), flush=True)
    # keep the highest version per day
    best={}
    for d,v in todo: 
        if d not in best or int(v)>int(best[d]): best[d]=v
    todo=sorted(best.items())
    print("total days: %d"%len(todo), flush=True)
    out={}; nb=0
    with ThreadPoolExecutor(max_workers=10) as ex:
        for k,(ymd,med,ns) in enumerate(ex.map(one, todo)):
            if med is not None: out[ymd]=med
            else: nb+=1
            if k%400==0 and k: print("  %d/%d  good=%d"%(k,len(todo),len(out)), flush=True)
    print("daily medians: %d   unusable days: %d"%(len(out), nb), flush=True)
    json.dump(out, open("/home/dxtra/wind_daily.json","w"))
    print("saved ~/wind_daily.json", flush=True)

main()
