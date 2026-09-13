"""
MAVEN/EUVM Level-2 band irradiances -- from PDS PPI at UCLA rather than LASP.

WHY THE SOURCE CHANGED. The LASP SDC throttles hard: the steady-paced fetcher
reached 826 of 1,945 days and then pinned at its 90 s ceiling with 85 refusals,
which is about 28 hours for the remainder. The identical files
(mvn_euv_l2_bands_YYYYMMDD_vNN_rNN.cdf) are mirrored in the PDS Planetary Plasma
Interactions node at UCLA, a different host, which served six consecutive
requests in 0.4 s each with no refusal. Coverage over the GOES span is essentially
complete: 2,452 files across 2019-2025, the only gap being March and most of
April 2022, which is the MAVEN safe-mode outage and is missing from the
instrument, not from the archive.

URLS ARE LISTED, NOT CONSTRUCTED. Version and revision suffixes vary by date and
some months are empty, so guessing "_v18_r01" returns 404s that look like missing
data. Each month directory is listed once and whatever it holds is taken.

Reduction is unchanged: 1-minute means per band, files discarded after use, so
the Earth and Mars series are directly comparable.
"""
import os, re, sys, time, json, datetime as dt
import numpy as np, requests, cdflib
from concurrent.futures import ThreadPoolExecutor

DATA = os.environ.get("BEACON_DATA") or os.path.expanduser("~")   # data root; see README

PPI = "https://pds-ppi.igpp.ucla.edu/data/maven-euv-calibrated/data/bands"
UA = {"User-Agent": "beacon-search/1.0 (rtg@dxtra.com; academic solar irradiance study)"}
OUT = os.path.join(DATA, "maven_l2_1m.json")
TMP = os.path.join(DATA, "mavtmp")
START, NDAYS = dt.date(2019, 12, 10), 1945
os.makedirs(TMP, exist_ok=True)
S = requests.Session(); S.headers.update(UA)

def month_index(y, m):
    try:
        r = S.get("%s/%04d/%02d/" % (PPI, y, m), timeout=90)
        if r.status_code != 200: return {}
    except Exception:
        return {}
    out = {}
    for fn in re.findall(r"mvn_euv_l2_bands_(\d{8})_v(\d+)_r(\d+)\.cdf", r.text):
        d, v, rr = fn
        key = (int(v), int(rr))
        if d not in out or key > out[d][0]:
            out[d] = (key, "mvn_euv_l2_bands_%s_v%s_r%s.cdf" % (d, v, rr))
    return {d: "%s/%04d/%02d/%s" % (PPI, y, m, f) for d, (_, f) in out.items()}

def reduce_day(path):
    c = cdflib.CDF(path); info = c.cdf_info()
    zv = list(getattr(info, "zVariables", []) or []) + \
         list(getattr(info, "rVariables", []) or [])
    tv = next((v for v in ("time_unix", "Epoch", "epoch") if v in zv), None)
    dv = next((v for v in ("data", "DATA", "irradiance") if v in zv), None)
    if tv is None or dv is None: return None
    t = np.asarray(c.varget(tv), float).ravel()
    d = np.asarray(c.varget(dv), float)
    if d.ndim == 1: d = d[:, None]
    if d.shape[0] != len(t): d = d.T
    if tv.lower() == "epoch": t = (t/1000.0) - 62167219200.0
    ok = np.isfinite(t) & (t > 1e9)
    t, d = t[ok], d[ok]
    if len(t) < 10: return None
    m = (t//60).astype(np.int64)
    um, inv = np.unique(m, return_inverse=True)
    out = np.full((len(um), d.shape[1]), np.nan)
    for b in range(d.shape[1]):
        v = d[:, b]
        g = np.isfinite(v) & (v > 0) & (v < 1e30)
        if g.sum() == 0: continue
        sm = np.bincount(inv[g], weights=v[g], minlength=len(um))
        n = np.bincount(inv[g], minlength=len(um))
        out[:, b] = np.where(n > 0, sm/np.maximum(n, 1), np.nan)
    return um, out

def one(args):
    key, url = args
    fp = os.path.join(TMP, "p_%s.cdf" % key)
    try:
        r = S.get(url, timeout=180)
        if r.status_code != 200 or len(r.content) < 5000: return key, None
        open(fp, "wb").write(r.content)
        red = reduce_day(fp)
        if red is None: return key, None
        um, out = red
        return key, {"min": um.tolist(),
                     "v": [[None if not np.isfinite(x) else round(float(x), 12)
                            for x in row] for row in out]}
    except Exception:
        return key, None
    finally:
        try: os.remove(fp)
        except Exception: pass

def main():
    store = json.load(open(OUT)) if os.path.exists(OUT) else {}
    have = {k for k, v in store.items() if v}
    print("resuming: %d days already reduced" % len(have), flush=True)
    want = [(START+dt.timedelta(days=i)).strftime("%Y%m%d") for i in range(NDAYS)]
    months = sorted({(int(k[:4]), int(k[4:6])) for k in want})
    index = {}
    for y, m in months:
        idx = month_index(y, m)
        index.update(idx)
        print("  %04d-%02d: %d files" % (y, m, len(idx)), flush=True)
        time.sleep(0.2)
    todo = [(k, index[k]) for k in want if k in index and k not in have]
    print("to fetch: %d (of %d wanted; %d not present in the archive)"
          % (len(todo), len(want), len(want)-len([k for k in want if k in index])), flush=True)
    t0 = time.time(); n = 0
    with ThreadPoolExecutor(max_workers=4) as ex:
        for key, rec in ex.map(one, todo):
            n += 1
            if rec: store[key] = rec
            if n % 50 == 0:
                json.dump(store, open(OUT, "w"))
                el = time.time()-t0
                print("  %d/%d  have=%d  %.0f d/h" % (n, len(todo),
                      sum(1 for v in store.values() if v), 3600.0*n/max(el, 1)), flush=True)
    json.dump(store, open(OUT, "w"))
    print("done: %d days with data" % sum(1 for v in store.values() if v), flush=True)

main()
