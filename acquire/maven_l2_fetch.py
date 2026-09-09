"""
MAVEN/EUVM Level-2 BAND IRRADIANCES -- the non-Earth-line series for item A.

WHY NOT THE PRODUCTS ALREADY PULLED. L3B is FISM-M *model* output (it carries a
MODEL_UNCERTAINTY array) and the model is partly driven by Earth-based inputs, so
using it for a viewpoint test would contaminate the non-Earth line with the very
line it is being compared against and manufacture a false "coherent at all
viewpoints" result. L2B is genuine measurement but orbit-averaged at 3.66 h,
which cannot reach the two-minute band where section 4.2 sets its deepest limits.
L2 `bands` is the real high-cadence diode product and is what the test needs.

SPAN. Matched to GOES-16 EUVS: 2019-12-10 to 2025-04-06, 1,945 days.

ON BEING A GOOD CITIZEN. LASP rate-limits hard and returns an HTML "Too Many
Requests" page with a 200 status -- so a naive fetcher silently writes 512 bytes
of HTML over every file and the corruption only shows up much later. This fetcher
therefore: goes one file at a time, sleeps between requests, DETECTS the HTML
page by magic number rather than status code, backs off exponentially when it
appears, and records progress after every file so it can resume.
"""
import os, sys, time, json, datetime as dt
import numpy as np, requests, cdflib

BASE = "https://lasp.colorado.edu/maven/sdc/public/data/sci/euv/l2"
UA = {"User-Agent": "beacon-search/1.0 (rtg@dxtra.com; academic solar irradiance study)"}
OUT = "/home/dxtra/maven_l2_1m.json"
TMP = "/home/dxtra/mavtmp"
START = dt.date(2019, 12, 10)
NDAYS = 1945
DELAY = 2.5           # polite baseline between requests
os.makedirs(TMP, exist_ok=True)
S = requests.Session(); S.headers.update(UA)

def looks_html(path):
    try:
        with open(path, "rb") as f: return f.read(16).lstrip()[:5].lower() == b"<html"
    except Exception:
        return True

def fetch(url, path, tries=6):
    """returns True on a real file; backs off on the HTML rate-limit page"""
    wait = 30.0
    for k in range(tries):
        try:
            r = S.get(url, timeout=180)
            if r.status_code == 404: return False
            open(path, "wb").write(r.content)
            if not looks_html(path): return True
            print("    rate-limited; sleeping %.0fs" % wait, flush=True)
            time.sleep(wait); wait = min(wait*2, 900)
        except Exception as e:
            print("    %s; sleeping %.0fs" % (type(e).__name__, wait), flush=True)
            time.sleep(wait); wait = min(wait*2, 900)
    return False

_names = {}
def reduce_day(path):
    """1-minute means of each band, so the series is directly comparable with
    the GOES EUVS 1-minute product the Earth line uses"""
    c = cdflib.CDF(path); info = c.cdf_info()
    zv = list(getattr(info, "zVariables", []) or []) + \
         list(getattr(info, "rVariables", []) or [])
    if not _names:
        _names["all"] = zv
        print("  L2 variables: %s" % zv, flush=True)
    tv = next((v for v in ("time_unix", "Epoch", "epoch", "TIME_UNIX") if v in zv), None)
    dv = next((v for v in ("data", "DATA", "irradiance", "IRRADIANCE") if v in zv), None)
    if tv is None or dv is None: return None
    t = np.asarray(c.varget(tv), float).ravel()
    d = np.asarray(c.varget(dv), float)
    if d.ndim == 1: d = d[:, None]
    if d.shape[0] != len(t): d = d.T
    if tv in ("Epoch", "epoch"):                       # CDF ms -> unix
        t = (t/1000.0) - 62167219200.0
    ok = np.isfinite(t) & (t > 1e9)
    t, d = t[ok], d[ok]
    if len(t) < 10: return None
    m = (t//60).astype(np.int64)
    um, inv = np.unique(m, return_inverse=True)
    nb = d.shape[1]
    out = np.full((len(um), nb), np.nan)
    for b in range(nb):
        v = d[:, b]
        g = np.isfinite(v) & (v > 0) & (v < 1e30)
        if g.sum() == 0: continue
        s = np.bincount(inv[g], weights=v[g], minlength=len(um))
        n = np.bincount(inv[g], minlength=len(um))
        out[:, b] = np.where(n > 0, s/np.maximum(n, 1), np.nan)
    return um, out

def main():
    print('cooling off for 5 min before first request', flush=True)
    time.sleep(300)
    store = json.load(open(OUT)) if os.path.exists(OUT) else {}
    print("resuming with %d days already reduced" % len(store), flush=True)
    nb_seen = None
    for i in range(NDAYS):
        day = START + dt.timedelta(days=i)
        key = day.strftime("%Y%m%d")
        if key in store: continue
        got = False
        for ver in ("v18_r01", "v18_r00", "v17_r01", "v19_r01", "v18_r02"):
            url = "%s/%04d/%02d/mvn_euv_l2_bands_%s_%s.cdf" % (BASE, day.year, day.month, key, ver)
            fp = os.path.join(TMP, "m_%s.cdf" % key)
            if fetch(url, fp):
                try:
                    red = reduce_day(fp)
                    if red is not None:
                        um, out = red
                        store[key] = {"min": um.tolist(),
                                      "v": [[None if not np.isfinite(x) else round(float(x), 12)
                                             for x in row] for row in out]}
                        nb_seen = out.shape[1]; got = True
                except Exception as e:
                    print("  %s reduce failed: %s" % (key, e), flush=True)
                finally:
                    try: os.remove(fp)
                    except Exception: pass
                if got: break
        if not got: store[key] = None
        if i % 25 == 0:
            json.dump(store, open(OUT, "w"))
            done = sum(1 for v in store.values() if v)
            print("  %s  %d/%d  reduced=%d  bands=%s" % (key, i+1, NDAYS, done, nb_seen), flush=True)
        time.sleep(DELAY)
    json.dump(store, open(OUT, "w"))
    print("done: %d days with data of %d" % (sum(1 for v in store.values() if v), NDAYS), flush=True)

main()
