"""
SELF-KEYED SEARCH: the detection procedure of section 3.5 applied to the
alignment-specific excess reported by ace_key.py (|rho| = 2.95, p = 0.005).

The excess was previously set aside because the correlation is NEGATIVE and a
beacon keyed to solar activity would not anticorrelate with it. That is using
the framework to interpret the result -- the error section 5.3 warns against.
This script instead runs the paper's own procedure: matched null, replication on
an independent instrument, then physics.

FOUR TESTS.

 1 ROTATION-MATCHED NULL. F10.7 and the interplanetary field both carry solar
   rotation: F10.7 because active regions rotate into and out of view, the field
   because corotating stream structure sweeps past the spacecraft. The published
   null shifts the key by arbitrary lags, which destroys the beacon alignment
   AND the rotational phase alignment together. It therefore cannot separate
   them. A null shifted by INTEGER CARRINGTON ROTATIONS preserves rotational
   phase while destroying arbitrary alignment. If the excess is rotational it
   vanishes against this null; if it survives, it is not rotation.

 2 SPLIT-HALF. Sign and magnitude in 2012-2015 against 2016-2019.

 3 INDEPENDENT MAGNETOMETER. Wind/MFI over the same span, same reduction, same
   key, same statistic. ACE and Wind are different spacecraft carrying different
   instruments at different points in the solar wind.

 4 TRIALS. The look-elsewhere denominator, stated explicitly.
"""
import os, json, glob, datetime as dt
import numpy as np

DATA = os.environ.get("BEACON_DATA") or os.path.expanduser("~")   # data root; see README

CARR = 27.2753          # synodic Carrington rotation, days

def detrend(y, w=61):
    n = len(y); out = np.empty(n); h = w // 2
    for i in range(n):
        lo, hi = max(0, i-h), min(n, i+h+1)
        out[i] = np.median(y[lo:hi])
    return y - out

def chips(flux):
    d = np.diff(flux); c = np.sign(d); c[c == 0] = 1
    return c

def corr(x, c):
    n = min(len(x), len(c)); a, b = x[:n], c[:n]
    a = a - a.mean(); sa = a.std()
    if sa <= 0: return 0.0
    return float(np.dot(a/sa, b)/np.sqrt(n))

def load_f107(path=None):
    if path is None:
        path = os.path.join(DATA, "f107.csv")
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    jd, flux = d[:, 0], d[:, 1]
    ok = np.isfinite(flux) & (flux > 0)
    return jd[ok], flux[ok]

def build(daily_jd, daily_b, jd_f, flux):
    """the ace_key.py pipeline, verbatim in its arithmetic"""
    good = np.isfinite(daily_b)
    fk = np.interp(daily_jd, jd_f, flux)
    filled = daily_b.copy(); filled[~good] = np.median(daily_b[good])
    x = detrend(filled); x[~good] = 0.0
    c = chips(fk)
    n = min(len(x)-1, len(c))
    x, c, gm = x[1:n+1], c[:n], good[1:n+1]
    return x[gm], c[gm]

def nulls(x, c, mode, step=7):
    n = len(c)
    if mode == "arbitrary":
        sh = [s for s in range(30, n-30, step)]
    else:                                   # integer Carrington rotations
        sh, k = [], 1
        while True:
            s = int(round(k*CARR))
            if s > n-30: break
            sh += [s, n-s]
            k += 1
        sh = sorted(set(s for s in sh if 30 <= s <= n-30))
    v = np.array([abs(corr(x, np.roll(c, s))) for s in sh])
    return v[np.isfinite(v)]

def report(tag, x, c, out):
    obs = corr(x, c)
    rec = {"n": int(len(x)), "rho": float(obs), "sign": "neg" if obs < 0 else "pos"}
    print("\n--- %s\n    n=%d   rho=%+.5f" % (tag, len(x), obs), flush=True)
    for mode in ("arbitrary", "carrington"):
        v = nulls(x, c, mode)
        if len(v) < 20:
            print("    %-11s too few shifts (%d)" % (mode, len(v)), flush=True); continue
        p = (1 + (v >= abs(obs)).sum())/(1 + len(v))
        print("    %-11s %3d shifts   med %.4f  95%% %.4f  max %.4f   p=%.4f  %s"
              % (mode, len(v), np.median(v), np.percentile(v, 95), v.max(), p,
                 "<-- EXCESS" if p < 0.05 else "no excess"), flush=True)
        rec[mode] = {"p": float(p), "n_shift": int(len(v)),
                     "med": float(np.median(v)), "p95": float(np.percentile(v, 95)),
                     "max": float(v.max())}
    out[tag] = rec
    return rec

def to_series(dd):
    ymd = sorted(dd)
    jd = np.array([dt.date(int(s[:4]), int(s[4:6]), int(s[6:8])).toordinal()+1721424.5
                   for s in ymd])
    b = np.array([dd[s] for s in ymd])
    full = np.arange(jd.min(), jd.max()+1)
    out = np.full(len(full), np.nan)
    out[np.searchsorted(full, jd)] = b
    return full, out

def ace_daily():
    cache = os.path.join(DATA, "ace_daily.json")
    if os.path.exists(cache):
        return json.load(open(cache))
    import cdflib
    d = {}
    files = sorted(glob.glob(os.path.join(DATA, "ace", "*.cdf")))
    print("reducing %d ACE files" % len(files), flush=True)
    for i, f in enumerate(files):
        try:
            c = cdflib.CDF(f); inf = c.cdf_info()
            zv = list(getattr(inf, "zVariables", []) or []) + \
                 list(getattr(inf, "rVariables", []) or [])
            bn = next((v for v in ("Magnitude", "BGSEc", "B_gse", "Bmag") if v in zv), None)
            if bn is None: continue
            b = np.asarray(c.varget(bn), float)
            if b.ndim == 2: b = np.linalg.norm(b, axis=1)
            m = np.isfinite(b) & (b > 0) & (b < 1e3)
            if "Q_FLAG" in zv:
                q = np.asarray(c.varget("Q_FLAG"), float).ravel()
                if len(q) == len(m): m &= (q == 0)
            if m.sum() < 100: continue
            d[os.path.basename(f).split("_")[3]] = float(np.median(b[m]))
        except Exception:
            pass
        if i % 500 == 0 and i: print("  %d" % i, flush=True)
    json.dump(d, open(cache, "w"))
    return d

def main():
    out = {}
    jd_f, flux = load_f107()

    d = ace_daily()
    print("ACE daily medians: %d" % len(d), flush=True)
    jd_a, b_a = to_series(d)
    xa, ca = build(jd_a, b_a, jd_f, flux)
    report("ACE 2012-2019 (published result)", xa, ca, out)

    h = len(xa)//2
    report("ACE first half", xa[:h], ca[:h], out)
    report("ACE second half", xa[h:], ca[h:], out)

    wf = os.path.join(DATA, "wind_daily.json")
    if os.path.exists(wf):
        w = json.load(open(wf))
        print("\nWind daily medians: %d" % len(w), flush=True)
        jd_w, b_w = to_series(w)
        xw, cw = build(jd_w, b_w, jd_f, flux)
        report("WIND/MFI 2012-2019 (independent instrument)", xw, cw, out)
        hw = len(xw)//2
        report("WIND first half", xw[:hw], cw[:hw], out)
        report("WIND second half", xw[hw:], cw[hw:], out)
    else:
        print("\nWind daily medians not present -- skipped", flush=True)

    json.dump(out, open(os.path.join(DATA, "selfkey_replicate.json"), "w"), indent=1)
    print("\nsaved ~/selfkey_replicate.json", flush=True)

if __name__=="__main__": main()
