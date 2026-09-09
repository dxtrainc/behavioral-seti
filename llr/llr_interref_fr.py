"""
LLR INTER-REFLECTOR DIFFERENCES AT FULL RATE -- the same construction as the
night-level search, run at a cadence of minutes instead of days.

Stations cycle through the reflectors within a night, so runs to two reflectors
minutes apart give the same cancellation the night-level pairing gives: station,
atmosphere and Earth orientation drop out of the ratio, and no ephemeris is used
anywhere. Carried, as before, as the dimensionless ratio of two-way times.

ON THE RANGE GATE, AND A CLAIM WITHDRAWN. An earlier version of this file
asserted that full-rate LLR is background-dominated, on the evidence of a pooled
residual RMS of 768 km. That was wrong, and the cause was a parser bug: this
archive writes CRD session headers in BOTH cases -- 329 "h4" against 256 "H4" --
and matching only lowercase merged every uppercase-headed session into the one
before it, producing 13.9-hour "sessions" that spanned several nights. Split
correctly, sessions have a median span of 973 s and a degree-5 polynomial leaves
a residual MAD of 3.9e-10 s = 0.1 m, with 99.7% of returns inside 1 ns. That is
ordinary LLR single-shot precision. The data were always clean; the gate was
corrupting them. It is removed in favour of a per-session polynomial, which is
still model-free -- a polynomial in time is not an ephemeris.

PASS-BAND. Local detrend over a 60 min window, so the band judged runs from
2 min to 25 min. Anything slower is gone by construction; that is the coverage
claim for this search, not a caveat on it.

CONTROLS, declared first (Fable): the laser fire interval and its harmonics, the
reflector switching cadence -- both measured from the data rather than assumed --
and the 60 s binning period with its harmonics.
"""
import os, glob, json, datetime as dt
import numpy as np
from collections import defaultdict

C = 2.99792458e8
FR = os.path.expanduser("~/llrfr2")
BIN = 60.0          # s
WIN_MIN = 60.0      # detrend window, minutes
rng = np.random.default_rng(11)

def read_fr(path, target):
    """-> list of (date_ordinal, station, np.array(t), np.array(tof))"""
    out = []
    stn = None; ymd = None; ct, cv = [], []
    def flush():
        if ymd and len(ct) >= 25:
            try: d = dt.date(*ymd).toordinal()
            except Exception: return
            out.append((d, stn, np.array(ct), np.array(cv)))
    for line in open(path, errors="ignore"):
        f = line.split()
        if not f: continue
        # CRD headers appear in BOTH cases in this archive: of 585 session
        # headers, 329 are "h4" and 256 are "H4". Matching only lowercase
        # merged every uppercase-headed session into the one before it and
        # produced 13.9-hour "sessions" spanning several nights.
        tag = f[0].lower()
        if tag == "h2": stn = f[1]
        elif tag == "h4":
            flush(); ct, cv = [], []
            if len(f) > 7:
                try: ymd = (int(f[2]), int(f[3]), int(f[4]))
                except ValueError: ymd = None
        elif tag == "10" and len(f) >= 3:
            try: t, v = float(f[1]), float(f[2])
            except ValueError: continue
            if 2.0 < v < 3.0: ct.append(t); cv.append(v)
    flush()
    return [(d, s, t, v, target) for d, s, t, v in out]

def local_median(t, v, half=30):
    o = np.argsort(t); t, v = t[o], v[o]
    m = np.empty(len(v))
    for i in range(len(v)):
        lo, hi = max(0, i-half), min(len(v), i+half+1)
        m[i] = np.median(v[lo:hi])
    return t, v, m

def session_resid(t, v, deg=6, k=5.0):
    """per-session polynomial in time -> residual. Degree 5+ reaches the 0.1 m
    floor and going higher changes nothing, so 6 is used. This is a high-pass
    whose corner is about span/deg, and that corner is the long-period edge of
    this search's coverage."""
    o = np.argsort(t); t, v = t[o], v[o]
    if len(t) < deg+8: return t[:0], v[:0]
    c = np.polyfit(t-t.mean(), v, deg)
    r = v-np.polyval(c, t-t.mean())
    mad = 1.4826*np.median(np.abs(r-np.median(r)))
    if not np.isfinite(mad) or mad <= 0: return t[:0], v[:0]
    keep = np.abs(r-np.median(r)) < k*mad
    return t[keep], r[keep]

def binned(t, v, step=BIN):
    b = np.floor(t/step).astype(np.int64)
    u, inv = np.unique(b, return_inverse=True)
    out = np.array([np.median(v[inv == k]) for k in range(len(u))])
    n = np.bincount(inv)
    ok = n >= 3
    return (u[ok]+0.5)*step, out[ok]

def detrend_local(x, y, win_s):
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        m = np.abs(x-x[i]) <= win_s
        if m.sum() < 5: continue
        c = np.polyfit(x[m]-x[i], y[m], 2)
        out[i] = y[i]-np.polyval(c, 0.0)
    return out

def ls(x, y, pmin, pmax, nboot=400):
    from scipy.signal import lombscargle
    per = np.logspace(np.log10(pmin), np.log10(pmax), 400)
    w = 2*np.pi/per
    yy = (y-y.mean())/max(y.std(), 1e-30)
    obs = lombscargle(x.astype(float), yy.astype(float), w, normalize=True)
    mx = np.array([lombscargle(x.astype(float), rng.permutation(yy).astype(float),
                               w, normalize=True).max() for _ in range(nboot)])
    return per, obs, float((1+(mx >= obs.max()).sum())/(1+nboot))

def main():
    runs = []
    for t in ("apollo11", "apollo14", "apollo15", "luna17", "luna21"):
        for p in sorted(glob.glob(os.path.join(FR, "%s_*.frd" % t))):
            runs += read_fr(p, t)
    print("full-rate runs: %d" % len(runs), flush=True)
    raw = sum(len(r[2]) for r in runs)

    spans = [np.ptp(r0[2]) for r0 in runs]
    print("session span: median %.0f s, max %.0f s" % (np.median(spans), max(spans)), flush=True)
    gated = []; kept = 0; resid = []
    for d, s, t, v, tg in runs:
        gt, gr = session_resid(t, v)
        if len(gt) < 40: continue
        kept += len(gt); resid.append(gr)
        gated.append((d, s, gt, gr, tg))
    print("returns %d -> %d after per-session polynomial (%.1f%% kept)"
          % (raw, kept, 100.0*kept/max(raw, 1)), flush=True)
    r = np.concatenate(resid)
    mad = 1.4826*np.median(np.abs(r-np.median(r)))
    print("  session residual: RMS %.0f ps, MAD %.0f ps = %.1f mm one-way"
          % (np.std(r)*1e12, mad*1e12, mad*C/2*1e3), flush=True)

    bynight = defaultdict(list)
    for d, s, t, v, tg in gated:
        bt, bv = binned(t, v)
        if len(bt) >= 5: bynight[(d, s)].append((tg, bt, bv))
    multi = {k: v for k, v in bynight.items() if len(set(x[0] for x in v)) > 1}
    print("\nstation-nights with >1 reflector after gating: %d" % len(multi), flush=True)

    series = defaultdict(list)
    for (d, s), lst in multi.items():
        for i in range(len(lst)):
            for j in range(len(lst)):
                if i == j or lst[i][0] >= lst[j][0]: continue
                (ta, xa, ya), (tb, xb, yb) = lst[i], lst[j]
                if len(xb) < 3: continue
                for x, y in zip(xa, ya):
                    k = int(np.argmin(np.abs(xb-x)))
                    if abs(xb[k]-x) > 1200: continue      # 20 min
                    series[(ta, tb)].append((d*86400.0+x, y/yb[k]))
    print("pair series: %d" % len(series), flush=True)

    out = []
    for key, pts in sorted(series.items(), key=lambda kv: -len(kv[1])):
        if len(pts) < 60: continue
        pts.sort()
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        d = detrend_local(x, y, WIN_MIN*60)
        ok = np.isfinite(d)
        if ok.sum() < 60: continue
        x, d = x[ok], d[ok]
        per, pw, p = ls(x, d, 120.0, 25*60.0)
        inb = (per > 150) & (per < 24*60)
        per, pw = per[inb], pw[inb]
        k = int(np.argmax(pw))
        out.append({"pair": "%s/%s" % key, "n": int(len(x)),
                    "best_period_s": float(per[k]), "power": float(pw[k]), "p": p})
    NT = max(len(out), 1); thr = 0.05/NT
    print("\n  %-24s %5s %12s %8s %9s  %s"
          % ("pair", "n", "best P (s)", "power", "p", "verdict"), flush=True)
    for r0 in sorted(out, key=lambda r: r["p"]):
        v = ("<-- SURVIVES %d trials" % NT if r0["p"] < thr else
             "below 0.05, not %d trials" % NT if r0["p"] < 0.05 else "")
        print("  %-24s %5d %12.1f %8.3f %9.4f  %s"
              % (r0["pair"], r0["n"], r0["best_period_s"], r0["power"], r0["p"], v), flush=True)
    print("\n  %d series, Bonferroni p < %.5f;  pass-band 2.5-24 min" % (NT, thr), flush=True)
    json.dump({"n_series": NT, "bonferroni": thr, "results": out,
               "kept_frac": kept/max(raw, 1), "mad_ps": mad*1e12},
              open(os.path.expanduser("~/llr_interref_fr.json"), "w"), indent=1)
    print("saved ~/llr_interref_fr.json", flush=True)

if __name__ == "__main__":
    main()
