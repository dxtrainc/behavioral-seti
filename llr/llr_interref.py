"""
LLR INTER-REFLECTOR DIFFERENCES -- the combination-space version, at both scales.

Fable's design, and it is a better search than five separate residual
periodograms. Ranging to ONE reflector gives a number contaminated by station
coordinates, atmosphere, Earth orientation and the lunar ephemeris. Ranging to
TWO from the same station on the same night and taking the difference cancels
all of those to first order. What survives is lunar orientation and anything
acting differently on two points of the Moon. Five reflectors give ten pairs.

DIMENSIONLESS CARRIER, per section 3.1: the pair is carried as the RATIO of the
two two-way times of flight, which is dimensionless by construction and needs no
unit shared with anyone.

TWO SCALES FROM ONE CODE PATH, per Fable's follow-up:
  night   normal points, one pairing per station-night, cadence of days
  session full-rate runs, stations cycle reflectors within a night, so the same
          construction runs at a cadence of minutes

NO DYNAMICAL MODEL IS USED ANYWHERE. Both scales detrend locally. That matters:
published LLR residuals are post-fit against DE or INPOP with tidal
dissipation, libration, station coordinates and Earth orientation all adjusted,
and a signal shaped like any fitted parameter is removed before anyone sees it.
Nothing here is fitted, so nothing here is absorbed.

AND THE PRICE OF THAT: a local detrend is a high-pass filter. Anything longer
than the detrend window is gone by construction, so the window IS the coverage
claim and is stated with every number below.

PHYSICS-SUPPLIED CONTROLS, declared before the search rather than after
(sections 4.4 and 4.11). A survivor at any of these is instrument or celestial
mechanics, not a sender:
    anomalistic month   27.5545 d      draconic month   27.2122 d
    synodic month       29.5306 d      sidereal month   27.3217 d
    206-day term        205.9 d        annual           365.25 d
    nodal cycle         18.613 yr      full-moon return deficit (sampling)
    laser fire interval and harmonics  reflector switching cadence
    libration-driven pulse spread (varies at the libration periods above)
"""
import os, glob, json
import numpy as np

C = 2.99792458e8
REFL = ["apollo11", "apollo14", "apollo15", "luna17", "luna21"]
CONTROLS = {"anomalistic": 27.5545, "draconic": 27.2122, "sidereal": 27.3217,
            "synodic": 29.5306, "206-day": 205.9, "annual": 365.25,
            "half-anomalistic": 27.5545/2, "half-synodic": 29.5306/2,
            "nodal": 18.613*365.25}

def read_npt(path, target):
    """(mjd_like_day, sec_of_day, tof, station) for every normal point"""
    out = []
    stn, y, mo, dy = None, None, None, None
    for line in open(path, errors="ignore"):
        f = line.split()
        if not f: continue
        tag = f[0].lower()          # the archive mixes header case
        if tag == "h2": stn = f[1]
        elif tag == "h4" and len(f) > 7:
            try: y, mo, dy = int(f[2]), int(f[3]), int(f[4])
            except ValueError: pass
        elif tag == "11" and len(f) >= 3 and y:
            try: t, v = float(f[1]), float(f[2])
            except ValueError: continue
            if not (2.0 < v < 3.0): continue
            import datetime as dt
            try: day = dt.date(y, mo, dy).toordinal()
            except Exception: continue
            out.append((day, t, v, stn, target))
    return out

def load_np(root="~/llrnp"):
    root = os.path.expanduser(root)
    rows = []
    for t in REFL:
        for p in sorted(glob.glob(os.path.join(root, "%s_*.npt" % t))):
            rows += read_npt(p, t)
    return rows

def pair_series(rows, max_sep_s=3600.0):
    """for each station-night, pair every two reflectors observed close in time
    and carry the DIMENSIONLESS ratio of their times of flight"""
    from collections import defaultdict
    byday = defaultdict(list)
    for day, t, v, stn, tgt in rows:
        byday[(day, stn)].append((t, v, tgt))
    series = defaultdict(list)
    npairs = 0
    for (day, stn), obs in byday.items():
        obs.sort()
        for i in range(len(obs)):
            for j in range(i+1, len(obs)):
                if obs[j][0]-obs[i][0] > max_sep_s: break
                a, b = obs[i], obs[j]
                if a[2] == b[2]: continue
                key = tuple(sorted((a[2], b[2])))
                r = (a[1]/b[1]) if key[0] == a[2] else (b[1]/a[1])
                series[(key, stn)].append((day + a[0]/86400.0, r))
                npairs += 1
    return series, npairs

def detrend_local(x, y, win_days):
    """local linear over +/- win_days. The window is the high-pass corner and
    therefore the long-period edge of this search's coverage."""
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        m = np.abs(x-x[i]) <= win_days
        if m.sum() < 4: continue
        c = np.polyfit(x[m]-x[i], y[m], 1)
        out[i] = y[i]-np.polyval(c, 0.0)
    return out

def ls_search(x, y, pmin, pmax, nboot=400, seed=1):
    from scipy.signal import lombscargle
    rng = np.random.default_rng(seed)
    per = np.logspace(np.log10(pmin), np.log10(pmax), 600)
    w = 2*np.pi/per
    yy = (y-y.mean())/max(y.std(), 1e-30)
    obs = lombscargle(x.astype(float), yy.astype(float), w, normalize=True)
    # NULL: keep the real epochs, permute the values. Preserves the sampling
    # pattern exactly -- circular shift would assume even sampling, which LLR
    # emphatically is not.
    mx = np.empty(nboot)
    for k in range(nboot):
        mx[k] = lombscargle(x.astype(float), rng.permutation(yy).astype(float),
                            w, normalize=True).max()
    p = (1 + (mx >= obs.max()).sum())/(1+nboot)
    return per, obs, float(p), float(np.percentile(mx, 99))

def near_control(P, tol=0.02, kmax=14):
    """A control is not just its fundamental. The first run flagged a survivor
    at 3.426 d as clean; the sidereal month over 8 is 3.415 d, 0.3% away. Every
    lunar period is checked to its 14th harmonic, because a high-pass filter
    with a 40 d corner suppresses the fundamentals and passes exactly these."""
    best = None
    for name, Pc in CONTROLS.items():
        for k in range(1, kmax+1):
            Pk = Pc/k
            if Pk < 0.5: break
            d = abs(P-Pk)/Pk
            if d < tol and (best is None or d < best[0]):
                best = (d, "%s/%d = %.4f d" % (name, k, Pk) if k > 1 else name)
    return best[1] if best else None

def main():
    rows = load_np()
    print("normal points loaded: %d" % len(rows), flush=True)
    if not rows:
        print("no data yet"); return
    from collections import Counter
    print("  by reflector: %s" % dict(Counter(r[4] for r in rows)), flush=True)
    print("  by station:   %s" % dict(Counter(r[3] for r in rows).most_common(6)), flush=True)

    series, npairs = pair_series(rows)
    print("\ninter-reflector pairings: %d, across %d (pair, station) series"
          % (npairs, len(series)), flush=True)

    WIN = 40.0     # days: the high-pass corner
    print("detrend window %.0f d  ->  pass-band judged 1.5 d to %.0f d\n" % (WIN, WIN/2-1.0),
          flush=True)
    res = []
    for (key, stn), pts in sorted(series.items(), key=lambda kv: -len(kv[1])):
        if len(pts) < 60: continue
        pts.sort()
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        d = detrend_local(x, y, WIN)
        ok = np.isfinite(d)
        if ok.sum() < 60: continue
        x, d = x[ok], d[ok]
        # Search 1.2-20 d but JUDGE only 1.5-19 d. A Lomb-Scargle peak that
        # lands on the first or last frequency of the band is an edge artefact,
        # not a detection -- the first run returned two "survivors" at exactly
        # P = 2.000 d, which was its lower limit.
        per, pw, p, thr99 = ls_search(x, d, 1.2, WIN/2)
        inb = (per > 1.5) & (per < WIN/2-1.0)
        per, pw = per[inb], pw[inb]
        k = int(np.argmax(pw))
        ctl = near_control(per[k])
        res.append({"pair": "%s/%s" % key, "station": stn, "n": int(len(x)),
                    "best_period_d": float(per[k]), "power": float(pw[k]),
                    "p": p, "control": ctl})
    # Trials: one test per (pair, station) series, counted and applied.
    NT = len(res)
    thr = 0.05/max(NT, 1)
    print("  %-22s %-5s %5s %10s %8s %9s  %s"
          % ("pair", "stn", "n", "best P (d)", "power", "p", "verdict"), flush=True)
    for r in sorted(res, key=lambda r: r["p"]):
        if r["control"]:
            v = "CONTROL: %s" % r["control"]
        elif r["p"] < thr:
            v = "<-- SURVIVES %d trials" % NT
        elif r["p"] < 0.05:
            v = "below 0.05 but not %d trials" % NT
        else:
            v = ""
        print("  %-22s %-5s %5d %10.3f %8.3f %9.4f  %s"
              % (r["pair"], r["station"], r["n"], r["best_period_d"],
                 r["power"], r["p"], v), flush=True)
    print("\n  %d series tested, Bonferroni threshold p < %.5f" % (NT, thr), flush=True)
    json.dump({"window_days": WIN, "n_pairings": npairs, "n_trials": NT,
           "bonferroni": thr, "results": res},
              open(os.path.expanduser("~/llr_interref.json"), "w"), indent=1)
    print("\nsaved ~/llr_interref.json", flush=True)

if __name__ == "__main__":
    main()
