"""
RSTN SPECTRAL-INDEX SEARCH AT SECOND CADENCE -- row 2 of section 3.6.

THE CARRIER. Not the flux at any one frequency, which is in solar flux units and
so needs a shared unit, but the RATIO of fluxes at two frequencies, which is
dimensionless (section 3.1). Seven frequencies are populated, giving 21 pairs.
A spectral index is also the natural thing for a redistributor to move: shifting
flux between frequencies changes the ratio while leaving the total alone.

THE BAND. 1-second sampling gives a Nyquist of 2 s. Each day is a contiguous
block of roughly ten daylight hours, detrended locally over a stated window, so
the pass-band runs from 2 s to that window. This is the fastest band anywhere in
the paper -- section 4.2's deepest limit is at two minutes.

CONTROLS, declared before the search:
    the 1-second sampling grid and its harmonics
    the receiver calibration cycle, measured from the data rather than assumed
    integer quantisation -- fluxes are reported as whole SFU, so at 21 SFU one
      count is 4.8%, and the low-frequency channels are quantisation-limited
    solar elevation and atmospheric opacity, both slow and removed by detrending

THE TWIN GATE. Learmonth, San Vito and Palehua are different hardware at widely
separated longitudes. A solar signal is present at all of them; interference, a
calibration cycle or a local artefact is not. This is the gate that worked for
Oulu/Kiel in section 4.4 and failed for GOES-16/17, which share a design.
"""
import os, gzip, glob, json, sys
import numpy as np
from scipy.ndimage import median_filter

D = os.path.expanduser("~/rstn")
FREQ = [245, 410, 610, 1415, 2695, 4995, 8800, 15400]
rng = np.random.default_rng(17)

def read_day(path):
    """-> seconds-of-day, flux array (n, 8)"""
    t, v = [], []
    with gzip.open(path, "rt", errors="ignore") as f:
        for line in f:
            if len(line) < 60: continue
            try:
                hh = int(line[12:14]); mm = int(line[14:16]); ss = int(line[16:18])
            except ValueError:
                continue
            row = []
            ok = True
            for k in range(8):
                fld = line[18+7*k:25+7*k].strip()
                if not fld: row.append(np.nan); continue
                try: row.append(float(fld))
                except ValueError: ok = False; break
            if not ok: continue
            t.append(hh*3600+mm*60+ss); v.append(row)
    if not t: return None
    t = np.asarray(t, float); v = np.asarray(v, float)
    o = np.argsort(t)
    return t[o], v[o]

def detrend(y, win=601):
    """running median: the high-pass corner, and so the long-period edge"""
    return y - median_filter(y, size=win, mode="nearest")

def spec(r):
    n = len(r)//2*2
    w = np.blackman(n)
    P = np.abs(np.fft.rfft((r[:n]-r[:n].mean())*w))**2
    return np.fft.rfftfreq(n, d=1.0)[1:], P[1:]

def cont(P, w=301):
    return np.exp(median_filter(np.log(np.maximum(P, 1e-300)), size=w, mode="nearest"))

def thr_of(R, alpha=0.05, ndays=1):
    """Threshold for the MEAN of `ndays` independent periodograms.

    A single continuum-normalised periodogram is roughly exponential, and the
    familiar mu*ln(N/alpha) follows. This search averages about ninety daily
    spectra, and the mean of n exponentials is Gamma(n, mu/n) -- a far lighter
    tail. Using the single-spectrum threshold on an averaged spectrum sets the
    bar far too high: the first run of this search returned exactly ZERO peaks
    at all three stations across 84 pair-stations, where chance alone should
    give a few, which is the tell. The scale mu is still estimated empirically
    from the data, as everywhere else in this paper."""
    mu = np.median(R)/np.log(2.0)
    if ndays <= 1:
        return mu*np.log(len(R)/alpha), mu
    from scipy.special import gammaincinv
    return (mu/ndays)*gammaincinv(ndays, 1.0 - alpha/len(R)), mu

# Every observing day has a different length -- 25 consecutive Learmonth days
# gave 25 distinct sample counts -- so spectra cannot simply be added. An
# earlier version guarded the accumulation with len(P) == len(acc), which
# accepted the first day of each pair and silently dropped all the rest; every
# pair then failed the n >= 5 test and the search reported zero peaks at three
# stations. Twice. Days are now truncated to a common length instead.
NSAMP = 36000                     # 10 h, the shortest useful observing day

def main():
    WIN = 601                     # s, the detrend window
    PMAX = WIN/2.0
    files = sorted(glob.glob(os.path.join(D, "*.gz")))
    bystn = {}
    for p in files:
        st = os.path.basename(p).split("_")[0]
        bystn.setdefault(st, []).append(p)
    print("stations: %s" % {k: len(v) for k, v in bystn.items()}, flush=True)
    print("detrend %d s -> pass-band 2 s to %.0f s\n" % (WIN, PMAX), flush=True)

    acc = {}
    for st, paths in sorted(bystn.items()):
        pooled = {}
        nday = 0
        for p in paths:
            got = read_day(p)
            if got is None: continue
            t, v = got
            if len(t) < NSAMP: continue
            t, v = t[:NSAMP], v[:NSAMP]
            good = [k for k in range(8) if np.isfinite(v[:, k]).mean() > 0.9
                    and np.nanmedian(v[:, k]) > 5]
            if len(good) < 2: continue
            nday += 1
            for a in range(len(good)):
                for b in range(a+1, len(good)):
                    ka, kb = good[a], good[b]
                    x = v[:, ka]/np.maximum(v[:, kb], 1e-9)
                    m = np.isfinite(x) & (x > 0)
                    if m.mean() < 0.9: continue
                    y = np.log(np.where(m, x, np.nanmedian(x[m])))
                    r = detrend(y, WIN)
                    f, P = spec(r)
                    key = "%d/%d" % (FREQ[ka], FREQ[kb])
                    if key not in pooled: pooled[key] = [np.zeros(len(f)), 0, f]
                    assert len(P) == len(pooled[key][0])   # guaranteed by NSAMP
                    pooled[key][0] += P; pooled[key][1] += 1
        nacc = max((v[1] for v in pooled.values()), default=0)
        print("  %-14s %3d days read, %d pairs, %d days accumulated per pair"
              % (st, nday, len(pooled), nacc), flush=True)
        acc[st] = pooled

    print("\n=== candidates per station (pass-band only) ===", flush=True)
    cands = {}
    for st, pooled in sorted(acc.items()):
        hits = []
        for key, (P, n, f) in sorted(pooled.items()):
            if n < 5: continue
            P = P/n   # incoherent average over n days
            C = cont(P); R = P/C
            inb = (1.0/f >= 2.0) & (1.0/f <= PMAX)
            thr, mu = thr_of(R[inb], ndays=n)
            idx = np.where((R > thr) & inb)[0]
            for i in idx[np.argsort(-R[idx])][:5]:
                hits.append((float(1.0/f[i]), key, float(R[i]), float(thr)))
        hits.sort(key=lambda h: -h[2])
        cands[st] = hits
        print("  %-14s %d peaks above threshold" % (st, len(hits)), flush=True)
        for P0, key, R0, thr in hits[:8]:
            print("      P=%8.3f s  %-12s R=%7.3f (thr %.3f)" % (P0, key, R0, thr), flush=True)

    print("\n=== TWIN GATE: which periods appear at more than one station? ===", flush=True)
    allp = {}
    for st, hits in cands.items():
        for P0, key, R0, thr in hits:
            allp.setdefault(round(P0, 2), set()).add(st)
    multi = {k: v for k, v in allp.items() if len(v) > 1}
    print("  distinct periods: %d;  at >1 station: %d" % (len(allp), len(multi)), flush=True)
    for P0, sts in sorted(multi.items(), key=lambda kv: -len(kv[1]))[:12]:
        print("    P=%8.3f s  at %s" % (P0, ", ".join(sorted(sts))), flush=True)
    if not multi:
        print("  none -- every peak is confined to a single station, which is the"
              " signature of local artefact rather than of the Sun", flush=True)
    json.dump({"window_s": WIN,
               "candidates": {k: v for k, v in cands.items()},
               "multi_station": {str(k): sorted(v) for k, v in multi.items()}},
              open(os.path.expanduser("~/rstn_search.json"), "w"), indent=1)
    print("\nsaved ~/rstn_search.json", flush=True)

if __name__ == "__main__":
    main()
