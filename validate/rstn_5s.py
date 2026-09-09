"""
The ~5.9 s multi-station peak: solar, or a shared receiver design?

Three questions decide it.
 1 band edge? 1 s sampling puts Nyquist at exactly 2.000 s, so the 2.006 and
   2.007 s peaks are rejected outright.
 2 ratio or channel? The carrier is a ratio. A solar modulation of the spectral
   index moves the ratio; a receiver artefact -- gain switching, calibration,
   digitiser -- shows in the CHANNELS.
 3 is the twin gate even valid? RSTN is a network of standardised receivers, so
   the three sites are not independent designs. A per-design artefact appears at
   all three exactly as a solar signal would. This is what defeated the
   GOES-16/17 gate in section 4.4.
"""
import sys, glob, os
import numpy as np
sys.path.insert(0, "/home/dxtra")
import rstn_search as R

TARGETS = [5.900, 5.910, 5.555, 2.975, 2.721]

def pool(station, k=None, ratio=None, nmax=36):
    acc = None; f = None; n = 0
    for p in sorted(glob.glob("/home/dxtra/rstn/%s_*.gz" % station))[:nmax]:
        g = R.read_day(p)
        if g is None: continue
        t, v = g
        if len(t) < R.NSAMP: continue
        v = v[:R.NSAMP]
        if ratio is None:
            x = v[:, k]
        else:
            a, b = ratio
            x = v[:, a]/np.maximum(v[:, b], 1e-9)
        y = np.log(np.clip(x, 1e-9, None))
        m = np.isfinite(y)
        if m.sum() < 0.5*len(y): continue
        y = np.where(m, y, np.median(y[m]))
        if np.std(y) == 0: continue
        r = R.detrend(y, 601)
        ff, P = R.spec(r)
        if not np.isfinite(P).all(): continue
        acc = P if acc is None else acc + P
        f = ff; n += 1
    if acc is None or n == 0: return None
    P = acc/n
    C = R.cont(P); Rr = P/C
    inb = (1.0/f >= 2.2) & (1.0/f <= 300.0)     # 2.2 s: clear of Nyquist
    thr, mu = R.thr_of(Rr[inb], ndays=n)
    return f, Rr, thr, n

def show(tag, got):
    if got is None:
        print("  %-26s no usable data" % tag, flush=True); return
    f, Rr, thr, n = got
    cells = []
    for P0 in TARGETS:
        kk = int(np.argmin(np.abs(f-1.0/P0)))
        loc = Rr[max(0, kk-3):kk+4].max()
        cells.append("%6.2f%s" % (loc, "*" if loc > thr else " "))
    print("  %-26s n=%2d thr=%5.2f  %s" % (tag, n, thr, " ".join(cells)), flush=True)

print("targets (s):        " + "  ".join("%6.3f" % p for p in TARGETS))
print("* = above that series' own threshold\n")
for st in ("learmonth", "palehua", "san-vito"):
    print("%s:" % st, flush=True)
    days0 = None
    for p in sorted(glob.glob("/home/dxtra/rstn/%s_*.gz" % st))[:1]:
        g = R.read_day(p)
        if g: days0 = g[1]
    good = [k for k in range(8) if days0 is not None
            and np.isfinite(days0[:, k]).mean() > 0.9 and np.nanmedian(days0[:, k]) > 5]
    for k in good:
        show("channel %5d MHz" % R.FREQ[k], pool(st, k=k))
    if len(good) >= 2:
        show("ratio %d/%d" % (R.FREQ[good[0]], R.FREQ[good[-1]]),
             pool(st, ratio=(good[0], good[-1])))
        if len(good) >= 4:
            show("ratio %d/%d" % (R.FREQ[good[1]], R.FREQ[good[-2]]),
                 pool(st, ratio=(good[1], good[-2])))
    print(flush=True)
