"""
Can the RSTN search fire at all?

Two runs have returned exactly zero peaks at three stations across 84
pair-stations, once with a single-spectrum threshold and once with a Gamma(n)
threshold matched to the n-day average. Zero twice is not obviously a result
about the Sun; it may be a result about the threshold. This asks the question
directly: where does the observed statistic sit relative to its threshold, and
does an injected tone of known depth get recovered?
"""
import os, glob, json
import numpy as np
import rstn_search as R

D = os.path.expanduser("~/rstn")

def main():
    WIN = 601
    paths = sorted(glob.glob(os.path.join(D, "learmonth_*.gz")))[:40]
    acc = None; n = 0; F = None
    days = []
    for p in paths:
        got = R.read_day(p)
        if got is None: continue
        t, v = got
        if len(t) < 20000: continue
        # 2695 / 610 MHz -- both strong at Learmonth
        x = v[:, 4]/np.maximum(v[:, 2], 1e-9)
        m = np.isfinite(x) & (x > 0)
        if m.mean() < 0.9: continue
        y = np.log(np.where(m, x, np.nanmedian(x[m])))
        days.append(y)
    print("days usable: %d, length %d" % (len(days), min(len(d) for d in days)), flush=True)
    L = min(len(d) for d in days)
    days = [d[:L] for d in days]

    def run(inject_amp=0.0, P=97.0):
        acc = None
        for k, y in enumerate(days):
            yy = y.copy()
            if inject_amp:
                t = np.arange(L, dtype=float)
                yy = yy + inject_amp*np.std(R.detrend(y, WIN))*np.cos(2*np.pi*t/P + 0.7*k)
            r = R.detrend(yy, WIN)
            f, Pw = R.spec(r)
            acc = Pw if acc is None else acc + Pw
        Pm = acc/len(days)
        C = R.cont(Pm); Rr = Pm/C
        inb = (1.0/f >= 2.0) & (1.0/f <= WIN/2.0)
        thr, mu = R.thr_of(Rr[inb], ndays=len(days))
        return f, Rr, inb, thr, mu

    f, Rr, inb, thr, mu = run(0.0)
    print("\nno injection:", flush=True)
    print("  mu = %.4f   threshold = %.4f" % (mu, thr), flush=True)
    print("  R in band: median %.4f  max %.4f  99.99pct %.4f"
          % (np.median(Rr[inb]), Rr[inb].max(), np.percentile(Rr[inb], 99.99)), flush=True)
    print("  headroom: max/thr = %.3f  %s"
          % (Rr[inb].max()/thr,
             "<- the statistic never approaches the bar" if Rr[inb].max() < 0.6*thr else ""), flush=True)

    print("\ninjected tone at P = 97 s:", flush=True)
    for a in (0.002, 0.005, 0.01, 0.02, 0.05):
        f2, R2, inb2, thr2, _ = run(a)
        k = int(np.argmin(np.abs(f2-1.0/97.0)))
        loc = R2[max(0, k-2):k+3].max()
        print("  depth %5.3f sigma -> R(97 s) = %8.3f   thr %7.3f   %s"
              % (a, loc, thr2, "DETECTED" if loc > thr2 else "missed"), flush=True)

main()
