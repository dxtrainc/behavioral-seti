"""
SENSITIVITY BOUND FOR THE COSMOGENIC ARCHIVES -- COMPUTED BEFORE ACQUIRING ANYTHING.

Three of my sensitivity estimates today came in optimistic against measurement
(peak-bagging 6x -> 1.33x; a 45-day detrend "opening" a band it did not open; three
§3.6 margins). So this is measured on the records already downloaded rather than
reasoned from instrument specifications, and it is measured with beacon.inject, whose
scan must begin at zero.

TWO STATISTICS, because the literature answer differs sharply between them:

  PERIODIC. Fifteen significant periodicities between 40 and 2,320 years are already
  known in these records (Peristykh & Damon 2003), and whether they are real cycles or
  stochastic forcing is disputed (Cameron & Schussler 2019). A periodic search walks
  into a contested null -- the condition that voided three results in this programme.

  IMPULSIVE. The Miyake events were found by looking for sharp features, not periods;
  an FFT missed them for decades. This paper's statistic is a kurtosis of increments,
  built for impulsive and quantised structure, and it does not depend on the disputed
  periodicity question at all. That is the axis worth bounding.

The bound is reported as a FRACTION OF THE MEAN LEVEL -- the units a sender works in
and the units §4.2 quotes -- not as a fraction of the record's own scatter.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.expanduser("~/beacon-repo"))
from beacon.inject import recovery_curve, staircase, impulses
from beacon import surrogates as S

ICE = os.path.expanduser("~/ice")

def load(fn, skip, col, year_col=0):
    yr, v = [], []
    for line in open(os.path.join(ICE, fn), errors="ignore").read().split("\n")[skip:]:
        p = line.split()
        if len(p) <= max(col, year_col): continue
        try:
            y = float(p[year_col]); x = float(p[col])
        except ValueError: continue
        if not np.isfinite(x): continue
        yr.append(y); v.append(x)
    o = np.argsort(yr)
    return np.array(yr)[o], np.array(v)[o]

def stat_kurt(x):
    d = np.diff(x)
    s = d.std()
    if s <= 0: return np.nan
    from scipy.stats import kurtosis
    return float(kurtosis(d, fisher=True) - 4.0*np.median(np.abs(d))/s)

if __name__ == "__main__":
    recs = [("Dome Fuji 10Be conc (Antarctica)", "domefuji-10be2008.txt", 124, 2),
            ("Antarctic 10Be stack (Delaygue)",  "delaygue2010be10.txt",  144, 1)]
    rng = np.random.default_rng(99)
    out = {}
    for name, fn, skip, col in recs:
        yr, v = load(fn, skip, col)
        if len(v) < 50:
            print("  %-34s could not parse" % name); continue
        span = yr.max()-yr.min(); cad = np.median(np.abs(np.diff(yr)))
        rel = v/np.mean(v) - 1.0                       # fractional residual
        print("\n=== %s ===" % name)
        print("  %d points, %.0f yr span, median cadence %.1f yr, Nyquist period %.1f yr"
              % (len(v), span, cad, 2*cad))
        print("  mean %.4g   fractional scatter %.3f   (i.e. %.1f%% point-to-point)"
              % (np.mean(v), rel.std(), 100*np.std(np.diff(rel))/np.sqrt(2)))

        # the surrogate: circular shift, which preserves the record's own spectrum
        # (including the fifteen known periodicities) and destroys only alignment
        ok = np.ones(len(rel), bool)
        def sur():
            return S.roll_masked(rel, ok, rng, margin=5)

        # threshold from the surrogate distribution of the statistic itself
        null = np.array([stat_kurt(sur()) for _ in range(400)])
        null = null[np.isfinite(null)]
        thr = float(np.quantile(null, 0.95))
        obs = stat_kurt(rel)
        print("  impulsive statistic: observed %.3f, threshold %.3f (%.0f pctile)"
              % (obs, thr, 100*(null < obs).mean()))

        # IMPULSIVE injection: a quantised staircase, the structure the statistic is for
        P = max(8*cad, 60.0)
        t = np.arange(len(rel), dtype=float)*cad
        # WAVEFORM MATCHED TO THE STATISTIC. A cosine drove this statistic backwards;
        # the kurtosis of increments is for quantised and impulsive structure.
        for wname, wf in (("staircase", staircase(P, levels=6)),
                          ("impulses", impulses(rate=0.03))):
          res = recovery_curve(rel, [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.80],
                               P, stat_kurt, sur, thr, n_rep=120, rng=rng, dt=cad,
                               waveform=wf)
          print("  --- waveform: %s ---" % wname)
          print("    zero arm %.3f   monotonic %s" % (res["zero_arm"]["rate"], res["monotonic"]["passed"]))
          print("    " + "  ".join("%.2f:%.0f%%" % (a, 100*r) for a, r in res["curve"]))
          lim = res["limit95"]
          print("    --> 95%% recovery at %s (fraction of the mean level)"
                % (("%.3f" % lim) if lim else "NOT REACHED in 0.80"))
          out.setdefault(name, {})[wname] = lim
        out.setdefault(name, {}).update(n=len(v), span_yr=float(span),
                         cadence_yr=float(cad), frac_scatter=float(rel.std()))
    json.dump(out, open(os.path.expanduser("~/ice_bound.json"), "w"), indent=1, default=str)

    print("\n" + "="*70)
    print("COMPARISON, in the units the paper quotes")
    print("  GOES EUVS fast band, injection-verified   1.4e-06   (2 min cadence, 5 yr)")
    print("  VIRGO SPM coherent, injection-verified    2.0e-07   (60 s cadence, 27 yr)")
    print("  neutron monitor, same physical channel    7.6e-05   (1 min cadence, 70 yr)")
    for k, d in out.items():
        for w in ("staircase", "impulses"):
            if d.get(w):
                print("  %-30s %-10s %.1e   (%.0f yr cadence, %.0f yr)"
                      % (k[:30], w, d[w], d["cadence_yr"], d["span_yr"]))
