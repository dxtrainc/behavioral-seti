"""
AN ESTIMATOR BENCH FOR ROW 7 -- built so that "more sensitive" cannot mean "tuned".

THE HAZARD. Improving a detector while looking at the data you intend to search is how
a search manufactures its own result. Every knob turned in the direction of a nicer
number is a degree of freedom that no trials correction accounts for.

THE GUARD. Candidates are ranked on a signal WE DID NOT PUT THERE and whose size comes
from outside this programme: the solar-cycle frequency shift, about 0.4 uHz
peak-to-peak at nu_max, known from decades of helioseismology. You cannot tune your way
to a better activity R^2 without genuinely measuring mode frequencies better. So:

    1. GATE     the 135 uHz comb must still come out. An estimator that loses its
                control is not sharper, only less checked.
    2. SELECT   on the KNOWN activity term -- R^2, and the recovered amplitude against
                the literature. This is the only criterion allowed to choose.
    3. REPORT   the injected 95% recovery amplitude. Reported, never optimised.

Two estimators today already show why the bench is needed rather than nice: a centroid
window of +/-8 uHz scored 0.87x (worse than the baseline) because it straddles the
9 uHz l=0/l=2 small separation, and the same window at +/-4 uHz scores 1.33x. That
swing is physics, and the activity R^2 tracked it -- 0.249 against 0.775 -- which is
exactly the behaviour that makes R^2 trustworthy as a selector.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import golf_row7 as G, golf_peakbag as PB, golf_search as SS

NU_MAX = 3.09e-3
LIT_PP = 0.4e-6          # literature solar-cycle shift, p-p at nu_max

def shift_centroid(X, seg_d, half_uHz):
    nss = int(seg_d*86400/G.DT)
    fb, S, duty = PB.seg_spectra(X, nss)
    good = [s for s in S if s is not None]
    if len(good) < 8: return None, None
    mean_P = np.mean(good, axis=0)
    modes = PB.find_modes(fb, mean_P)
    df = fb[1]-fb[0]; half = max(2, int(round(half_uHz*1e-6/df)))
    ref = {}
    for j in modes:
        c, w = PB.centroid(fb, mean_P, j, half)
        if np.isfinite(c): ref[j] = c
    sh = np.full(len(S), np.nan)
    for k, P in enumerate(S):
        if P is None: continue
        d, w = [], []
        for j, c0 in ref.items():
            c, ww = PB.centroid(fb, P, j, half)
            if np.isfinite(c) and ww > 0: d.append(c-c0); w.append(ww)
        if len(d) < 8: continue
        d, w = np.array(d), np.array(w); m = np.abs(d) < 5e-6
        if m.sum() < 8: continue
        sh[k] = np.sum(w[m]*d[m])/np.sum(w[m])
    return sh, seg_d

def score(sh, seg_d, label):
    """gate is checked separately; here: the KNOWN activity term only"""
    act = SS.activity_on_grid(len(sh), seg_d)
    m = np.isfinite(sh) & np.isfinite(act)
    if m.sum() < 12: return None
    A = np.vstack([np.ones(m.sum()), act[m]]).T
    co, *_ = np.linalg.lstsq(A, sh[m], rcond=None)
    fit = A@co; r = sh[m]-fit
    R2 = 1.0 - np.var(r)/np.var(sh[m])
    amp_pp = (fit.max()-fit.min())
    return dict(label=label, seg_d=seg_d, n=int(m.sum()), R2=float(R2),
                amp_pp_uHz=float(amp_pp*1e6),
                amp_ratio=float(amp_pp/LIT_PP),
                resid_rms_uHz=float(r.std()*1e6))

if __name__ == "__main__":
    t0 = time.time()
    X = {t: G.load(t) for t in ("MEAN", "PM1", "PM2")}
    rows = []

    print("CANDIDATES -- selected on the KNOWN solar-cycle term, not on sensitivity\n")
    for seg in (45, 90, 180):
        for hw in (2.0, 4.0, 6.0):
            sh, sd = shift_centroid(X["MEAN"], seg, hw)
            if sh is None: continue
            s = score(sh, sd, "centroid %gd +/-%.0fuHz" % (seg, hw))
            if s: rows.append((s, sh, sd))
    # PM1+PM2 averaged, at the best segment length found above
    best_seg = max(rows, key=lambda r: r[0]["R2"])[0]["seg_d"]
    s1, _ = shift_centroid(X["PM1"], best_seg, 4.0)
    s2, _ = shift_centroid(X["PM2"], best_seg, 4.0)
    if s1 is not None and s2 is not None:
        sh = np.nanmean(np.vstack([s1, s2]), axis=0)
        s = score(sh, best_seg, "PM1+PM2 mean %gd +/-4uHz" % best_seg)
        if s: rows.append((s, sh, best_seg))

    rows.sort(key=lambda r: -r[0]["R2"])
    print("  %-28s %5s %6s %9s %9s %10s" % ("estimator", "n", "R2", "amp p-p", "vs lit", "resid rms"))
    for s, _, _ in rows:
        print("  %-28s %5d %6.3f %8.3f  %8.2f  %9.4f"
              % (s["label"], s["n"], s["R2"], s["amp_pp_uHz"], s["amp_ratio"], s["resid_rms_uHz"]))
    print("\n  baseline (band cross-correlation, 90 d): R2 = 0.517, resid rms 0.0760 uHz")
    best = rows[0][0]
    print("\n  SELECTED on activity R^2: %s" % best["label"])
    print("    R^2 %.3f, recovered amplitude %.3f uHz = %.2fx the literature 0.4 uHz"
          % (best["R2"], best["amp_pp_uHz"], best["amp_ratio"]))
    print("    residual rms %.4f uHz -> implied limit %.2e fractional at nu_max"
          % (best["resid_rms_uHz"], 1.94e-5*best["resid_rms_uHz"]/0.0760))
    print("\n  NOTE: the implied limit is an extrapolation from the baseline's injected")
    print("  figure. It is NOT a limit until injected directly -- that is step 3.")
    json.dump([r[0] for r in rows], open(os.path.expanduser("~/golf_bench.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))
