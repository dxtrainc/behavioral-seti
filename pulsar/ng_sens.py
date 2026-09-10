"""
Sensitivity for the NANOGrav step search, on a ladder that actually brackets it.

The first attempt ran 1e-13 to 1e-11 and returned 100% recovery at every depth,
which says only that the ladder began far above the threshold. A fractional step
dnu/nu produces a residual ramp of dnu/nu * (t-t0), so over a 15 yr baseline
1e-13 is a 24 us signal against a 0.1-3 us rms -- enormous. The scale that
matters is 2*rms/T reduced by the sqrt(N) gain of the matched filter, which is
of order 1e-17.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.expanduser("~"))
import ng_step as NG

# Bracketed by two failed ladders: 1e-13 upward gave 100% everywhere, 1e-15
# downward gave the 5% false-alarm floor everywhere. The limit is between them,
# and it is ~3 orders WORSE than the naive 2*rms/T/sqrt(N) estimate because the
# timing model absorbs most of a ramp -- a linear ramp from t0 is strongly
# degenerate with the fitted nu and nudot. That degeneracy is the real limit
# here, not the measurement noise.
DEPTHS = (1e-15, 3e-15, 1e-14, 3e-14, 1e-13, 3e-13)

def main():
    d = json.load(open(os.path.expanduser("~/ng_step.json")))
    best = sorted(d, key=lambda x: x["rms_us"])[:4]
    print("95% recovery in fractional frequency step, dnu/nu\n", flush=True)
    print("  %-14s %7s %6s  %s" % ("pulsar", "rms/us", "thr", "recovery ladder"), flush=True)
    out = []
    for b in best:
        f = os.path.join(NG.DIR, b["name"]+".pkl")
        if not os.path.exists(f): continue
        name, t, r, e, M = NG.load(f)
        thr, curve, lim = NG.sensitivity(name, t, r, e, M, depths=DEPTHS, N=40)
        print("  %-14s %7.3f %6.0f  %s  ->  95%% at %s" %
              (name, b["rms_us"], thr,
               " ".join("%.0e:%3.0f%%" % c for c in curve),
               ("%.1e" % lim) if lim else ">%.0e" % DEPTHS[-1]), flush=True)
        out.append({"name": name, "rms_us": b["rms_us"], "thr": thr,
                    "curve": curve, "lim95": lim,
                    "span_yr": (t.max()-t.min())/NG.YR, "ntoa": len(t)})
    json.dump(out, open(os.path.expanduser("~/ng_sens.json"), "w"), indent=1)
    print("\nsaved ~/ng_sens.json", flush=True)

main()
