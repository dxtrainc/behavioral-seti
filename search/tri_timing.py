"""
WHAT WOULD THE TRIPLE SPACE ACTUALLY COST UNDER THE PAIRWISE-PRESERVING NULL?

The estimate so far was scaled from a PAIR measurement (46.5 ms/draw at
n=16,050, two channels) and spanned 9 to 25 days depending on an unknown: how
long the three-way overlaps actually are. Pairs run to 16,050 days, but
requiring three channels to overlap simultaneously cuts that hard -- demanding
all thirty collapses to 463 days -- and the cost is superlinear in n, so the
overlap distribution is the whole uncertainty.

This measures it instead:

  1. Enumerate every triple and record its REAL overlap. No sampling -- this is
     cheap, and it is the quantity the extrapolation is most sensitive to.
  2. Time actual surrogate draws on a sample of triples spanning that range.
  3. Fit cost = a*n*log2(n) + b, the scaling FFT and argsort both follow, and
     predict the total over every test at its own n rather than at an average.

A cost model fitted at one n and applied at another is how the 9-to-25-day
spread arose in the first place; the per-test prediction is the point.

MUST RUN ON AN UNLOADED BOX. A per-draw time measured under contention is a
measurement of the contention. The launcher waits for the sweeps to finish.
"""
import itertools, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "validate"))
from tripsurr import Plan                                    # noqa: E402

MINOV = 2000
ALPHA = 0.05
ITERS = 12
NSAMP = 14           # triples to time
NDRAW = 25           # draws each


def load_sweepT(path):
    src = open(path, encoding="utf-8").read()
    marker = "# ---- build the job list"
    saved = sys.argv
    sys.argv = ["sweepT.py", "--mode", "triples", "--shifts", "1", "--pool", "1"]
    ns = {"__name__": "sweepT_channels", "__file__": path}
    try:
        exec(compile(src[:src.index(marker)], path, "exec"), ns)
    finally:
        sys.argv = saved
    return ns


S = load_sweepT(os.path.join(HERE, "sweepT.py"))
M, names, N, prep = S["M"], S["names"], S["N"], S["prep"]
fin = [np.isfinite(M[i]) for i in range(N)]

print("enumerating triples (overlap >= %d d)..." % MINOV, flush=True)
tri = []
for i, j, k in itertools.combinations(range(N), 3):
    n = int((fin[i] & fin[j] & fin[k]).sum())
    if n >= MINOV:
        tri.append((i, j, k, n))
ns = np.array([t[3] for t in tri])
# 1 resid-triple + 3 log-2nd-diff per surviving triple, as sweepT builds them
ntest = 4 * len(tri)
thr = ALPHA / ntest
need = int(np.ceil(1.0 / thr))

print("  triples with overlap : %d of %d" % (len(tri), 4060))
print("  tests (4 per triple) : %d" % ntest)
print("  Bonferroni threshold : %.3e   -> draws required M > N/alpha = %d" % (thr, need))
print("\n  overlap distribution (days):")
for q in (0, 10, 25, 50, 75, 90, 100):
    print("    %3d%%  %6d" % (q, int(np.percentile(ns, q))))
print("    mean  %6d" % ns.mean())

# --- time real draws across the observed range ---------------------------
order = np.argsort(ns)
picks = [tri[order[int(x)]] for x in np.linspace(0, len(tri) - 1, NSAMP)]
print("\n  timing %d draws on each of %d triples:" % (NDRAW, len(picks)), flush=True)
rows = []
for (i, j, k, n) in picks:
    ov = fin[i] & fin[j] & fin[k]
    # log-2nd-diff reads the ratio r; resid-triple reads the residual d. Both are
    # three-channel surrogates, so either times the same -- prep once, outside.
    P = [prep(M[c][ov]) for c in (i, j, k)]
    X = np.vstack([p[0] for p in P])
    plan = Plan(X)
    rng = np.random.default_rng(3)
    plan.draw(rng, ITERS)                       # warm
    t0 = time.time()
    for _ in range(NDRAW):
        plan.draw(rng, ITERS)
    ms = (time.time() - t0) / NDRAW * 1e3
    rows.append((n, ms))
    print("    n=%6d  %7.2f ms/draw" % (n, ms), flush=True)

# --- fit cost = a*n*log2(n) + b -----------------------------------------
nn = np.array([r[0] for r in rows], float)
tt = np.array([r[1] for r in rows], float)
x = nn * np.log2(nn)
a, b = np.polyfit(x, tt, 1)
pred = a * x + b
r2 = 1.0 - np.sum((tt - pred) ** 2) / max(np.sum((tt - tt.mean()) ** 2), 1e-12)
print("\n  fit: ms = %.3e * n*log2(n) + %.3f   (R^2 = %.4f)" % (a, b, r2))

# --- predict the whole space, each test at its own n --------------------
cost_ms = a * (ns * np.log2(ns)) + b
core_s = float((cost_ms * 1e-3 * need * 4).sum())      # 4 tests per triple
ch = core_s / 3600.0
print("\n  PREDICTED COST OF THE FULL TRIPLE SPACE")
print("    draws per test : %d" % need)
print("    core-hours     : %.0f" % ch)
for w, lab in ((44, "44 physical cores"), (88, "88 threads")):
    d = ch / w / 24.0
    print("    on %-18s %8.1f h  = %5.1f days" % (lab, ch / w, d))
scr = 25000
ch_s = float((cost_ms * 1e-3 * scr * 4).sum()) / 3600.0
print("\n    a %d-draw SCREEN would be %.0f core-hours = %.1f days on 88 threads"
      % (scr, ch_s, ch_s / 88 / 24))
print("    (a screen cannot fire at threshold; survivors still need the full depth)")

json.dump(dict(triples=len(tri), ntests=ntest, threshold=thr, draws_needed=need,
               overlap_pct={str(q): int(np.percentile(ns, q)) for q in (0,25,50,75,100)},
               fit_a=a, fit_b=b, r2=r2, core_hours=ch, samples=rows),
          open(os.path.expanduser("~/tri_timing.json"), "w"), indent=1)
print("\n  saved ~/tri_timing.json")
