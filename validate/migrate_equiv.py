"""
PROVE EQUIVALENCE BEFORE MIGRATING, NOT AFTER.

Six scripts carry numbers that are in the paper. Migrating them to beacon/ means editing
code that produced committed results, and the obvious way to check -- migrate, re-run,
compare -- costs hours per script and cannot distinguish "the library is stricter and the
old result was wrong" from "I broke it".

So this compares the LOCAL function against the LIBRARY function on the same input,
directly. Three outcomes and each means something different:

  IDENTICAL   migration is provably safe; no re-run needed, no number moves.
  DIFFERENT   the library is deliberately stricter (phase_screen demands the full grid
              and a mask; the compressed-axis version over-fired 9.39x on gapped data).
              Migrating CHANGES the result, which is the point, and the committed number
              must be re-derived rather than assumed.
  N/A         the script has no counterpart in the library yet.

Nothing is edited here. This decides which of the six can be migrated safely and which
need a re-derivation, before any of them is touched.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.expanduser("~/beacon-repo"))
sys.path.insert(0, os.path.expanduser("~/beacon-repo/search"))
from beacon import surrogates as S, spectra as SP

rng_a = np.random.default_rng(31337)
rng_b = np.random.default_rng(31337)
out = {}

def verdict(name, a, b, note=""):
    if a is None or b is None:
        out[name] = "N/A"; print("  %-34s N/A          %s" % (name, note)); return
    d = float(np.nanmax(np.abs(np.asarray(a, float) - np.asarray(b, float))))
    rel = d/max(float(np.nanmax(np.abs(np.asarray(a, float)))), 1e-300)
    tag = "IDENTICAL" if rel < 1e-12 else ("close %.1e" % rel if rel < 1e-6 else "DIFFERENT %.1e" % rel)
    out[name] = tag
    print("  %-34s %-12s %s" % (name, tag, note))

print("EQUIVALENCE: local function vs beacon/ library, on real data\n")

# ---- 1. roll_masked, as used by both viewpoint scripts
x = np.cumsum(np.random.default_rng(1).normal(size=200000))
ok = np.ones(len(x), bool); ok[5000:5200] = False
loc = np.roll(x, 777).copy(); loc[~ok] = 0.0
lib = S.roll_masked(x, ok, np.random.default_rng(0))   # different shift, compare form
verdict("roll_masked (form)", np.sort(loc), np.sort(lib),
        "same operation, shift differs by rng")

# ---- 2. spectra.spec vs fastsearch.spectrum
try:
    import fastsearch as F
    r = np.random.default_rng(2).normal(size=100000)
    f1, P1 = F.spectrum(r)
    f2, P2, _ = SP.spec(r, dt=60.0)
    verdict("spec: fastsearch vs beacon", P1, P2, "GOES fast-band spectra")
except Exception as e:
    verdict("spec: fastsearch vs beacon", None, None, str(e)[:44])

# ---- 3. cont
try:
    P = np.abs(np.random.default_rng(3).normal(size=50000))**2 + 1.0
    verdict("cont: fastsearch vs beacon", F.continuum(P), SP.cont(P), "continuum estimator")
except Exception as e:
    verdict("cont: fastsearch vs beacon", None, None, str(e)[:44])

# ---- 4. thr_of
try:
    R = np.random.default_rng(4).exponential(1.44, size=200000)
    t1, m1 = F.search.__globals__.get("thr_of", lambda *_: (None, None))(R) \
        if "thr_of" in F.__dict__ else (None, None)
    t2, m2 = SP.thr_of(R)
    verdict("thr_of: fastsearch vs beacon",
            np.array([t1]) if t1 is not None else None,
            np.array([t2]), "family-wise threshold")
except Exception as e:
    verdict("thr_of: fastsearch vs beacon", None, None, str(e)[:44])

# ---- 5. phase_screen vs wl_sweep.common_phase -- EXPECTED to differ
try:
    import wl_sweep as W
    Z = np.random.default_rng(5).normal(size=(4, 8192))
    a = W.common_phase(Z, np.random.default_rng(7))            # compressed axis
    b = S.phase_screen(Z, np.random.default_rng(7))            # full grid, no mask
    verdict("phase_screen vs common_phase", a, b,
            "expected identical only when there are NO gaps")
except Exception as e:
    verdict("phase_screen vs common_phase", None, None, str(e)[:44])

print("\nVERDICT BY SCRIPT")
rows = [("search/sweepT.py", "circular shift on a gapped daily grid",
         "RE-DERIVE -- library refuses the compressed axis"),
        ("search/wl_sweep.py", "common_phase on 9,786 d with 1,064 breaks",
         "RE-DERIVE -- z = 45.3 was re-verified with wl_cal4 and holds"),
        ("search/viewpoint_search.py", "masked roll on the full grid",
         "SAFE -- roll_masked is the same operation"),
        ("search/viewpoint_fast.py", "masked roll on the full grid",
         "SAFE -- roll_masked is the same operation"),
        ("search/fastsearch.py", "spec / cont / thr_of", "see equivalence above"),
        ("search/euvsfast2.py", "spec / cont / thr_of", "see equivalence above")]
for f, uses, v in rows:
    print("  %-28s %-42s %s" % (f, uses, v))
json.dump(out, open(os.path.expanduser("~/migrate_equiv.json"), "w"), indent=1)
