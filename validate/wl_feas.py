"""Feasibility probe for a MATCHED higher-order null.

The circular shift voided the triple sweep because it destroys pairwise
structure along with three-way structure (section 4.9). The surrogate we need
holds every pairwise cross-correlation FIXED and destroys only the higher-order
alignment. Schreiber-style constrained annealing does that; the question here is
only whether it converges fast enough to afford 15 keys x N surrogates.
"""
import numpy as np, time, os
M=np.load(os.path.expanduser("~/wl_M.npy"))
names=open(os.path.expanduser("~/wl_names.txt")).read().split("\n")
ix={n:i for i,n in enumerate(names)}
WL=["TSI","F10.7","cosmic ray","sunspot","X-ray bg","MgII"]
print("wordlist:", WL)
for n in WL: assert n in ix, n

# 4-way overlap census over all C(6,4)=15 quadruples
import itertools
print("\n=== 4-way overlap, all 15 quadruples ===")
rows=[]
for q in itertools.combinations(WL,4):
    ov=np.all(np.isfinite(M[[ix[n] for n in q]]),axis=0)
    rows.append((q,int(ov.sum())))
for q,c in sorted(rows,key=lambda r:-r[1]):
    flag="" if c>=2000 else "  <-- BELOW MINOV"
    print("  %-44s %6d d  %5.1f yr%s"%(" / ".join(q),c,c/365.25,flag))
ok=[r for r in rows if r[1]>=2000]
print("\n  %d of 15 quadruples clear the 2000-day gate"%len(ok))
