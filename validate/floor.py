"""
WHAT SETS THE DETECTION FLOOR FOR THREE-WAY-ONLY STRUCTURE?

Four statistics have now been measured and the trade is systematic: the better
calibrated a statistic is, the blinder it is. rank3d passes calibration cleanly
(KS p 0.50) and recovers nothing at 10% modulation; the incumbent kurtosis has
power and a construction artifact 3.85x its own null width. That pattern says
the limit may not belong to any particular estimator.

This measures the limit itself. It is a statement about the INSTRUMENT -- what
this construction, on these archives, at these cadences can reach -- and belongs
in the coverage ledger. It is not a verdict on the framework and does not bear on
whether a gate exists at any amplitude.

THREE CANDIDATE GOVERNORS, which point at different archives:

  estimator noise      floor ~ 1/sqrt(n)   -> long baselines help: Kp to 1932,
                                             decay-rate series, clock ensembles
  natural third-order  floor flat in n     -> length does not help; quieter
     variability                             channels do (in-situ, laboratory)
  cadence              floor tracks the     -> the fast-band searches of 5.5
                       signal timescale

The scan varies n and the injected signal's timescale and finds, per cell, the
amplitude at which recovery reaches 50%. a50 is BRACKETED by the amplitude grid
and reported by interpolation between measured points -- never extrapolated past
the end of the data, which is the practice this project removed from the paper.

Injection is the same three-way-only construction used in the gates: a *= 1+A*u,
b *= 1+A*v, c *= 1+A*u*v with u,v independent, so no pair acquires correlation.
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, "/home/dxtra/beacon-repo/validate")
from phase0 import prep, make_triple                        # noqa: E402
from tripsurr import Plan                                   # noqa: E402
from d_bed import STATS                                     # noqa: E402


def pm1(n, rng, blk):
    nb = int(np.ceil(n / blk))
    return np.repeat(rng.choice((-1.0, 1.0), size=nb), blk)[:n]


def one(job):
    t, chans, n, amp, blk = job
    sf = STATS[CFG["stat"]]
    r = np.random.default_rng(9000 + t * 131 + int(amp * 1e6) + blk)
    tri = make_triple(chans, n, CFG["ltrend"], CFG["lres"], np.random.default_rng(1000 + t))
    if amp > 0:
        u, v = pm1(n, r, blk), pm1(n, r, blk)
        tri = [tri[0] * (1 + amp * u), tri[1] * (1 + amp * v), tri[2] * (1 + amp * u * v)]
    D = [prep(x)[1] for x in tri]
    obs = sf(D)
    if not np.isfinite(obs):
        return None
    plan = Plan(np.vstack(tri))
    null = []
    for _ in range(CFG["shifts"]):
        Q = [prep(x)[1] for x in plan.draw(r, CFG["iters"])]
        w = sf(Q)
        if np.isfinite(w):
            null.append(w)
    if len(null) < CFG["shifts"] // 4:
        return None
    a = np.asarray(null)
    ne = int((np.abs(a - a.mean()) >= abs(obs - a.mean())).sum())
    return (1.0 + ne) / (1.0 + a.size)


CFG = {}


def a50(amps, rec):
    """Interpolate the 50% crossing BETWEEN measured points. None if unbracketed."""
    for i in range(1, len(rec)):
        if rec[i] >= 0.5 > rec[i - 1]:
            f = (0.5 - rec[i - 1]) / max(rec[i] - rec[i - 1], 1e-9)
            return float(np.exp(np.log(amps[i - 1]) + f * (np.log(amps[i]) - np.log(amps[i - 1]))))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntrip", type=int, default=20)
    ap.add_argument("--shifts", type=int, default=1200)
    ap.add_argument("--iters", type=int, default=6)
    ap.add_argument("--ltrend", type=int, default=1460)
    ap.add_argument("--lres", type=int, default=180)
    ap.add_argument("--stat", default="rank3", choices=list(STATS))
    ap.add_argument("--pool", type=int, default=80)
    ap.add_argument("--out", default="~/floor.json")
    A = ap.parse_args()

    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M = z["M"]; fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    global CFG
    CFG = dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts, iters=A.iters, stat=A.stat)

    NS = [2000, 4500, 9000]
    BLKS = [3, 10, 30]
    AMPS = [3e-2, 6e-2, 1e-1, 2e-1, 3e-1]

    print("\nDETECTION FLOOR FOR THREE-WAY-ONLY STRUCTURE  (stat=%s)" % A.stat)
    print("  n in %s, signal block in %s d, %d triples x %d draws per cell\n"
          % (NS, BLKS, A.ntrip, A.shifts), flush=True)
    print("  %6s %7s   %s" % ("n", "block", "  ".join("%7.0e" % a for a in AMPS) + "     a50"))
    print("  " + "-" * 72)

    from multiprocessing import Pool
    out = []
    for n in NS:
        usable = [i for i in range(M.shape[0]) if fin[i].sum() >= n]
        rng = np.random.default_rng(424242)
        base = []
        for t in range(A.ntrip):
            pick = rng.choice(usable, size=3, replace=False)
            base.append((t, [M[i][fin[i]][:n] for i in pick], n))
        for blk in BLKS:
            rec = []
            for amp in AMPS:
                jobs = [(t, c, nn, amp, blk) for (t, c, nn) in base]
                with Pool(A.pool) as p:
                    ps = [x for x in p.map(one, jobs) if x is not None]
                rec.append(float(np.mean(np.array(ps) < 0.05)) if ps else 0.0)
            a = a50(AMPS, rec)
            out.append(dict(n=n, blk=blk, rec=rec, a50=a))
            print("  %6d %7d   %s   %s"
                  % (n, blk, "  ".join("%6.0f%%" % (100 * r) for r in rec),
                     ("%.3f" % a) if a else "  >0.3"), flush=True)

    print("\n  --- what does the floor track? ---")
    for blk in BLKS:
        row = [o for o in out if o["blk"] == blk and o["a50"]]
        if len(row) >= 2:
            ns = np.array([o["n"] for o in row], float)
            aa = np.array([o["a50"] for o in row])
            sl = np.polyfit(np.log(ns), np.log(aa), 1)[0]
            print("  block %2d d: a50 vs n has slope %+.2f   (-0.50 = estimator noise, "
                  "0.00 = natural variability)" % (blk, sl))
        else:
            print("  block %2d d: a50 unbracketed at %d of %d lengths -- floor above 0.3"
                  % (blk, len(NS) - len(row), len(NS)))
    for n in NS:
        row = [o for o in out if o["n"] == n and o["a50"]]
        if len(row) >= 2:
            bs = np.array([o["blk"] for o in row], float)
            aa = np.array([o["a50"] for o in row])
            sl = np.polyfit(np.log(bs), np.log(aa), 1)[0]
            print("  n=%5d  : a50 vs block slope %+.2f   (nonzero = cadence matters)" % (n, sl))
    json.dump(dict(stat=A.stat, ns=NS, blks=BLKS, amps=AMPS, cells=out),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("\n  saved %s" % A.out)


if __name__ == "__main__":
    main()
