"""
CAN rank3 ACTUALLY FIRE? The other half of the validation.

rank3 passes most of the phase-0 gate: it does not fire on triples built to
contain no three-way structure. That is necessary and it is the half that is
easy to achieve BY ACCIDENT -- a statistic that can never fire would pass phase 0
perfectly. This measures the other half.

THE INJECTED SIGNAL MUST BE THREE-WAY AND NOT PAIRWISE, or the test is trivial:

    a *= (1 + A*u)      u, v independent, zero-mean, +/-1 in ~10-day blocks
    b *= (1 + A*v)
    c *= (1 + A*u*v)

Pairwise terms carry E[uv] = 0, so no pair acquires any correlation; the
three-way term carries E[u^2 v^2] = 1, so the co-moment acquires A^3. This is
the canonical XOR-like dependence that is invisible to every pairwise statistic
and is exactly what the triple search exists to find. The blocks sit well inside
the 181-day detrend window so the injection survives prep().

THE FIRST AMPLITUDE IS ZERO. Row 7's original injection reported 100% recovery
at its smallest amplitude because it had no zero arm and was measuring its own
residual; that mistake is not repeated.

The pairwise correlations are measured at every amplitude and reported, so the
claim that the injection is three-way-only is checked rather than asserted.
"""
import argparse, json, os, sys
import numpy as np
from scipy import stats as sps

sys.path.insert(0, "/home/dxtra/beacon-repo/validate")
from phase0 import prep, make_triple                       # noqa: E402
from tripsurr import Plan                                  # noqa: E402
from d_bed import STATS                                    # noqa: E402

BLK = 10          # days a sign is held; well inside the 181-day detrend window


def pm1(n, rng, blk=BLK):
    nb = int(np.ceil(n / blk))
    return np.repeat(rng.choice((-1.0, 1.0), size=nb), blk)[:n]


def inject(tri, amp, rng):
    n = len(tri[0])
    u, v = pm1(n, rng), pm1(n, rng)
    return [tri[0] * (1 + amp * u), tri[1] * (1 + amp * v), tri[2] * (1 + amp * u * v)]


def one(tri, amp, nsh, iters, rng, sf):
    t = inject(tri, amp, rng) if amp > 0 else [x.copy() for x in tri]
    D = [prep(x)[1] for x in t]
    obs = sf(D)
    if not np.isfinite(obs):
        return None
    plan = Plan(np.vstack(t))
    null = []
    for _ in range(nsh):
        Q = [prep(x)[1] for x in plan.draw(rng, iters)]
        v = sf(Q)
        if np.isfinite(v):
            null.append(v)
    if len(null) < nsh // 4:
        return None
    a = np.asarray(null)
    ne = int((np.abs(a - a.mean()) >= abs(obs - a.mean())).sum())
    # pairwise correlations, to verify the injection added none
    pc = [float(np.corrcoef(D[i], D[j])[0, 1]) for i, j in ((0, 1), (0, 2), (1, 2))]
    return dict(p=(1.0 + ne) / (1.0 + a.size), amp=amp, pc=pc,
                z=float((obs - a.mean()) / max(a.std(), 1e-12)))


CFG = {}


def _work(job):
    t, chans, n, amp = job
    r = np.random.default_rng(5000 + t * 31 + int(amp * 1e7))
    tri = make_triple(chans, n, CFG["ltrend"], CFG["lres"], np.random.default_rng(1000 + t))
    return one(tri, amp, CFG["shifts"], CFG["iters"], r,
               STATS[CFG.get("stat", "rank3")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntrip", type=int, default=32)
    ap.add_argument("--shifts", type=int, default=2000)
    ap.add_argument("--iters", type=int, default=6)
    ap.add_argument("--ltrend", type=int, default=1460)
    ap.add_argument("--lres", type=int, default=180)
    ap.add_argument("--pool", type=int, default=64)
    ap.add_argument("--stat", default="rank3", choices=list(STATS))
    ap.add_argument("--out", default="~/d_inject.json")
    A = ap.parse_args()

    z = np.load(os.path.expanduser("~/beacon_channels.npz"), allow_pickle=True)
    M = z["M"]; fin = [np.isfinite(M[i]) for i in range(M.shape[0])]
    usable = [i for i in range(M.shape[0]) if fin[i].sum() >= 4600]
    rng = np.random.default_rng(20260915)

    base = []
    for t in range(A.ntrip):
        pick = rng.choice(usable, size=3, replace=False)
        n = 4500
        base.append((t, [M[i][fin[i]][:n] for i in pick], n))

    amps = [0.0, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
    global CFG
    CFG = dict(ltrend=A.ltrend, lres=A.lres, shifts=A.shifts,
               iters=A.iters, stat=A.stat)

    print("\n%s INJECTION TEST -- three-way signal, no pairwise signature" % A.stat)
    print("  %d triples per amplitude, %d draws each\n" % (A.ntrip, A.shifts), flush=True)
    print("  %10s %9s %9s %10s   %s" % ("amplitude", "p<0.05", "p<0.001", "median z", "max |d pairwise r|"))
    print("  " + "-" * 68)

    from multiprocessing import Pool
    curve, pc0 = [], None
    for amp in amps:
        jobs = [(t, c, n, amp) for (t, c, n) in base]
        with Pool(A.pool) as p:
            R = [x for x in p.map(_work, jobs) if x]
        ps = np.array([r["p"] for r in R]); zs = np.array([r["z"] for r in R])
        pcs = np.array([r["pc"] for r in R])
        if pc0 is None:
            pc0 = pcs.mean(axis=0)
        dpc = np.abs(pcs.mean(axis=0) - pc0).max()
        rec = float((ps < 0.05).mean())
        curve.append(dict(amp=amp, rec05=rec, rec001=float((ps < 0.001).mean()),
                          med_z=float(np.median(zs)), dpc=float(dpc)))
        print("  %10.1e %8.1f%% %8.1f%% %10.2f   %.4f"
              % (amp, 100 * rec, 100 * (ps < 0.001).mean(), np.median(zs), dpc), flush=True)

    r05 = [c["rec05"] for c in curve]
    mono = all(r05[i] >= r05[i - 1] - 0.10 for i in range(1, len(r05)))
    print("\n  zero-amplitude recovery: %.1f%%  %s"
          % (100 * r05[0], "(correct -- does not fire on nothing)" if r05[0] <= 0.10
             else "*** FIRES ON NOTHING ***"))
    print("  monotonic in amplitude : %s" % mono)
    a95 = None
    for i in range(1, len(r05)):
        if r05[i] >= 0.95 > r05[i - 1]:
            a95 = amps[i]; break
    print("  95%% recovery at        : %s" % ("%.1e" % a95 if a95 else "not reached on this grid"))
    print("  max pairwise-r drift   : %.4f  (injection is three-way only)"
          % max(c["dpc"] for c in curve))
    ok = r05[0] <= 0.10 and mono and r05[-1] > 0.5
    print("\n  VERDICT: %s" % ("rank3 CAN FIRE -- usable as a detector" if ok
                               else "rank3 does not recover an injected signal"))
    json.dump(dict(curve=curve, zero_arm=r05[0], monotonic=bool(mono), a95=a95),
              open(os.path.expanduser(A.out), "w"), indent=1)
    print("  saved %s" % A.out)


if __name__ == "__main__":
    main()
