"""
THE DICTIONARY ATTACK -- 15 quadruples x 3 four-body forms on the six most
plausible long-record channels, with a MATCHED four-body null.

WHY A DICTIONARY. The key a sender chooses is (which channels, which form). Brute
force over C(30,4) = 27,405 quadruples treats every key as equally likely, but the
sender is constrained: the combination must be one the receiver plausibly MEASURES,
on records long enough to hold a message. That concentrates the key distribution on
the best-instrumented, longest, most obvious observables -- cryptographically, weak
keys. No password cracker starts with brute force; it starts with a wordlist.

THE NULL (section 4.9's requirement, met by construction). One random phase screen
phi(f) applied to ALL channels at once: every pairwise cross-spectrum picks up
exp(i phi) exp(-i phi) = 1 and is preserved EXACTLY (measured: 6.4e-17), while the
bi- and trispectrum pick up a random phase and four-way alignment is destroyed.
Channels are normal-score transformed first so the marginal is preserved too.

THE GAP TRAP. The six-way overlap has 1,064 breaks -- 10.9% of consecutive days are
not consecutive. Compressing to an index series puts a JUMP at every break. Those
jumps carry kurtosis, and the surrogate -- being phase randomised -- spreads that
energy out instead of concentrating it, so the observed statistic would sit above
its own null for a reason that has nothing to do with the Sun. Increments spanning
a break are therefore DROPPED, at identical index positions in data and surrogate.
Leaving them in over-fires by 9.39x (validate/wl_cal3.py).

THE THRESHOLD. The 45 tests are strongly dependent (six channels recur across all
fifteen quadruples), so Bonferroni over-corrects -- the look-elsewhere trials factor
is well below 45. The surrogates give the joint null directly: one phase screen
yields a whole surrogate DATASET, all 45 statistics are computed on it, and the
distribution of the MAXIMUM is the exact FWER threshold under the true correlation
(Westfall-Young max-T). Holm-Bonferroni is reported alongside as the
assumption-free fallback. Sidak is NOT used: it buys power only under independence,
which is exactly what these tests lack.

OUTCOME. Void. See validate/wl_cal4.py for the calibrated null and the note in
docs/: against a null calibrated to 1.02x with 0/45 uniformity failures, real solar
data sits at median z = 45.3. That is intrinsic solar nonlinearity, not a signal --
an omnibus test for four-way structure cannot work on a channel that was never
random. Retained because the null itself is a result.
"""
import io, os, sys, json, time, itertools
import numpy as np
from scipy.stats import norm, kurtosis
from multiprocessing import Pool

def _p(name):
    """Resolve a data or result path: as given, then ./results, then $HOME.
    Same resolver the rest of the repository uses, so this runs from a clone."""
    for c in (name, os.path.join("results", os.path.basename(name)),
              os.path.join(os.path.dirname(__file__), "..", "results", os.path.basename(name)),
              os.path.join(os.path.expanduser("~"), os.path.basename(name))):
        if os.path.exists(c): return c
    return os.path.expanduser(name)

# Channels come from the sweep loader rather than a cached array, so a clone needs
# no pre-built artefact. sweepT.py parses its arguments at import, so only the part
# up to the combination section is executed.
_ST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweepT.py")
_src = io.open(_ST, encoding="utf-8").read()
_g = {"__file__": _ST}
exec(compile(_src[:_src.index("# ---------- dimensionless combinations ----------")],
             "sweepT_prefix", "exec"), _g)
CH = _g["CH"]; names = _g["names"]
M = np.vstack([CH[n] for n in names])
ix = {n: i for i, n in enumerate(names)}

WL = ["TSI", "F10.7", "cosmic ray", "sunspot", "X-ray bg", "MgII"]
OV = np.all(np.isfinite(M[[ix[n] for n in WL]]), axis=0)
DAYS = np.where(OV)[0]
KEEP = np.diff(DAYS) == 1          # increments that do NOT span a break

def normal_score(x):
    r = np.argsort(np.argsort(x)); return norm.ppf((r + 0.5) / len(x))

X = np.vstack([normal_score(M[ix[n]][OV]) for n in WL])
K, N = X.shape

def common_phase(Z, rng):
    F = np.fft.rfft(Z, axis=1)
    ph = rng.uniform(0, 2*np.pi, F.shape[1]); ph[0] = 0.0
    if N % 2 == 0: ph[-1] = 0.0
    return np.fft.irfft(F * np.exp(1j*ph)[None, :], n=N, axis=1)

def stat(v):
    d = np.diff(v)[KEEP]                    # gap-spanning increments dropped
    s = d.std()
    if s <= 0: return np.nan
    return float(kurtosis(d, fisher=True) - 4.0*np.median(np.abs(d))/s)

# 3rd-diff and cross are LINEAR and provably carry nothing a pair cannot see:
# adding c_k*s(t) to channel k gives a four-body pickup of (c_a-3c_b+3c_c-c_d)s and
# a pair pickup of (c_i-c_j)s, so invisibility to every pair forces all c_k equal,
# whereupon (1-3+3-1) = 0. Only the multiplicative form carries a four-way-exclusive
# signal, which lives in the trispectrum. Kept here for the record; the measured
# z-scores corroborate the theorem.
FORMS = {"3rd-diff":  lambda a, b, c, d: a - 3*b + 3*c - d,
         "cross":     lambda a, b, c, d: a - b - c + d,
         "quad-prod": lambda a, b, c, d: a * b * c * d}
QUADS = list(itertools.combinations(range(K), 4))
TESTS = [(q, fn) for q in QUADS for fn in FORMS]

def all_stats(Z):
    return np.array([stat(FORMS[fn](*[Z[i] for i in q])) for q, fn in TESTS])

def sur_stats(seed):
    return all_stats(common_phase(X, np.random.default_rng(seed)))

if __name__ == "__main__":
    NSUR = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    print("channels: %s" % ", ".join(WL))
    print("common overlap %d d (%.1f yr); %d of %d increments kept (%d breaks dropped)"
          % (N, N/365.25, KEEP.sum(), len(KEEP), (~KEEP).sum()))
    print("tests: %d quadruples x %d forms = %d\n" % (len(QUADS), len(FORMS), len(TESTS)))
    obs = all_stats(X)
    t0 = time.time()
    with Pool(80) as p: S = np.array(p.map(sur_stats, range(1, NSUR+1), chunksize=25))
    print("%d surrogates in %.1f s\n" % (NSUR, time.time()-t0))

    pmarg = np.array([(1.0+(S[:, t] >= obs[t]).sum())/(1.0+NSUR) for t in range(len(TESTS))])
    mx = np.nanmax(S, axis=1)
    pmaxT = np.array([(1.0+(mx >= obs[t]).sum())/(1.0+NSUR) for t in range(len(TESTS))])
    z = (obs - np.nanmean(S, 0)) / np.nanstd(S, 0)
    order = np.argsort(pmarg); m = len(TESTS); holm = np.empty(m)
    for r, t in enumerate(order): holm[t] = min(1.0, pmarg[t]*(m-r))
    holm = np.maximum.accumulate(holm[order])[np.argsort(order)]

    print("  %-46s %-10s %8s %9s %9s %9s" % ("quadruple", "form", "z", "p(marg)", "p(maxT)", "p(Holm)"))
    for t in np.argsort(pmarg):
        q, fn = TESTS[t]
        print("  %-46s %-10s %8.1f %9.5f %9.5f %9.5f"
              % (" / ".join(WL[i] for i in q), fn, z[t], pmarg[t], pmaxT[t], holm[t]))
    print("\n  median z of observed against the null: %.1f" % np.nanmedian(z))
    print("  tests with |z| > 20: %d of %d" % ((np.abs(z) > 20).sum(), len(TESTS)))
    print("\n  NOTE: a large z here is NOT a detection. The Sun has intrinsic")
    print("  higher-order structure; see the module docstring and validate/wl_cal4.py.")
    json.dump({"channels": WL, "n_days": int(N), "nsur": NSUR,
               "tests": [{"quad": [WL[i] for i in q], "form": fn, "obs": float(obs[t]),
                          "z": float(z[t]), "p_marg": float(pmarg[t]),
                          "p_maxT": float(pmaxT[t]), "p_holm": float(holm[t])}
                         for t, (q, fn) in enumerate(TESTS)]},
              open(_p("wl_sweep.json"), "w"), indent=1)
    print("\nsaved wl_sweep.json")
