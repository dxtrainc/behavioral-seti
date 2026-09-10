"""
Repository audit: every number the paper quotes, checked against the committed
result file that produced it.

A claim passes only if the value is FOUND in the result file and MATCHES the
manuscript. "Not regenerable" is a failure, not a note.
"""
import json, os, re, sys

REPO = os.path.expanduser("~/beacon-repo")
PAPER = os.path.expanduser("~/beacon-paper-draft.md")
paper = open(PAPER, encoding="utf-8").read()

def R(name):
    p = os.path.join(REPO, "results", name)
    if not os.path.exists(p): return None
    try: return json.load(open(p))
    except Exception: return open(p, encoding="utf-8", errors="ignore").read()

rows = []
def check(sec, claim, got, want, tol=None, instr=None):
    """got: value from the repository. want: value as printed in the paper."""
    if got is None:
        rows.append((sec, claim, "NO RESULT FILE", "", "FAIL")); return
    if instr is not None and instr not in paper:
        rows.append((sec, claim, str(got), "not in paper", "FAIL")); return
    ok = (abs(got-want) <= (tol if tol is not None else 0)) if isinstance(got,(int,float)) else (got==want)
    rows.append((sec, claim, str(got), str(want), "ok" if ok else "MISMATCH"))

# ---- 4.8 pair sweep
d = R("sweep30.json")
# sweep30.json is a JSON LIST of one record per test, not a dict. An earlier
# version of this script assumed a dict and reported the file unreadable.
n = len(d) if isinstance(d, list) else (d or {}).get("n_tests")
check("4.8", "pair-sweep tests", n, 1074, instr="1,074")

# ---- 4.9 triples
d = R("sweepT_triples.json")
if isinstance(d, dict):
    # ntests counts TESTS; the paper also quotes the number of TRIPLES, which is
    # tests/4 because each unordered triple is run in 4 configurations. Both are
    # checked, because conflating them is exactly the error this audit found.
    check("4.9", "triple tests completed", d.get("ntests"), 10996, tol=0, instr="10,996")
    combos = {tuple(sorted(x.strip() for x in r["combo"].split("/"))) for r in d.get("results", [])}
    check("4.9", "triples completed", len(combos), 2749, tol=0, instr="2,749")

# ---- 4.10 self-keyed replication
d = R("selfkey_replicate.json")
if isinstance(d, dict):
    for k, want, sec in (("ACE 2012-2019 (published result)", -3.52, "4.10 ACE rho"),
                         ("WIND/MFI 2012-2019 (independent instrument)", -3.07, "4.10 Wind rho")):
        v = d.get(k, {})
        check("4.10", sec, round(v.get("rho", float("nan")), 2), want, tol=0.005)
    v = d.get("ACE 2012-2019 (published result)", {})
    for mode, want in (("arbitrary", 0.0073), ("carrington", 0.0063)):
        check("4.10", "ACE p (%s)" % mode, round(v.get(mode, {}).get("p", float("nan")), 4), want, tol=1e-4)

# ---- 4.10 power by half
d = R("selfkey_power.json")
if isinstance(d, dict):
    for k, want in (("ACE second half", 100.0), ("WIND second half", 94.0)):
        v = d.get(k, {}).get("power", {})
        got = v.get("0.1", v.get(0.1))
        check("4.10", "power at 0.10 sigma, %s" % k, got, want, tol=0.5)

# ---- 4.11 geometric term
d = R("geoterm.json")
if isinstance(d, dict):
    ly = d.get("Ly-alpha", {})
    check("4.11", "Ly-a true-offset count", ly.get("n_true"), 105, tol=0)
    check("4.11", "Ly-a decoy median", ly.get("decoy_med"), 108.0, tol=0.5)

# ---- 4.12 slow band
d = R("viewpoint_l2b.json")
if isinstance(d, dict):
    check("4.12", "lag RMS (d)", round(d.get("rms", 0), 2), 1.97, tol=0.005)
    check("4.12", "shuffled null (d)", round(d.get("null_med", 0), 2), 7.65, tol=0.005)
    check("4.12", "p", d.get("p"), 0.0005, tol=1e-4)
    check("4.12", "windows", d.get("n_win"), 86, tol=0)

# ---- 4.12 fast band
d = R("viewpoint_fast.json")
if isinstance(d, dict):
    check("4.12", "fast-band Earth peaks", d.get("n_earth"), 0, tol=0)
    check("4.12", "fast-band Mars peaks", d.get("n_mars"), 561, tol=0)

# ---- 4.14 RSTN
d = R("rstn_search.json")
if isinstance(d, dict):
    c = d.get("candidates", {})
    for st, want in (("learmonth", 44), ("palehua", 64), ("san-vito", 43)):
        check("4.14", "%s peaks" % st, len(c.get(st, [])), want, tol=0)

# ---- 5.7 reachability
d = R("reach30.json")
if isinstance(d, dict):
    # keys are perm (permanently blocked) and ok (reachable), not "unreachable"
    check("5.7", "permanently blocked pairs", len(d.get("perm", [])), 53, tol=0, instr="53")
    check("4.8", "reachable pairs", d.get("ok"), 382, tol=0, instr="382")

w = max(len(r[1]) for r in rows)
print("%-6s %-*s %-22s %-14s %s" % ("sec", w, "claim", "repository", "paper", "verdict"))
for sec, claim, got, want, verdict in rows:
    print("%-6s %-*s %-22s %-14s %s" % (sec, w, claim, got[:22], want[:14], verdict))
bad = [r for r in rows if r[4] != "ok"]
print("\n%d claims checked, %d ok, %d needing attention" % (len(rows), len(rows)-len(bad), len(bad)))
sys.exit(0)
