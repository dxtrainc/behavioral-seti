"""
SWEEP EVERY SEARCH SCRIPT FOR THE FAILURE MODES ALREADY FOUND.

Seven controls failed in one session and FOUR WERE THE SAME ERROR, made three more
times after the first was found and fixed. The defect is not that the errors are subtle;
it is that each was fixed where it surfaced and never swept for elsewhere. This is that
sweep, and it should be run after any new search script is written.

Each check corresponds to a failure that actually happened and what it cost:

  W  smoothing/masking window vs the scale of the feature it operates on
     GOLF continuum 501 bins = 0.61 uHz against a mode linewidth -> comb ACF r = 0.004.
     VIRGO envelope 2001 bins = 2.3 uHz against a ~1000 uHz envelope -> gate 1 failed.
     centroid +/-8 uHz across a 9 uHz mode separation -> peak-bagging scored 0.87x.

  Z  injection scan with no zero-amplitude arm
     row 7 injection reported 100% recovery at ZERO signal and a 12x sensitivity gain.

  S  index-axis surrogate on a gapped series
     measured at 9.39x over-firing; fixed in wl_cal4.py, NOT propagated to wl_twin.py.

  T  pass threshold defined as a fraction of the expected value rather than against a
     measured null
     row 5 gate v1 reported 6 of 6 transits recovered; all six numbers were noise.

  V  frequency veto whose width is asserted rather than measured
     notched +/-0.03% around 180 s; real sidebands reach +/-2%, sixty times wider.
"""
import os, re, sys, glob

ROOT = os.path.expanduser("~/beacon-repo")
FINDINGS = []

def flag(code, path, line, msg):
    FINDINGS.append((code, os.path.relpath(path, ROOT), line, msg))

def scan(path):
    src = open(path, encoding="utf-8", errors="ignore").read()
    lines = src.split("\n")
    for i, L in enumerate(lines, 1):
        s = L.strip()
        if s.startswith("#") or s.startswith('"'): continue

        # W -- a literal window size with no nearby statement of the feature scale
        m = re.search(r"(median_filter|uniform_filter|convolve)\s*\([^)]*size\s*=\s*(\d{2,})", L) \
            or re.search(r"\bsize\s*=\s*(\d{3,})\b", L) \
            or re.search(r"\bw(?:in|min)?\s*=\s*(\d{3,})\b", L)
        if m and "DETREND_MIN" not in L:
            ctx = "\n".join(lines[max(0, i-8):i+3])
            if not re.search(r"uHz|linewidth|separation|envelope|feature|corner|Nyquist|scale", ctx):
                flag("W", path, i, "window %s with no stated feature scale within 8 lines" % m.groups()[-1])

        # Z -- amplitude list that does not start at zero
        m = re.search(r"(AMPS|EPS|amps|eps)\s*=\s*[\(\[]\s*([0-9.]+)", L)
        if m and float(m.group(2)) != 0.0:
            flag("Z", path, i, "injection amplitudes start at %s, not 0 -- no control arm" % m.group(2))

        # T -- threshold as a fraction of the expected value
        if re.search(r"[<>]\s*[0-9.]+\s*\*\s*(depth|expected|predicted|exp_)", L):
            flag("T", path, i, "pass threshold is a fraction of the expected value, not a measured null")

        # V -- asserted veto width
        m = re.search(r"halfwidth\w*\s*=\s*([0-9.eE-]+)", L)
        if m:
            ctx = "\n".join(lines[max(0, i-12):i+3])
            if not re.search(r"measured|observed|sideband", ctx):
                flag("V", path, i, "veto half-width %s asserted, not measured against observed sidebands" % m.group(1))

    # S -- index-axis surrogate on a series known to be gapped
    if re.search(r"np\.fft\.rfft\(", src) and re.search(r"uniform\(0,\s*2\*np\.pi", src):
        if not re.search(r"SUB|full grid|calendar|subsample", src):
            flag("S", os.path.join(ROOT, os.path.basename(path)), 0,
                 "phase-randomises on the compressed axis; no full-grid/subsample path")

if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(ROOT, "search", "*.py")) +
                   glob.glob(os.path.join(ROOT, "validate", "*.py")))
    for f in files:
        if os.path.islink(f): continue
        scan(f)
    print("swept %d scripts\n" % len(files))
    by = {}
    for c, p, l, m in FINDINGS: by.setdefault(c, []).append((p, l, m))
    names = {"W": "window vs feature scale", "Z": "no zero-amplitude arm",
             "S": "index-axis surrogate on gapped data",
             "T": "threshold as a fraction of the expected value",
             "V": "asserted veto width"}
    for c in "WZSTV":
        rows = by.get(c, [])
        print("[%s] %-42s %d" % (c, names[c], len(rows)))
        for p, l, m in rows[:14]:
            print("      %-34s %5s  %s" % (p, l or "-", m))
        if len(rows) > 14: print("      ... %d more" % (len(rows)-14))
        print()
    print("total: %d" % len(FINDINGS))
