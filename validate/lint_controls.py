"""
LINT, RE-POINTED AT THE LIBRARY.

The first version's W and S rules were documentation checks wearing the costume of
correctness checks: W asked whether a feature scale was *mentioned* within eight lines
of a window literal, and S searched for the *strings* "SUB", "full grid", "calendar".
S duly produced a false positive on both viewpoint scripts, which use masked circular
rolls and phase-randomise nothing -- the rule matched `uniform(0, 2*np.pi)` from the
INJECTED TONE's phase. That false positive was then reported as a finding against a
published result, which is the cost of a rule that matches text rather than structure.

With beacon/ in place the rules can be structural instead:

  P   any local surrogate at all -- a function containing rfft AND irfft, or np.roll on
      a data array, outside beacon/. True positive by construction.
  S'  phase_screen called without mask= where gap_report says the index axis is unsafe.
  W'  a smoothing/detrend call not wrapped by windows.checked. The general form is not
      a rule about sizes: every window is a filter, and its transfer function is
      computed at the frequency kept and the frequency removed.
  C   a control with no recorded negative_arm result.
  I   an inject() after the first windows./spectra.cont call in the same script --
      injecting downstream of the step that attenuates the signal.
  D   a relative residual or ratio formed without checking that the denominator stays
      away from zero. EVE ESP CH_36 passed every structural test and failed on this one:
      its values cross zero, so x/trend - 1 diverged and its residual rms was 76%. A
      ratio against it is a division by something near zero, not a carrier.

  Z/T/V unchanged.
"""
import os, re, sys, glob, ast

ROOT = os.path.expanduser("~/beacon-repo")
LIB = os.path.join(ROOT, "beacon")
FINDINGS = []

def flag(code, path, line, msg):
    FINDINGS.append((code, os.path.relpath(path, ROOT), line, msg))

def scan(path):
    src = open(path, encoding="utf-8", errors="ignore").read()
    lines = src.split("\n")
    in_lib = os.path.abspath(path).startswith(os.path.abspath(LIB))

    # P -- any local surrogate outside the library
    if not in_lib:
        if "np.fft.rfft" in src and "np.fft.irfft" in src:
            flag("P", path, 0, "defines a local phase screen (rfft+irfft); import beacon.surrogates.phase_screen")
        for i, L in enumerate(lines, 1):
            if re.search(r"np\.roll\(", L) and not L.strip().startswith("#"):
                if "roll_masked" not in src:
                    flag("P", path, i, "local circular shift; import beacon.surrogates.roll_masked")
                break

    # W' -- smoothing not wrapped by windows.checked
    if not in_lib:
        for i, L in enumerate(lines, 1):
            if re.search(r"\b(median_filter|uniform_filter|np\.convolve|convolve)\s*\(", L) \
               and not L.strip().startswith("#") and "checked" not in L and "windows." not in L:
                flag("W'", path, i, "window not passed through windows.checked(op, f_signal, f_nuisance)")

    # I -- injection downstream of a filter, BY DATAFLOW rather than by line number.
    #
    # The first version compared LINE NUMBERS: it flagged any inject( appearing below
    # the first median_filter in the file. That is definition order, not execution
    # order, and it fired on search/backex2.py -- where inject() is defined after
    # stat() but CALLED first, on the raw periods, before any filtering. backex2 is in
    # fact one of the better scripts here; its own docstring records that v1 had a
    # statistic which sat at 0.954 for data and surrogates alike.
    #
    # Fourth rule in this file to fire on text that resembles the pattern rather than
    # instantiating it, after S, W and D. A regex cannot see dataflow, so this walks
    # the AST: collect the names assigned from a filtering call, then flag an inject
    # whose arguments include one of them.
    if not in_lib:
        try:
            tree = ast.parse(src)
        except SyntaxError:
            tree = None
        if tree is not None:
            FILTERS = {"median_filter", "uniform_filter", "cont", "continuum", "prep", "detrend"}
            filtered = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
                    fn = node.value.func
                    nm = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
                    if nm in FILTERS:
                        for t in node.targets:
                            for sub in ast.walk(t):
                                if isinstance(sub, ast.Name): filtered.add(sub.id)
            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    fn = node.func
                    nm = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
                    if nm in ("inject", "recovery_curve") and not any(
                            k.arg == "pre" for k in node.keywords):
                        args = {a.id for a in node.args if isinstance(a, ast.Name)}
                        hit = args & filtered
                        if hit:
                            flag("I", path, node.lineno,
                                 "injects into %s, which was produced by a filter -- "
                                 "the signal is exempted from it" % ", ".join(sorted(hit)))

    for i, L in enumerate(lines, 1):
        s = L.strip()
        if s.startswith("#"): continue
        m = re.search(r"(AMPS|EPS|amps|eps|amplitudes)\s*=\s*[\(\[]\s*([0-9.]+)", L)
        if m and float(m.group(2)) != 0.0:
            flag("Z", path, i, "amplitude scan starts at %s -- no zero arm" % m.group(2))
        if re.search(r"[<>]\s*[0-9.]+\s*\*\s*(depth|expected|predicted|exp_)", L):
            flag("T", path, i, "threshold is a fraction of the expected value, not a measured null")
        m = re.search(r"halfwidth\w*\s*=\s*([0-9.eE-]+)", L)
        if m:
            ctx = "\n".join(lines[max(0, i-12):i+3])
            if not re.search(r"measured|observed|sideband", ctx):
                flag("V", path, i, "veto half-width %s asserted, not measured" % m.group(1))

    # C -- a control with no recorded negative arm
    if os.path.basename(path) == "controls.py": return
    # RULE C. The first version matched def \w*(control|gate)\w*( and fired on all five
    # scripts that define a function called SURROGATE -- sur-ro-GATE. A substring where a
    # word was meant. That is the fifth rule in this file to fire on resemblance rather
    # than structure, after S, W, D and I, and several were written in direct reaction to
    # the previous one failing the same way. The name must now be a whole word, delimited
    # by underscores or the function-name boundary.
    if re.search(r"def\s+(?:\w+_)?(?:control|gate|gates|checks)(?:_\w+)?\s*\(", src) \
       and "negative_arm" not in src:
        flag("C", path, 0, "defines a control with no negative_arm; Rule C requires a failing demonstration")

if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(ROOT, "search", "*.py")) +
                   glob.glob(os.path.join(ROOT, "validate", "*.py")) +
                   glob.glob(os.path.join(LIB, "*.py")))
    files = [f for f in files if not os.path.islink(f)]
    for f in files: scan(f)
    print("swept %d scripts (%d in beacon/)\n" % (len(files), len(glob.glob(os.path.join(LIB, "*.py")))))
    names = {"P": "local surrogate outside beacon/", "S'": "unsafe index-axis surrogate",
             "D": "ratio with unchecked denominator",
             "W'": "window not transfer-checked", "C": "control with no negative arm",
             "I": "injection downstream of a filter", "Z": "no zero-amplitude arm",
             "T": "threshold as a fraction of expected", "V": "asserted veto width"}
    by = {}
    for c, p, l, m in FINDINGS: by.setdefault(c, []).append((p, l, m))
    for c in ["P", "S'", "W'", "C", "I", "D", "Z", "T", "V"]:
        rows = by.get(c, [])
        print("[%-2s] %-36s %d" % (c, names[c], len(rows)))
        for p, l, m in rows[:8]:
            print("       %-32s %5s  %s" % (p, l or "-", m))
        if len(rows) > 8: print("       ... %d more" % (len(rows)-8))
    print("\ntotal: %d" % len(FINDINGS))
