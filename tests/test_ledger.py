"""One regression test per ledger entry: each bug that happened, run against the
library function that now prevents it. A fix that is not a test is a fix that will be
made again -- which is the finding these tests exist to close."""
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from beacon import surrogates as S, windows as W, controls as C
from beacon.inject import recovery_curve

def test_gap_report_rejects_the_daily_series():
    """10.9% gaps with 1,064 breaks is where the compressed-axis screen over-fired 9.39x."""
    m = np.ones(16802, bool); rng = np.random.default_rng(0)
    m[rng.choice(16802, 1830, replace=False)] = False
    r = S.gap_report(m, shortest_period_samples=2)
    assert not r["index_axis_ok"], r

def test_gap_report_allows_a_clean_series():
    m = np.ones(100000, bool); m[::5000] = False          # 0.02% gaps, all length 1
    assert S.gap_report(m, shortest_period_samples=125)["index_axis_ok"]

def test_phase_screen_refuses_1d():
    try:
        S.phase_screen(np.zeros(64), np.random.default_rng(0)); assert False, "should raise"
    except ValueError: pass

def test_phase_screen_preserves_pairwise_cross_correlation():
    """6.4e-17 was the measured figure; anything above 1e-10 means the construction broke."""
    rng = np.random.default_rng(1)
    X = np.cumsum(rng.normal(size=(3, 4096)), axis=1)
    X = (X - X.mean(1, keepdims=True))/X.std(1, keepdims=True)
    def xc(Z): return np.array([np.mean(Z[i]*Z[j]) for i in range(3) for j in range(i+1, 3)])
    assert np.abs(xc(S.phase_screen(X, rng)) - xc(X)).max() < 1e-10

def test_roll_masked_reapplies_the_mask():
    ok = np.ones(1000, bool); ok[400:450] = False
    y = np.ones(1000)
    s = S.roll_masked(y, ok, np.random.default_rng(0))
    assert np.all(s[~ok] == 0.0)

def test_checked_rejects_a_detrend_that_eats_its_signal():
    """Row 5: a one-day median against a six-hour transit absorbed the feature."""
    fs = 1/60.0                                   # per minute
    op = W.running_median(1441)
    try:
        W.checked(op, 20000, fs, f_signal=1/(6*3600.0), f_nuisance=1/(30*86400.0), label="row5")
        assert False, "should have rejected the one-day detrend"
    except ValueError: pass

def test_injection_scan_must_start_at_zero():
    try:
        recovery_curve(np.zeros(128), [0.1, 0.2], 10.0, lambda y: 0.0, lambda: np.zeros(128), 1.0)
        assert False, "should raise"
    except ValueError: pass

def test_negative_arm_catches_a_threshold_that_always_fires():
    """The row 7 scan reported 100% recovery at zero injected signal."""
    res = C.negative_arm(lambda y: 10.0, lambda: np.zeros(8), threshold=1.0, n=50)
    assert not res, res

def test_monotonic_catches_the_row7_curve():
    assert not C.monotonic_in_amplitude(
        [(0.0,1.0),(0.005,1.0),(0.02,1.0),(0.03,0.695),(0.045,0.835),(0.07,1.0)])

def test_recovers_known_catches_the_row5_gate():
    """v1 reported Mercury at 53 ppm against a predicted 12.3 and called it a pass."""
    assert not C.recovers_known(measured=53.0, predicted=12.3)
    assert C.recovers_known(measured=90.4, predicted=75.6)

def test_denominator_safe_catches_ch36():
    """CH_36 crosses zero: median 2.6e-4, min -5.2e-4. x/trend - 1 gave rms 76%."""
    rng = np.random.default_rng(3)
    ch36 = rng.normal(2.6e-4, 3.0e-4, 4000)      # crosses zero, as the real channel does
    ch18 = rng.normal(1.23e-3, 2.0e-5, 4000)     # well away from zero
    assert not C.denominator_safe(ch36, "CH_36")
    assert C.denominator_safe(ch18, "CH_18")


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    bad = 0
    for k, v in fns:
        try:
            v(); print("  pass  %s" % k)
        except Exception as e:
            bad += 1; print("  FAIL  %s -- %s" % (k, e))
    print("\n%d of %d passed" % (len(fns)-bad, len(fns)))
    sys.exit(1 if bad else 0)
