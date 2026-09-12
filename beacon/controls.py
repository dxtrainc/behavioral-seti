"""
CONTROLS, EACH AS A FUNCTION THAT RETURNS PASS/FAIL AND RECORDS A NEGATIVE ARM.

Rule C: a control ships with a demonstration of it failing, or it does not ship.

Seven controls failed in one session and NOT ONE was a well-designed control rejecting
a good result -- every one was a control that could not return "no". A pass threshold at
0.4x the expected depth against a residual rms fifteen times larger. A synthetic null
with no gaps to mask. An injection scan with no zero-amplitude arm, which reported 100%
recovery at zero signal. A depth-0 control at a single bin against a family-wise
threshold set for 5e5 bins. In each case the control was present, looked reasonable,
and was structurally incapable of failing.
"""
import numpy as np


class ControlResult(dict):
    def __bool__(self):
        return bool(self["passed"])


def negative_arm(statistic, surrogate_fn, threshold, n=200, expected=0.05, tol=3.0):
    """THE arm that must be run: the statistic on data containing no signal.

    `surrogate_fn()` returns one realisation with the signal absent. The observed
    exceedance rate must sit near `expected`; if it is more than `tol` times that, the
    threshold is not calibrated and no limit from it is meaningful.
    """
    hits = sum(1 for _ in range(n) if statistic(surrogate_fn()) >= threshold)
    rate = hits/float(n)
    return ControlResult(name="negative_arm", rate=rate, expected=expected, n=n,
                         passed=bool(rate <= tol*expected),
                         note=("false-alarm rate %.3f against an expected %.3f" % (rate, expected)))


def recovers_known(measured, predicted, tol=0.5, name="known_signal"):
    """An AMPLITUDE anchor: a signal of known size, measured.

    Injection validates recovery; it cannot validate calibration, because a signal
    injected downstream of an attenuating step recovers perfectly while a real one does
    not. Row 5's one-day detrend was absorbing the transit it searched for, and only the
    predicted 75.6 ppm Venus depth could catch that.
    """
    ratio = float(measured)/float(predicted) if predicted else np.nan
    return ControlResult(name=name, measured=measured, predicted=predicted, ratio=ratio,
                         passed=bool(np.isfinite(ratio) and abs(np.log(max(ratio, 1e-12))) <= abs(np.log(1+tol))),
                         note=("recovered %.4g against a predicted %.4g (ratio %.2f)"
                               % (measured, predicted, ratio)))


def monotonic_in_amplitude(curve, slack=0.08):
    """Recovery must rise with injected amplitude. The row 7 scan went
    100, 100, 100, 100, 69.5, 83.5, 100 -- which is what a curve looks like when it is
    measuring the residual rather than the injection."""
    rec = [r for _, r in curve]
    ok = all(rec[i] <= rec[i+1] + slack for i in range(len(rec)-1))
    return ControlResult(name="monotonic", curve=curve, passed=bool(ok),
                         note="recovery %s monotonic in amplitude" % ("is" if ok else "is NOT"))


def independent_of(param_name, values, run_fn, tol=0.25):
    """A feature that tracks an analysis parameter belongs to the analysis.

    Used to test whether a spectral peak moves with the detrend window; the VIRGO
    slow-band peak was checked this way before the injection showed the band had no
    sensitivity at all.
    """
    out = [(v, run_fn(v)) for v in values]
    vals = np.array([o[1] for o in out], float)
    spread = float(np.nanmax(vals)/max(np.nanmin(vals), 1e-30) - 1.0)
    return ControlResult(name="independent_of_"+param_name, results=out, spread=spread,
                         passed=bool(spread <= tol),
                         note=("feature varies by %.0f%% across %s" % (100*spread, param_name)))


def denominator_safe(x, name="channel", margin=5.0):
    """RULE D: a relative residual or ratio is only dimensionless if its denominator
    stays away from zero.

    EVE ESP CH_36 passed every structural test -- right instrument, right cadence, right
    length, correct dimensionless construction -- and failed on the one property nothing
    checked: its values cross zero (median 2.6e-4, minimum -5.2e-4), so x/trend - 1
    diverges and its residual rms came out at 765,522 ppm, i.e. 76%. A ratio formed
    against it is not a carrier, it is a division by something near zero, and any
    survivor in it would be that and nothing else.

    The test is whether the distribution's distance from zero is large compared with its
    own spread: |median| must exceed `margin` robust scatters. A channel that fails this
    may still be searched on its own; it may not be used as a denominator.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return ControlResult(name="denominator_safe:"+name, passed=False,
                             note="too few finite samples to judge")
    med = float(np.median(x))
    mad = float(1.4826*np.median(np.abs(x - med)))
    ratio = abs(med)/mad if mad > 0 else np.inf
    crosses = bool(np.nanmin(x) < 0 < np.nanmax(x))
    ok = (ratio >= margin) and not crosses
    return ControlResult(name="denominator_safe:"+name, median=med, mad=mad,
                         median_over_mad=ratio, crosses_zero=crosses, passed=bool(ok),
                         note=("median %.3g is %.1f robust scatters from zero%s"
                               % (med, ratio, "; RANGE CROSSES ZERO" if crosses else "")))
