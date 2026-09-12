"""
INJECTION, WITH THE ZERO ARM AND THE UNIT CONVENTION ENFORCED.

THE UNIT CONVENTION IS FRACTIONAL, FULL STOP. A line-of-sight modulator imposes the
same FRACTIONAL amplitude at every viewpoint. viewpoint_fast.py computed the Earth
amplitude as a fraction of EARTH's scatter and injected it as a fraction of MARS's,
inflating the tone by sdM/sdE -- 0.84 against 3.6e-3, a factor of about 230 -- which
overstates Mars power in the direction of declaring Earth-only detections.

THE SCAN STARTS AT ZERO. A scan that starts above zero cannot tell recovery from a
residual that already exceeds the threshold. One did exactly that: 100% "recovery" at
the smallest amplitude tested, which was a 12x sensitivity gain that did not exist.
"""
import numpy as np
from .controls import negative_arm, monotonic_in_amplitude


def tone(period):
    """A coherent sinusoid. Matched to a periodogram or matched-filter statistic."""
    def w(t, rng):
        return np.cos(2*np.pi*t/period + rng.uniform(0, 2*np.pi))
    w.kind = "tone"
    return w


def staircase(period, levels=6):
    """A QUANTISED waveform: the structure a kurtosis-of-increments statistic is for.

    THE WAVEFORM MUST MATCH THE STATISTIC, and this module's first version did not
    offer the choice -- it hardcoded a cosine. Injected into a kurtosis statistic that
    cost nothing at small amplitude and went BACKWARDS at large: a sinusoid's increments
    are sinusoidal, so a bigger tone lowers the kurtosis and raises the scatter, driving
    the statistic away from its threshold. Recovery fell to zero as amplitude rose.
    The paper states the rule for surrogates -- match the null to the statistic -- and
    it applies identically to the injection.
    """
    def w(t, rng):
        ph = rng.uniform(0, 2*np.pi)
        return np.round(levels*np.cos(2*np.pi*t/period + ph))/levels
    w.kind = "staircase"
    return w


def impulses(rate, width=1):
    """Sparse sharp events -- the Miyake-event shape. Also for a kurtosis statistic."""
    def w(t, rng):
        y = np.zeros(len(t))
        k = max(1, int(rate*len(t)))
        for i in rng.choice(len(t)-width, size=k, replace=False):
            y[i:i+width] = 1.0
        return y - y.mean()
    w.kind = "impulses"
    return w


def recovery_curve(series, amplitudes, period, statistic, surrogate_fn, threshold,
                   n_rep=200, rng=None, dt=1.0, pre=None, waveform=None):
    """95%-recovery curve, with the zero arm first and monotonicity checked.

    `amplitudes` are FRACTIONAL and in the units of `series`; they are never rescaled
    by anyone's scatter. `pre` is any processing the real signal would have to survive
    (a detrend, a regression) -- the injection goes in BEFORE it, so the tone meets the
    same filter a real one would rather than being exempted from it.
    """
    rng = rng or np.random.default_rng(0)
    if amplitudes[0] != 0.0:
        raise ValueError("the amplitude scan must begin at 0.0: without a zero arm a "
                         "residual already above threshold reads as perfect recovery")
    t = np.arange(len(series), dtype=float)*dt
    wf = waveform if waveform is not None else tone(period)
    curve = []
    for a in amplitudes:
        hits = 0
        for _ in range(n_rep):
            y = surrogate_fn() + a*wf(t, rng)
            if pre is not None:
                y = pre(y)
            if statistic(y) >= threshold:
                hits += 1
        curve.append((float(a), hits/float(n_rep)))
    zero = negative_arm(statistic, surrogate_fn, threshold, n=n_rep)
    mono = monotonic_in_amplitude(curve)
    lim = None
    if zero and mono:
        for a, r in curve:
            if a > 0 and r >= 0.95:
                lim = a
                break
    return dict(curve=curve, zero_arm=dict(zero), monotonic=dict(mono), limit95=lim,
                waveform=getattr(wf, "kind", "custom"), valid=bool(zero and mono))
