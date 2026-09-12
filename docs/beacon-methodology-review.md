# Methodology review request — where my analysis kept going wrong

**For Fable.** Written 2026-09-12 after a day of running §3.6's unsearched rows. Everything
below is from that one session. I am asking for a review of **method**, not of results; the
results are in `beacon-paper-draft-final.md` and the repository at commit `1c8005f`.

I am writing this because the errors formed two clear patterns rather than being scattered,
and I do not trust my own read on how much they should worry us.

---

## 1. The two patterns

### Pattern A — every control that failed was a control that could not fail

Seven controls failed today. **Not one was a well-designed control rejecting a good result.**
All seven were controls that were incapable of returning "no":

| # | control | how it could not fail | what it nearly certified |
|---|---|---|---|
| 1 | row 5 transit gate v1 | pass threshold set at 0.4× the predicted depth, i.e. 4.9 ppm for Mercury, against a residual rms of 62 ppm | "6 of 6 transits recovered" — all six numbers were noise |
| 2 | gap control (`wl_cal2.py`) | synthetic data phase-randomised on the compressed index axis, so it had no gap-jumps to mask | a pipeline over-firing **9.39×**, 34/45 tests failing uniformity |
| 3 | GOLF continuum filter | 501 bins = 0.61 µHz, narrower than a mode linewidth | comb ACF of r = 0.004 — no comb, in the instrument p-modes are measured with |
| 4 | VIRGO envelope filter | 2001 bins = 2.3 µHz against a ~1000 µHz envelope | gate 1 failing on 2 of 3 channels, estimator railing at the band edge |
| 5 | centroid window | ±8 µHz, straddling the 9 µHz l=0/l=2 small separation | peak-bagging scoring **worse** than the method it replaced (0.87×) |
| 6 | row 7 injection scan | no zero-amplitude arm | **100% "recovery" at zero injected signal**; a reported 12× sensitivity gain and margin 0.618 |
| 7 | line veto | notched ±0.03% around 180 s; the real sidebands reach ±2% | twenty instrument lines presented as candidates |

**Every control that WORKED was one that did not depend on my reasoning:**

- the EVE ESP **dark diode** — same electronics, no photons
- the **cross-product gate** between the PMO6 radiometer and the three SPM photometers, which caught all twenty candidates the line veto missed (R ≈ 2000 vs R < 5.1)
- the **planetary transit depth**, (R_planet/R_☉)², two radii and no fit
- the **literature amplitude anchor** for the solar-cycle frequency shift

The distinction is sharp and I did not design for it: controls anchored in an external fact
worked; controls anchored in my own analysis of what could go wrong did not.

### Pattern B — five predictions, all wrong in the flattering direction

| prediction | reasoned | measured |
|---|---|---|
| peak-bagging would improve row 7 | ~6× | **1.33×** |
| long-period blindness was VIRGO's 2-month L2 highpass | archive's fault | **my own 1-day detrend**; the archive's filter is never reached |
| a 45-day detrend would open the 1 d – 2 month band | would open it | **did not**; zero recovery at 10 d+ at any amplitude |
| §3.6's margin column (rows 5, 6, 7) | ~1, ~1, 0.1–1 | **unreachable, 0.5, 0.069** |
| the twin-instrument verdict was wrong because of index-axis surrogates | surrogates to blame | **arms A and B unchanged** under correct surrogates — my rejection may itself be unfounded |

Five for five. **None was wrong in the pessimistic direction.** That asymmetry is what I would
most like reviewed: a random error process should miss both ways.

---

## 2. What I think is happening, and where I am unsure

**My working diagnosis.** When I reason about what could go wrong with an analysis, I am
reasoning inside the same frame that produced the analysis, so I systematically fail to
anticipate the failure modes that frame does not contain. A specification says what an
instrument resolves; it does not say what survives a detrending, a continuum estimate and a
trials correction — and all four of my optimistic margin estimates came from specifications.

**A concrete instance worth checking.** Four of the seven control failures are literally the
same mistake: **a smoothing or masking window chosen without checking it against the scale of
the feature it operates on.** 501 bins reads as "wide" until you notice it is 0.61 µHz against
a mode linewidth. 2001 bins reads as wide until it is 2.3 µHz against a 1000 µHz envelope.
±8 µHz reads as narrow until it crosses a 9 µHz mode spacing. This should be a lint rule, not
a judgement call, and I do not know what the general form of that rule is.

**What I am unsure about.** Whether the corrections I applied are themselves sound, or merely
the next layer of the same error. I fixed the window sizes by measuring against the feature —
but I chose which feature to measure against, and that is a judgement of the same kind.

---

## 3. What I changed, and what I would like checked

**Adopted during the session:**

1. **Every injection scan now carries a zero-amplitude arm** and a monotonicity check, and the
   script refuses to print a limit if either fails.
2. **Estimators are selected on a known external signal, never on sensitivity.** The row 7
   estimator bench ranked ten candidates on how well they recover the solar-cycle frequency
   shift against its literature amplitude, and reports the injected limit only afterwards.
3. **A proposed calibration standard** (Part IV §17 of `beacon-note-gate-as-cipher.md`):
   injection validates *recovery*; an amplitude anchor validates *calibration*; they catch
   different failures, and injection alone cannot catch attenuation upstream of the injection
   point. Row 5 is the worked example — an injection into the detrended residual would have
   recovered perfectly and certified a search that cannot see a 76 ppm transit.
4. **Veto costs are reported.** The widened line veto excludes 12% of the band and says so; a
   notch broad enough to remove every candidate is not a null unless its cost is stated.

**Four questions I would value an outside answer to:**

- **Is the "controls that cannot fail" pattern a design problem or a review problem?** Every one
  was caught eventually, by running an arm I had not planned to run. Should there be a standing
  requirement that every control ship with a demonstration of it failing?
- **How much should the five-for-five optimistic bias discount the unmeasured claims still in
  the paper?** §3.6's margin column is now flagged, but the same reasoning produced other
  numbers — the energetics in §2.6, the compute estimates in §5.2 — and those have not been
  measured against anything.
- **Is the estimator bench sound, or is it laundering a selection?** Candidates are ranked on
  the solar-cycle term and the limit measured afterwards, which I believe separates selection
  from sensitivity. But both come from the same data, and I would like that checked by someone
  who did not design it.
- **Four artefact families appeared and three were unanticipated** — instrument-line sidebands
  at 60× the predicted width, solar-cycle harmonics, harmonics of my own detrend window, and
  the Nyquist frequency. Is there a principled enumeration of artefact families for a search of
  this kind, or is discovering them one at a time the normal state of affairs?

---

## 3b. A sweep of the whole codebase, and what it found

**The sharper diagnosis, after review.** §1 framed this as incomplete reasoning about failure
modes. That is too generous. **Four of the seven failures were the same error, and three of
them were made *after* the first had been found and fixed.** Each was corrected where it
surfaced and never swept for elsewhere. The same is true of the zero-amplitude arm — row 5's
gate taught "a control must be able to fail", and an injection scan without a zero arm was
written hours later — and of the calendar-symmetric surrogate, fixed in `wl_cal4.py` and never
propagated to `wl_twin.py`, which is why that verdict is now in doubt.

**The defect is not subtlety. It is no propagation.** That is mechanically fixable, so
`validate/lint_controls.py` now sweeps every search and validation script for the five failure
modes that actually occurred. Run against 68 scripts:

| code | sites | what it flags |
|---|---|---|
| **W** | 26 | a smoothing/masking window with no stated feature scale within eight lines |
| **S** | 11 | phase randomisation on the compressed axis with no full-grid/subsample path |
| Z | 0 | injection scan with no zero-amplitude arm |
| T | 0 | pass threshold set as a fraction of the expected value |
| V | 0 | frequency-veto width asserted rather than measured |

**These are 37 sites to check, not 37 bugs.** The lint is heuristic: a window of 801 bins may be
correctly sized, and the check only reports that nothing nearby records what it was sized
against. An index-axis surrogate is *correct* on an evenly sampled series and fails only on
gapped data. But that absence of recorded justification is precisely the state in which today's
four window errors were made.

**Three [S] hits matter more than the rest, because they sit under published results:**

- `wl_sweep.py` — the quadruple sweep, on a 9,786-day series with 1,064 breaks, which is the
  exact configuration measured at **9.39× over-firing**. The conclusion survives: z = 45.3 was
  re-verified with calendar-symmetric surrogates in `wl_cal4.py`. But the committed script
  would reproduce the result through the wrong null.
- `viewpoint_fast.py` and `viewpoint_search.py` — these produced §4.12's Earth-vs-Mars result on
  MAVEN data with 94.3% coverage, so those series are gapped as well. **§4.12 carries a stated
  limit of 8.9×10⁻⁴ at 3 d and a geometry recovery at _p_ = 5×10⁻⁴. If its null is miscalibrated
  in the same way, that limit moves.** This is the one place in the paper where a published
  number could change, and it has not yet been checked.

**What I would value from a statistics review**, beyond the four questions in §3: whether the
five lint checks are the right five, and what a principled general form of check W would be. I
arrived at "compare the window to the feature scale" by generalising from four instances, which
is the same inductive move that produced the errors in the first place.

## 4. What this does not undermine

Stated plainly so the document is not read as more alarming than it is.

**No result was published in error.** Every failure above was caught before it reached the
manuscript, and several were caught only because an arm was run that did not have to be. The
paper's current numbers survive an audit of 22 of 22 claims against committed result files.

**The failures were asymmetric in a way that favours the conclusions.** All seven controls that
could not fail would have produced **false positives or false limits** — "6 of 6 transits
recovered", a 12× sensitivity gain, twenty candidates. Every one, on correction, moved the
result *toward* a null and *away* from a detection. Nothing was suppressed; things were
prevented from being invented.

**The day's actual output is four searched rows, none of them a detection**, two bands measured
as unsearchable rather than unsearched, and a limit at margin 0.5 which is the closest this
programme has come to the level its own logic predicts. That stands.
