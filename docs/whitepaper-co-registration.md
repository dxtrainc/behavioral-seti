# The co-registration gap in Sun-as-a-star observables

**A white paper on why simultaneity, not sensitivity, is the binding constraint on multi-channel solar time-series analysis — and on the null distribution that would make such analysis defensible**

Robert Griffin · Dxtra Inc. (dxtra.com) · rtg@dxtra.com

---

## Synopsis

Seventy years of solar monitoring has produced a record that is deep in a few
observables, blank in others, and — critically — **assembled from instruments
that never observed together**. Every long solar series is a concatenation across
handovers rather than a measurement, and the size of that effect is measurable:
across the GOES-15 to GOES-16 X-ray handover, detector output shifts by **27.7σ**
(Griffin 2026, search 25).

For any question about a *single* observable this is a manageable
cross-calibration problem. For any question about a **relationship between
observables** it is close to fatal, because the relationship is only measurable
if the quantities are simultaneous, co-registered, and share a calibration chain.

This paper argues that the binding constraint is neither sensitivity nor compute,
and that this is now demonstrable rather than assertable. It asks for three
things in a specific order: a simulation campaign that would give higher-order
solar statistics a calibrated null distribution for the first time; a single
Sun-as-a-star instrument measuring four presently-separated quantities on one
clock; and, on a decadal timescale, two in-situ measurements that no telescope
can supply.

All three are independently motivated by mainstream heliophysics.

---

## 1. The measurable claim

A recent search programme (Griffin 2026) enumerated the space of dimensionless
relationships among 30 solar and heliospheric observables and searched it
exhaustively. The result of interest here is not the search outcome —
twenty-five of twenty-nine searches were null — but what the exercise
**measured about the archive itself**:

- The pair space over 30 observables closes at **88%** — 382 of 435 pairs. The
  remaining 53 are unreachable **because no overlapping record exists**, not
  because of sensitivity.
- Acquiring seventeen further observables moved coverage from **18% to 88%**: a
  direct demonstration that coverage responds to the number of *jointly measured*
  quantities and to nothing else.
- The same instrument, days and pipeline are **245× more sensitive at two minutes
  than at one day**. Sensitivity is not what is limiting.
- Of 4,060 possible triples, 2,934 are reachable and 2,749 were completed. The
  result was **void, not null**, for want of a defensible three-body null.

The last point is addressed in §2, because it must be solved before any new
measurement can be interpreted.

---

## 2. The null distribution, which comes first

The triple and quadruple spaces returned void for a reason that generalises well
beyond this application:

> **Nobody knows the null distribution of higher-order statistics — bispectra,
> trispectra, cross-channel phase alignment — of solar output.**

Observation cannot supply it. There is one Sun and one realisation of it; a null
distribution requires an ensemble. Any published claim of non-linear coupling in
solar time series presently rests on a surrogate whose false-alarm rate has never
been characterised against a Sun-like system.

**Proposed:** a large ensemble of independent global solar convective-dynamo
simulations — Rayleigh (Featherstone & Hindman 2016), ASH (Brun, Miesch &
Toomre 2004) or MURaM (Vögler et al. 2005) class — each integrated over several
simulated activity cycles, from which synthetic disc-integrated irradiance and
velocity series are drawn and higher-order statistics measured across members.

**Scale.** Taking 5×10⁵ to 2×10⁶ core-hours per realisation as an order-of-magnitude
figure, an ensemble of 50–200 members is **0.4 to 1.8 million node-hours** — within
a single INCITE-class award, and modest by the standards of the simulations
themselves. *This bracket is an estimate and is not yet sourced to published
per-run costs; it requires a scoping study against the chosen code before it is
quoted in a proposal.*

**The gate this project must pass, stated first because it decides whether it is
worth running.** A null built from simulations that do not reproduce the Sun *in
the statistics being tested* is worse than no null: it would license exactly the
false positives such a null exists to prevent. The ensemble must be validated
against the real Sun in the same statistics before any limit derived from it
means anything.

**Why it comes first.** New instruments cannot be interpreted without it. A
multi-channel measurement with no characterised false-alarm rate produces void
results exactly as the existing archive does. This is also by far the cheapest of
the three requests, and the only one achievable within a single allocation cycle.

---

## 3. One instrument: co-registration is the requirement

Four quantities sit high in the ordering of what a Sun-as-a-star programme should
measure, and all four are either uninstrumented or short of the precision
required. Stated reaches below are measured by injection-recovery on real
archival data, not design goals:

| Quantity | Required | Presently achieved | Shortfall |
|---|---|---|---|
| Disc-integrated spectropolarimetry | 10⁻⁷ | 10⁻⁵–10⁻⁶ | **10–100×** |
| Line-profile ratios, disc-integrated | 10⁻⁷ | 10⁻⁶ | **10×** |
| p-mode frequency structure | 10⁻⁶ | 1.46×10⁻⁵ (SOHO/GOLF, 25.9 yr) | **15×** |
| Core g-modes | — | no detection, ever | **no instrument** |

**These are not four instruments.** They are four channels of one telescope, and
the argument for building them together rather than separately is the whole point
of this paper:

1. **A relationship is only measurable if its terms are simultaneous.** Two
   excellent instruments with independent clocks and separate calibration
   histories are *worse* for this class of question than one modest instrument
   measuring both quantities together.
2. **A common calibration chain removes the handover systematic entirely.** A
   27.7σ step in detector output across a single handover is the scale of the
   obstacle facing any claim about slow structure; a single long-lived instrument
   has none.
3. **The marginal cost of an additional channel on an existing solar feed is
   small** relative to a standalone mission, which is why the combination has
   never been costed as a unit.

**Precedent exists for every component.** HARPS-N and NEID both operate
Sun-as-a-star solar feeds (Dumusque et al. 2015; Lin et al. 2022); SOHO/VIRGO
performed multi-channel irradiance photometry (Fröhlich et al. 1995); BiSON and
GONG have run helioseismic networks for decades (Chaplin et al. 1996; Harvey et
al. 1996). What does
not exist is the **combination**, at these precisions, on a common clock and
calibration chain. That makes this a describable instrument rather than a
speculative one.

### Why each channel is independently wanted

**Spectropolarimetry at 10⁻⁷.** The Sun is the only star whose surface magnetism
we resolve, and the only one for which a disc-integrated polarimetric signal can
be checked against resolved truth. Closing that gap calibrates every inference
drawn from unresolved stellar polarimetry, including exoplanet host
characterisation.

**Line-profile ratios at 10⁻⁷.** Chromospheric and transition-region diagnostics
rest on line ratios whose decade-scale stability has never been established at
this level. A stable ratio measurement constrains chromospheric heating models
directly.

**p-modes at 10⁻⁶.** GOLF (Gabriel et al. 1995) reached 1.46×10⁻⁵ over 25.9
years, as measured in Griffin (2026). A successor a factor
of fifteen better resolves activity-cycle frequency shifts at a precision that
discriminates between competing models of the near-surface shear layer.

**Core g-modes.** The only direct probe of the solar core, sought for fifty years
and still without an uncontested detection (Appourchaux et al. 2010; cf. García
et al. 2007), and the outstanding unsolved objective of helioseismology
independent of any application here.

---

## 4. Two spacecraft: what no telescope can supply

Two of the permanently blank columns are not observable from a solar telescope at
any precision:

**The interplanetary electric field** has never been measured continuously. It is
fundamental to solar wind acceleration and requires in-situ instrumentation.

**High-latitude solar wind** existed only while Ulysses flew, 1990–2009 (Wenzel
et al. 1992), and has had no successor. It requires an out-of-ecliptic orbit;
Solar Orbiter reaches only ~33° inclination (Müller et al. 2020). The gap is stark: one
mission, one snapshot, no plan to replace it — and it bears directly on the polar
field reversal that sets the activity cycle.

These are decadal-scale requests and are stated here for completeness of the
argument rather than as near-term asks. The point worth carrying into
prioritisation is that **the cost of their absence is now calculable**: each
removes a column from every multi-observable analysis for as long as it persists.

---

## 5. What is *not* being asked for

This paper makes no request for computing time beyond §2, and that is deliberate.

The entire search programme — 29 searches, 117 million statistic evaluations —
cost **1.1×10¹⁵ FLOP**: eighty days on a Cray-1 of 1976, 1.7 hours on a single
88-core workstation, and **under a millisecond of Frontier**. An exhaustive
quadruple sweep over 30 channels is 17.8 days on that same workstation and 0.2
seconds on a flagship machine.

Compute has not been the constraint for some time. Saying so plainly is part of
the argument: when the cheap resource is already saturated, what remains is the
expensive one — and here that is **joint measurement**.

---

## 6. Ordering and cost class

| | Request | Scale | Timescale |
|---|---|---|---|
| 1 | Dynamo ensemble for the higher-order null | 0.4–1.8M node-hours | one allocation cycle |
| 2 | Co-registered Sun-as-a-star multi-channel monitor | one instrument, four channels | mid-term |
| 3 | Interplanetary E-field; out-of-ecliptic solar wind | two missions | decadal |

The ordering is not arbitrary. **(1) is a precondition for interpreting (2)**,
and is roughly three orders of magnitude cheaper. **(2) closes four gaps with one
instrument** and removes the handover systematic that limits the existing record.
**(3)** is stated for completeness and is where the real money is.

---

## 7. Relevance to the decadal process

The request is unusual in shape: not a mission concept so much as a case for
**breadth and simultaneity across capabilities that already exist in principle**.
Three implications for prioritisation:

1. **Instrument handovers are a scientific risk, not a programmatic
   inconvenience.** A slow signal spanning a handover can be destroyed by
   cross-calibration or manufactured by it, and the GOES-15/16 step shows the
   scale. Overlap periods deserve treatment as science requirements with stated
   durations — the 1,031 dual days of GOES overlap are what made that control
   possible at all.

2. **Continuity has a value that can now be quoted.** The high-latitude gap is
   not merely regrettable; its cost to multi-observable analysis is calculable,
   and the same method prices any other proposed gap.

3. **Sun-as-a-star measurements are undervalued relative to resolved ones.** The
   Sun is the calibration source for all unresolved stellar observation, and the
   disc-integrated quantities that would serve that role are precisely the ones
   left uninstrumented.

---

## 8. Statement of limitations

The search programme that motivated this paper returned **twenty-five nulls,
three void results, and one re-detection of a known signal used as a positive
control**. It has detected nothing. The framework that motivated it is not
supported by its results and is not invoked to explain them.

The claim advanced here is narrower and does not depend on that framework:

> Coverage across *jointly measured* observables is the binding constraint on a
> class of solar time-series questions; this is now measured rather than
> asserted; and the measurements that would relieve it are ones heliophysics has
> independent reason to want.

The four blank columns are blank regardless of whether anything is encoded in the
filled ones, and the 36 handovers complicate slow-structure claims regardless of
their origin.

---

## References

Prior art is credited below. Entries marked **[verify]** are cited from the
author's recollection of the canonical reference and must be checked against the
literature before this document is submitted anywhere.

**The search programme this paper draws on**

Griffin, R. (2026). *A combination-space search for embedded technosignatures in
solar and heliospheric archives.* Radio Club of America. Code and data:
https://github.com/dxtrainc/behavioral-seti ·
https://dxtra.com/static/galactic-dx/

**Instruments**

Dumusque, X., et al. (2015). HARPS-N observes the Sun as a star. *ApJL* **814**,
L21. **[verify]**

Fröhlich, C., et al. (1995). VIRGO: experiment for helioseismology and solar
irradiance monitoring. *Solar Physics* **162**, 101. **[verify]**

Gabriel, A. H., et al. (1995). Global oscillations at low frequency from the SOHO
mission (GOLF). *Solar Physics* **162**, 61. **[verify]**

Harvey, J. W., et al. (1996). The Global Oscillation Network Group (GONG)
project. *Science* **272**, 1284. **[verify]**

Chaplin, W. J., et al. (1996). BiSON performance. *Solar Physics* **168**, 1.
**[verify]**

Lin, A. S. J., et al. (2022). The NEID solar feed. **[verify — journal and volume
not confirmed]**

Müller, D., et al. (2020). The Solar Orbiter mission. *A&A* **642**, A1.
**[verify]**

Wenzel, K.-P., et al. (1992). The Ulysses mission. *A&AS* **92**, 207.
**[verify]**

Kopp, G., & Lean, J. L. (2011). A new, lower value of total solar irradiance.
*GRL* **38**, L01706. **[verify]**

**Simulation codes**

Featherstone, N. A., & Hindman, B. W. (2016). The emergence of solar
supergranulation as a natural consequence of rotationally constrained
interior convection [Rayleigh]. *ApJ* **818**, 32. **[verify]**

Brun, A. S., Miesch, M. S., & Toomre, J. (2004). Global-scale turbulent
convection and magnetic dynamo action in the solar envelope [ASH]. *ApJ*
**614**, 1073. **[verify]**

Vögler, A., et al. (2005). Simulations of magneto-convection in the solar
photosphere [MURaM]. *A&A* **429**, 335. **[verify]**

**Solar g-modes**

Appourchaux, T., et al. (2010). The quest for the solar g modes. *A&A Review*
**18**, 197. **[verify]**

García, R. A., et al. (2007). Tracking solar gravity modes: the dynamics of the
solar core. *Science* **316**, 1591. **[verify]**

**Cited in the source paper, and carried here with its citation**

Eden, T. D., et al. (2024). Solar atmospheric oscillations as measured by GOES-R
EXIS EUVS-C. *ApJL* **973**, L18.

Wright, J. T., Kanodia, S., & Lubar, E. (2018). How much SETI has been done?
*AJ* **156**, 260.

Sheikh, S. Z. (2020). The nine axes of merit for technosignature searches.
*Int. J. Astrobiology* **19**, 237.
