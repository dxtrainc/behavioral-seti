# Note for integration — the gate is a cryptographic primitive, and naming it correctly changes the search order

**For:** Fable, for integration into the main document
**From:** session of 2026-09-12
**Target sections:** §2.2 (primary), §3.x tier ladder, §5 discussion, §6 conclusions
**Status:** analysis complete, numbers computed and checked. Not yet in either paper.

---

## 1. The claim

The draft calls the gate **proof of work** and reaches for the CAPTCHA analogy (§2.2, and the abstract).
That analogy is doing useful work but it names the wrong primitive, and the mis-naming has hidden a
methodological consequence that bears directly on what we run next.

Stated precisely, the gate is a **symmetric cipher under brute-force keysearch**:

| cryptographic object | our object |
|---|---|
| key | the tuple *(which channels, which functional form)* |
| ciphertext | the public archive |
| plaintext | the imposed modulation |
| keyspace | C(N,k) × \|F\| |
| distinguisher | the detection statistic |

Everything is public except the key — the archives, the physics, the method, this paper. That is
**Kerckhoffs's principle** in its purest available form, and it is not a coincidence: a sender who
wanted concealment would not modulate a star.

The property that makes the whole scheme work is **self-identifying plaintext**. The receiver knows
when the right key has been tried because the statistic fires; no crib, no known-plaintext pair, and
no side channel is needed. This is exactly the property that made DES brute-forceable — one can
recognise English without knowing in advance which English.

## 2. Why "proof of work" is the wrong name

In hashcash and in Bitcoin the **sender** performs the work, in order to prove commitment or to price
a message. Here the **receiver** performs the work. The direction of payment is inverted, and with it
the design goal.

What the gate actually is, is a **proof-of-capability challenge** — a CAPTCHA turned inside out. A
CAPTCHA is a puzzle a computer is supposed to *fail*. This is a puzzle only a sufficiently advanced
computer can *pass*. The draft's CAPTCHA sentence in §2.2 is right about the mechanism (the poser need
not know who solves it) and wrong about the polarity. Suggest keeping the sentence and adding the
inversion explicitly — it is a memorable line and it sharpens rather than weakens the section.

## 3. The nearest true antecedent, and where it breaks

§2.2 currently cites Clarke's *The Sentinel* as the nearest antecedent and notes, correctly, that it is
an antecedent for the intuition rather than a derivation. There is a real technical antecedent as well,
and it should be added:

**Rivest, Shamir & Wagner (1996), "Time-lock puzzles and timed-release crypto."** RSW pose a puzzle
whose solution requires a chosen amount of computation, and — this is the part that matters — they
explicitly reason about *projected future hardware speed* in order to set the difficulty so the puzzle
opens at a chosen future date. That is precisely the tier-ladder-as-clock argument, arrived at
independently and fifty years earlier in a different field.

**But the construction differs in one load-bearing respect, and the difference is favourable to us.**

A time-lock puzzle, and its modern descendant the verifiable delay function, is **inherently
sequential**: repeated squaring in a group of unknown order cannot be parallelised, so you cannot buy
your way out of the delay with more machines. Our keysearch is embarrassingly parallel. Measured
consequence:

| resource | time for T6 (2.8×10²³ FLOP) |
|---|---|
| El Capitan, dedicated | 6.2 d |
| Frontier, dedicated | 9.8 d |
| **Folding@home at its 2020 peak (2.43 EFLOPS)** | **1.3 d** |
| BOINC aggregate, typical | 108 d |

A volunteer network beats the flagship machine, because the problem admits parallelism the VDF
construction forbids.

**A designer who wanted a calendar clock therefore got the construction wrong. A designer who wanted
to gate on total civilisation compute got it exactly right** — and that is the more sensible thing to
gate on. Calendar dates are parochial; you cannot fake Kardashev level by waiting. Recommend this as a
short subsection or note box after the two-sided difficulty bound in §2.2, since it strengthens the
lower-bound argument already there: the gate measures *aggregate capability*, which is the quantity the
lower bound was trying to threshold in the first place.

## 4. Keyspace entropy, and why the data requirement is mild

Key = (combination, form). Required record length is set by the need for the true key's statistic to
stand above the maximum of K draws from the null; that maximum grows as √(2 ln K), and detectable
amplitude falls as 1/√N, so **N_required ∝ ln K**.

| tier | combinations | forms | H(K) bits | record length, relative |
|---|---|---|---|---|
| T1 pairs / 30 | 4.35×10² | 4 | 10.8 | 1.00× |
| T2 triples / 30 | 4.06×10³ | 4 | 14.0 | 1.30× |
| T3 quadruples / 30 | 2.74×10⁴ | 4 | 16.7 | 1.56× |
| T4 triples / 60 | 3.42×10⁴ | 10 | 18.4 | 1.71× |
| T6 quadruples / 60 | 4.88×10⁵ | 10 | 22.2 | 2.06× |
| combinatorial max, k=30 of 60 | 1.18×10¹⁷ | 10 | 60.0 | 5.58× |

**The keyspace spans fifteen orders of magnitude; the required record length spans 5.6×.** Compute
scales with K, data scales with log K.

Two things follow, both worth stating in the paper:

- **The gate is compute-bound, not data-bound.** That is precisely how one would design a time-lock
  intended to be opened eventually: make the barrier the thing that falls predictably (compute) rather
  than the thing that does not (how long our instruments have been running).
- It supplies an independent argument for the §4 finding that sensitivity is not the binding
  constraint. We reached that empirically (245× more sensitive at two minutes than at one day, yet the
  limiting quantity is combination-space coverage). The entropy argument reaches the same place from
  first principles.

*(Housekeeping: C(60,30) = 118,264,581,564,861,424 = 1.18×10¹⁷. An earlier verbal discussion in
session quoted 1.2×10¹⁸, which was wrong by a factor of ten. It never reached either manuscript;
flagging only so the wrong value is not reintroduced.)*

## 5. Unicity distance ≡ our multiple-testing threshold

Shannon's **unicity distance** is the ciphertext length below which more than one key yields a
plausible decryption — i.e. the length below which spurious keys decrypt. That is the same quantity as
the multiple-testing threshold: below it, wrong combinations fire.

**Our Bonferroni correction is a unicity-distance calculation in disguise.** This is worth a sentence
or two because it reframes the N²/α scaling result (surrogate count must scale with test count) as an
instance of a known information-theoretic bound rather than as a bespoke statistical fix. It also
supplies the right vocabulary for the §4 void results: a void is not a failed search, it is a search
conducted below the unicity distance for its keyspace.

## 6. The consequence that changes what we run — weak keys

**This is the actionable item, and I would put it in §5 and in §6.**

A cipher is strong when keys are drawn uniformly from the keyspace. Ours are not. The sender is
constrained: the chosen combination must be one the receiver *plausibly measures*, on records long
enough to hold a message. The key distribution is therefore concentrated on the best-instrumented,
most obvious observables. In cryptographic terms these are **weak keys**, and against weak keys brute
force is the wrong attack.

**No password cracker starts with brute force. It starts with a wordlist.**

The wordlist writes itself from the designer's own constraints — the channels any civilisation would
measure early, precisely, and for a long time:

| channel | record | precision | why a designer picks it |
|---|---|---|---|
| TSI | 49 yr | 29 ppm | the single most obvious stellar observable |
| sunspot number | 270 yr | integer | longest continuous record humans have |
| F10.7 | 79 yr | 0.1 sfu | the canonical activity proxy |
| p-mode frequencies | 30 yr | µHz | the only channel probing the interior |
| MgII index | 45 yr | 10⁻⁴ | chromospheric, instrument-independent |
| SSI band ratios | 20 yr | 100 ppm | dimensionless by construction |

| attack | keys | cost on zeus |
|---|---|---|
| brute-force T3, all 30 channels | 27,405 quadruples | ~18 d |
| **top-6 wordlist: 15 pairs + 20 triples + 15 quadruples** | **50** | **seconds** |

A **548× reduction**, and it is not a heuristic shortcut — it is the correct attack against a key
distribution we have positive reason to believe is non-uniform.

**Therefore: running brute-force T3 before the wordlist is exhausted is a methodological error, not
merely an inefficient ordering.** This lands on the same recommendation §3.6 already makes on
designer-side grounds (the most plausible first-contact carriers sit in archives that exist, are
resolved at the required level, and have not been searched) — but it arrives from a different
direction, which is worth saying explicitly, because two independent arguments converging on the same
search order is itself evidence the ordering is right.

Rows 3–7 are not low-hanging fruit. **They are the dictionary attack.**

## 7. Suggested placement

| where | what |
|---|---|
| §2.2 | rename the primitive: keysearch under Kerckhoffs, not proof of work; keep the CAPTCHA sentence, add the polarity inversion (§2 above) |
| §2.2, after the two-sided bound | RSW time-lock antecedent + the sequential/parallel distinction; gate measures aggregate compute, not calendar time (§3) |
| §3, tier ladder | keyspace entropy table; compute ∝ K, data ∝ log K (§4) |
| §4 or §5 | unicity distance ≡ multiple-testing threshold; recasts void results (§5) |
| §5 and §6 | weak keys → dictionary attack → rows 3–7 before T3 (§6). **The one that changes what we do.** |

## 8. Caveats, stated plainly

- The weak-key argument depends on a prior over what a sender would choose. That prior is
  *our* reasoning about *their* design, and it is not testable. It justifies search **order**, which
  costs nothing if wrong. It would not justify restricting the search **space**, and the note should
  not be read as licence to drop the brute-force tiers.
- The cipher mapping is structural, not a claim that any cryptographic security proof transfers. No
  hardness result is being imported; the keyspace arithmetic and the √(2 ln K) threshold scaling stand
  on their own.
- Nothing in this note is a result. It is a reframing plus one scheduling consequence. It should not
  add to the searches count (25) or to any sensitivity claim.

---

# PART II — the dictionary attack was run. Result: void, with a precise diagnosis.

**Added 2026-09-12, same session.** Part I above was written before the run. Sections
9–12 report what happened when it was executed, and they change the recommendation in
§6. Fable: §9 and §11 are the ones that belong in the paper.

## 9. What was run

Two of the six wordlist channels do not exist as daily series — p-mode frequencies
(BiSON returns 403, GONG serves full-disk FITS behind an interactive form, VIRGO is
the IPv6-only FTP) and SSI bands (GOES EUVS is 1,945 days, below the 2,000-day
overlap floor). They were replaced by the two longest-record channels meeting the same
design criteria. Final wordlist, all ≥27 yr:

**TSI, F10.7, cosmic ray (Oulu NM), sunspot, X-ray background, MgII.**

All 15 quadruples are reachable on a 9,786-day (26.8 yr) six-way common overlap —
against 28% of brute-force triples that were not. Three four-body forms give 45 tests.

Of the 50 dictionary keys, the **pairs and triples were already committed and null**
(min p 0.0094 vs a 4.7×10⁻⁵ threshold; 0.0063 vs 4.5×10⁻⁶). They were not re-run. The
quadruples had never been tested, because no four-body form existed and because the
null that voided the triple sweep would have voided them too.

## 10. The matched null of §4.9 — built, and calibrated

§4.9 specified a surrogate that "preserves every pairwise cross-correlation and
destroys only the three-way alignment — the third-order analogue of IAAFT," and called
it "the concrete next step for this programme." **It is a one-line construction.**

Give every channel the *same* random phase screen φ(f). The cross-spectrum
S_ij = X_i X_j* picks up e^{iφ}e^{−iφ} = 1, so every pairwise cross-correlation at
every lag is preserved **exactly**; each power spectrum is untouched; but the
bispectrum picks up e^{i[φ(f₁)+φ(f₂)−φ(f₁+f₂)]}, which is random unless φ is linear in
f, so three- and four-way alignment is destroyed. Channels are normal-score
transformed first so the marginal is preserved too.

Measured preservation: **6.4×10⁻¹⁷**. Exact by construction, not annealed toward — a
Schreiber-style constrained-annealing attempt reached only 8%.

Calibration took three iterations, and the failures are worth recording because two of
them would have shipped:

| construction | over-firing | tests failing uniformity |
|---|---|---|
| index-axis surrogate, gap-spanning increments kept | **9.39×** | 34 / 45 |
| index-axis surrogate, gap-spanning increments dropped | 0.84× | 15 / 45 |
| **calendar-symmetric surrogate, increments dropped** | **1.02×** | **0 / 45** |

Two lessons, both generalisable:

- **The first gap control was a no-op.** Its synthetic data was phase-randomised on the
  compressed index axis, hence smooth in index space, hence had no gap-jumps to mask.
  Masking-on and masking-off gave identical results and the pipeline looked clean
  either way. Only rebuilding the control so the synthetic data carried real calendar
  structure exposed the 9.39× — which is the same order as the ninefold tail excess
  that voided the triples. *A control that cannot fail certifies nothing.*
- **Surrogates must share the data's calendar, not just its spectrum.** The data is a
  smooth process sampled on an irregular calendar with 1,064 breaks; randomising on the
  compressed index axis, where the calendar has been flattened away, leaves a residue
  that masking alone cannot remove (15/45 still failed). Randomising on the full grid
  and subsampling through the same mask closes it completely.

## 11. THE RESULT — void, and why. This is the part for the paper.

Against the calibrated null (1.02×, 0/45 uniformity failures), **real solar data sits
at median z = 45.3, maximum 426.9, with 30 of 45 tests beyond |z| > 20.**

| form | min z | median z | max z |
|---|---|---|---|
| 3rd-diff | 45.3 | 52.1 | 426.9 |
| cross | 35.5 | 48.2 | 222.0 |
| quad-prod | 1.6 | 5.1 | 12.0 |

**This is not a detection. It is the null being wrong for real data, and the diagnosis
is clean** because the same pipeline returns uniform p (1.02×) on gap-realistic
synthetic data and z = 45 on the Sun. The only difference between those two datasets is
intrinsic higher-order structure: the solar cycle is asymmetric, activity is
multiplicative, flares are intermittent. Common-phase randomisation destroys all of it,
so real solar data is extreme *by construction*.

> **The quadruple space is void — but for a sharper reason than the triples were.** The
> triples failed because the null was unmatched. This null *is* matched, and calibrates
> perfectly. It fails because **an omnibus test for four-way structure cannot work on a
> channel that was never random.** The Sun supplies the structure the test looks for.

**And the cryptographic framing of Part I predicts exactly this.** No cryptanalyst asks
whether a ciphertext is non-random — such tests fail on any real channel. The question
is always *does this key produce recognisable plaintext*, calibrated against **decoy
keys** run through the identical pipeline. The programme already does this everywhere
it works: the decoy-offset control in the annual-geometric-term search (adopted after a
Gamma(5) null over-fired 3–5×), the Earth/Mars cross-viewpoint test, the CH_D dark
channel in the EVE ESP acquisition. **The fix for the quadruple space is not a better
surrogate. It is replacing an omnibus question with a targeted one.**

### 11b. A structural theorem, found while designing the injection

Suppose a sender adds c_k·s(t) to channel k. A linear four-body form picks up
(c_a − 3c_b + 3c_c − c_d)s; a linear pair form picks up (c_i − c_j)s. For the signal to
be invisible to **every** pair, all c_k must be equal — and then the four-body form
picks up (1−3+3−1)s = 0.

**A linear four-body form cannot carry a signal that pairs cannot also see.** 3rd-diff
and cross are higher-order differences, not independent carriers. Only the
multiplicative form can: a four-way phase alignment lives in the trispectrum and leaves
every pairwise cross-correlation untouched.

This is not an aside — it is corroborated by the table above. The two linear forms,
which provably carry nothing new, are the two swamped at z ≈ 50. **quad-prod, the only
genuinely four-body form, is the only one anywhere near its null (median z 5.1).** The
quadruple space is genuinely new only through the multiplicative form, and any future
enumeration should drop the linear ones rather than spend trials on them.

## 12. Errors made in this session, for the record

- **The injection could never have fired.** With 400 surrogates the p-floor is 1/401 =
  2.5×10⁻³, *above* the Bonferroni α of 1.11×10⁻³. Every recovery fraction was
  structurally zero. **Surrogate count must exceed 1/α**, and by a margin. Belongs with
  the N²/α scaling discussion.
- **The first gap control was a no-op** (§10). Caught only by running the arm that was
  supposed to break.
- **C(60,30)** is 1.18×10¹⁷, not the 1.2×10¹⁸ quoted verbally earlier. Never reached
  either manuscript.

## 13. Revised recommendation, superseding §6

§6 said rows 3–7 are the dictionary attack rather than low-hanging fruit. **That still
holds for the pair space.** What Part II adds:

1. **Do not run brute-force T3.** Not because it is expensive — because it would push
   27,405 keys through an omnibus null that the Sun's own nonlinearity rejects at
   z ≈ 45. It would return 27,405 keys' worth of void. The dictionary's real value was
   never the 548× saving; it is that 15 keys is small enough to afford a null one can
   prove correct, and proving it correct is what revealed the space is not yet testable.
2. **The matched null of §4.9 is delivered and calibrated** and should be reported as
   such — it is a genuine methodological result independent of the void, and it applies
   to the triple sweep as well as the quadruple one.
3. **Retire the linear four-body forms** on the theorem in §11b.
4. **The next real step is a targeted test with decoy keys**, not a larger enumeration:
   the dark channel, the twin instrument (GOES-16/17), and the cross-viewpoint geometry
   are the controls that can distinguish an imposed alignment from solar nonlinearity.

---

# PART III — the obvious rescue was tried and rejected

**Added 2026-09-12, same session, after Part II.** §11 concluded the quadruple space
is void because the Sun's intrinsic higher-order structure swamps an omnibus test.
There is one obvious objection to that conclusion, it was tested, and it failed in a
way that is more instructive than the original void. §14 is the part for the paper;
§15 is the part that should change how the programme picks its next test.

## 14. The stationarity objection, tested

**The objection.** The common-phase null of §10 assumes the series are stationary.
The Sun plainly is not: activity is amplitude modulated by the eleven-year cycle, so
variance rises and falls. A stationary surrogate therefore looks unlike real data for
a reason having nothing to do with four-way alignment, and the median z = 45.3 of §11
might be measuring nothing more interesting than *the solar cycle exists*.

**The test.** Estimate a slowly varying robust scale s_k(t) per channel; standardise
u_k = d_k/s_k; apply the common phase screen to u — preserving every pairwise
cross-spectrum of the standardised series, exactly as before — then multiply the
envelope back, so the surrogate carries the data's own amplitude modulation.

The window is **scanned, not chosen**. It is a sensitivity limit rather than a tuning
knob: any modulation slower than the window is absorbed into the envelope and becomes
invisible, so a 365-day window means no sensitivity to modulation slower than a year.
Both arms — the z-reduction and the calibration — were specified before the run, with
the decision rule stated in advance: *a window counts only if it gives both a small z
and a clean calibration; a small z with a broken calibration means the surrogate is
absorbing signal rather than explaining structure.*

| window | median z | max z | \|z\|>20 | over-firing | failing uniformity |
|---|---|---|---|---|---|
| none (common-phase) | 45.3 | 426.9 | 30/45 | **1.02×** | **0/45** |
| 91 d | **0.8** | 97.7 | 2/45 | 4.25× | 32/45 |
| 365 d | 16.6 | 184.9 | 18/45 | 4.41× | 22/45 |
| 1461 d | 1.5 | 300.5 | 12/45 | 14.25× | 45/45 |

**Rejected.** The z collapse from 45.3 to 0.8 is not a result. Driving z down by a
factor of fifty *while breaking calibration by 4.25×* is the signature of a null given
enough freedom to absorb the data rather than account for it — and no window is
calibrated at all, the longest failing every one of the 45 tests.

> **Had the second arm not been run, the first arm alone would have been reported as
> the quadruple space reopening.** This is the same shape as the error that produced
> 63,856 spurious candidates under a Gamma(5) null, and the same shape as the gap
> control in §10 that was a no-op. Three times in one programme, a control that looked
> like a success was one unrun arm away from a false result.

**What this does NOT show.** It does not show that the Sun's four-way structure is
irreducible. The calibration failure is in the surrogate, not in the Sun, and a
construction bug in the envelope estimator cannot presently be distinguished from a
wrong hypothesis about solar variability. Further tuning was **declined deliberately**:
this was the third null construction on the same question, and continuing to adjust a
null until it passes is how one ends up with a null quietly fitted to produce the
wanted answer. The honest statement is that the stationarity objection is unresolved
and that resolving it is not on the critical path — see §15.

## 15. Three nulls, three false assumptions — the portable result

This is the most transferable thing to come out of the session, and it deserves a
place in the discussion rather than a footnote.

| null | load-bearing assumption about the Sun | how it failed |
|---|---|---|
| circular shift (§4.9) | pairwise structure is irrelevant to a three-body statistic | destroys it; ninefold tail excess; voided the triples |
| common-phase (§10) | the series are stationary | they are not; real data at median z = 45.3 |
| envelope (§14) | the nonlinearity is pure amplitude modulation | absorbs signal; 4.25× over-firing, 32/45 failing |

Each null required a modelling assumption about solar behaviour in order to be
constructed at all. Each assumption turned out to be false. And there is no reason to
expect the fourth attempt to be different, because the underlying difficulty is not
statistical: **we are trying to specify what the Sun would do in the absence of a
signal, for a star whose nonlinear behaviour is not understood well enough to
specify.**

**A differential test escapes the problem entirely.** Two instruments viewing the same
source share whatever the Sun does: solar nonlinearity is common-mode and cancels in
the difference, so it never has to be modelled. Nothing needs to be assumed stationary,
amplitude-modulated, or anything else.

This is not a new idea in the programme — it is a description of everything in it that
has already worked:

| test | what plays the role of the control | outcome |
|---|---|---|
| p-mode comb | GOES-16 **and** GOES-17 must both show it | detection, held up |
| annual geometric term | decoy offset, after Gamma(5) over-fired 3–5× | null, calibrated |
| cross-viewpoint | Earth line vs Mars line | null, with a stated limit |
| EVE ESP (in progress) | CH_D dark channel through the identical pipeline | pending |

**Every differential test in this programme has produced a usable result. Every
omnibus test has gone null-by-luck or void.** That is now four-for-four and
three-for-three, and it is a strong enough pattern to promote from observation to
design rule.

> **Proposed design rule, for §5 or §6: a search in this programme should be
> constructed so that the control is a second measurement rather than a model of the
> source.** Where no second measurement exists — a second instrument, a second
> viewpoint, a dark channel, a decoy key — the search should be deferred rather than
> run against a modelled null, because the modelled null will encode an assumption
> about solar behaviour that we are not in a position to defend.

This also explains, retrospectively, why the pair sweep was trustworthy and the triple
sweep was not. The circular shift is a *legitimate* null for a two-body statistic —
it destroys exactly the alignment being tested and nothing else. It stops being
legitimate at three bodies because there is now pairwise structure for it to damage.
The programme's reach in combination space is therefore bounded not by compute, as §4
argued, and not by coverage, but by **the order at which a defensible null stops
existing** — which is two for shift-based nulls, and appears to be two for everything
else tried so far.

## 16. Revised recommendation, superseding §13 item 4

§13 recommended a targeted test with decoy keys as the next real step. §14 and §15
sharpen the ordering:

1. **Rows 3–7 first.** They are pair searches. The pair null is the one construction
   in this programme that is defensible, the searches need no new methodology, and
   §3.6 has been pointing at them for four revisions. Nothing in Parts II or III
   touches them.
2. **The twin-instrument test (GOES-16 vs GOES-17) is the method for reopening the
   higher-order space**, and should be framed that way rather than as one more search:
   it is the first higher-order test in the programme that requires no model of solar
   variability. Note it is a fast-band test — the two-spacecraft overlap is 1,945 days
   at one-minute cadence — so it is not a substitute for the daily quadruple sweep but
   a different and better-controlled question.
3. **Do not run brute-force T3**, for the reason in §13 item 1 and now more strongly:
   there is no defensible four-body null to run it against.
4. **Retire the linear four-body forms** (§11b theorem, corroborated by the z-table).
5. **The stationarity question stays open and is not on the critical path.** If it is
   ever worth settling, settle it with a second instrument rather than a better
   surrogate.

---

# PART IV — a calibration standard the paper does not yet set

**Added 2026-09-12, after rows 5, 6 and 7 were run.** §17 is the part for §5.6. It comes
out of a specific failure today that the paper's current standard would have let through.

## 17. Injection validates recovery. It does not validate calibration.

The paper's standard for a limit is §4.2's: *the amplitude recovered 95% of the time
under injection*. That is the right standard for one failure mode and silent about
another, and the two are worth separating explicitly.

**Injection asks:** if a signal of size _A_ is put into the pipeline **at this point**,
does it come out? That catches a dead detector, a mis-set threshold, a statistic with
no power. It is necessary.

**It cannot catch a pipeline that attenuates the signal upstream of the injection
point.** If the injection goes in after the step that destroys the signal, recovery is
perfect and means nothing.

> **Row 5 is precisely that case, and it is the reason to write this down.** The dip
> search detrended against a running median of one day — only four times the duration
> of the transit it was meant to find — so the detrend was absorbing the very feature
> the statistic looked for. **An injection into the detrended residual would have
> recovered beautifully and certified a broken search.** What caught it was Venus:
> a real occultation whose depth is (R_venus/R_☉)² = 75.6 ppm, fixed by two radii
> before anyone touches the data. The search could not see it, and said so.

### Not all anchors are equal, and the distinction is worth a sentence in the paper

| kind of anchor | what it fixes | what it proves |
|---|---|---|
| **position** — a known signal at a known frequency | where | the detector can *fire* |
| **amplitude** — a known signal of known size | where **and how big** | the detector fires at the *right size*, so a limit is in real units |

Of the four rows run today:

| row | anchor | kind |
|---|---|---|
| 5 — occulter transits | Venus 75.6 ppm, Mercury 12.3 ppm, from disk-area ratios | **position + amplitude** |
| 7 — p-mode frequencies | solar-cycle shift ≈ 0.4 µHz p-p, from the literature | **position + amplitude** |
| 6 — VIRGO SPM coherent | p-mode envelope at 3.09 mHz | position only |
| 3 — EVE ESP sub-minute | p-modes in EUV irradiance | position only |

Rows 3 and 6 can prove their pipelines fire. **Their limits rest entirely on injection
and inherit its blind spot**; nothing in either search establishes that a quoted ppm is
a real ppm. Rows 5 and 7 can establish that, and row 5's anchor is the strongest control
anywhere in this programme because it contains **no free parameter at all** — the depth
is two radii, not a fit, not a calibration, not a quantity derived from the data.

### The proposed standard

> **Where a channel contains a signal of known amplitude, the search must be calibrated
> against it before any limit is quoted, and the recovered amplitude reported alongside
> the predicted one. Where no such signal exists, the limit should say so** — because a
> limit resting on injection alone is a statement about the detector, not about the sky,
> whenever any upstream step can attenuate.

This is sharper than §5.6's current position, which credits physics-supplied controls
with removing two apparent detections but treats them as a way of **killing false
positives**. Today's row 5 shows the other half: an amplitude anchor also kills a
**false limit**, which is the more dangerous error because nothing downstream ever
contradicts it. A wrong null is quietly wrong forever.

It also costs nothing to adopt retrospectively. §4.6's p-mode re-detection already
functions as a position anchor for the fast band; what it does not do is fix an
amplitude, and saying so in one sentence would tell a reader exactly how far to trust
the fast-band limits.

## 18. Three estimates that came in optimistic against measurement

Recorded together because they are one failure shape, and it is not the window-size
class of §15 — it is **a number reasoned about rather than measured, wrong in the
flattering direction**:

| claim | reasoned | measured |
|---|---|---|
| §3.6's reach for row 7 | 10⁻⁵–10⁻⁶, margin 0.1–1 | **1.94×10⁻⁵, margin 0.05** |
| peak-bagging would gain ~6× | from published mode-fitting precision | **1.33×** |
| row 5 gate v1 | "6 of 6 transits recovered" | all six were noise; 0 of 6 survive a control-day null |

The first is the paper's own, and §3.6's footnote already admits those columns are
order-of-magnitude except for rows 1, 8 and 9. **Row 7 is now the first of the estimated
rows to receive a measurement, and it came in worse than the estimate.** That is one
data point, not a pattern — but it is the only data point there is, and the remaining
estimates in that column should be read with it in mind rather than at face value.

---

# PART V — scored against Sheikh's nine axes

**Added 2026-09-12.** §19 is for §1.4 or §5. Sheikh (2020), *Int. J. Astrobiology* **19**, 237,
"Nine axes of merit for technosignature searches" — a framework **developed at the 2018 NASA
Technosignatures Workshop** and now the field's standard instrument for comparing search
strategies. The paper cites it in the references and does not score against it. It should,
because the exercise is not flattering everywhere and the places it is not are the places a
reviewer from that community will press.

## 19. The nine axes

Four axes are functions of *us*; five are functions of the technology sought.

### The practical four

| axis | score | why |
|---|---|---|
| **Cost** | **very high** | Zero new observations. Twenty-nine searches, every one against a public archive collected for unrelated purposes, run on one workstation. No telescope time was requested and none is needed. This is the single strongest axis and the paper under-sells it |
| **Observing capability** | **high, with a correction** | The data exists now. But §3.6's capability estimates have been measured three times and came in short each time — row 5 proved **unreachable** (6 h noise floor 31–108 ppm against a 76 ppm transit), rows 6 and 7 at margins 0.5 and 0.069 against an estimated ~1. Capability is real but narrower than the table claimed |
| **Ancillary benefits** | **moderate** | The sweeps recover genuine heliophysics unprompted — the four surviving solar-wind pairs are textbook, the p-mode comb is a known signal re-detected. Two bands are now *measured* as unsearchable at ppm level (VIRGO 3 d–2 months; transit timescales at 6 h), which is a real space-weather statement. But no new astrophysics is produced, and the honest ancillary product is **methodological**: the matched four-body null, and the control-design ledger |
| **Detectability / sensitivity** | **mixed** | Deep in absolute terms — **1.4×10⁻⁶** in fractional Lyman-α at two minutes, **2.0×10⁻⁷** in VIRGO broadband. Marginal against the level the framework's own logic predicts: the best margin achieved anywhere is **0.5** |

### The technology-dependent five

| axis | score | why |
|---|---|---|
| **Duration (Drake's _L_)** | **very high** | A passive modulator redistributing flux needs no power, no consumables and no maintenance (§2.6). A radio beacon must be *run*; this must only persist. **_L_ could be geological**, which is the axis on which this architecture most clearly beats a transmitter, and the paper barely argues it |
| **Ambiguity** | **high in principle, poor in practice** | In principle low: a dimensionless relationship among independently measured observables is not something natural processes have reason to impose. In practice this programme found **four artefact families** in one day — instrument lines with ±2% sidebands, solar-cycle harmonics, harmonics of its own detrend window, and the Nyquist frequency — and three were unanticipated. Ambiguity is the axis where the gap between theory and measurement is widest |
| **Extrapolation** | **moderate** | §2.6 sizes the modulator at ~1.5×10¹² m² of film for 1 ppm. Far beyond us, but it is passive structure rather than exotic physics — easier than a Dyson sphere, harder than a radio beacon |
| **Inevitability** | **LOW — the weakest axis, and it should be conceded** | The argument requires a civilisation to want to signal, to prefer a capability gate to a broadcast, **and** to choose a combination we happen to measure. §2 argues each is the *efficient* choice, but efficient is not inevitable. Radio leakage and waste heat are far more inevitable because they require no intent at all. **Paper 2's own arithmetic puts the chance we measure the right pair at ~21%, and the right triple at ~9%** — that is an inevitability statement, and it is not a good one |
| **Information content** | **very high** | Unlike an infrared excess or a transit anomaly, which report only that something is there, this architecture carries a message by construction — the staged pointer → node → dictionary → content of §2.4. The decode *is* the demonstration of capability |

## 20. What the profile says, including the part that is unflattering

**Strong where transmitters are weak.** Cost, duration and information content are the three axes on which a passive stellar modulator beats a radio beacon outright, and the paper argues only the third of them properly. §2.6 has the persistence argument; it is not connected to Drake's _L_, and it should be.

**Weak exactly where a workshop reviewer will press.** Inevitability is the soft axis, and the honest response is to concede it rather than argue: this is a search for a *deliberate*, *designed* signature, and deliberateness is less inevitable than leakage. What the paper can say in exchange is that the search costs nothing — the archives are already collected — so a low inevitability prior is affordable in a way that a telescope-time proposal would not be.

**And ambiguity is where today's work should change the claim.** §1.4 argues the dimensionless carrier is hard for nature to imitate. That is true and it is not the binding problem: the binding problem is that the *analysis* generates ambiguity faster than the sky does. Four artefact families in one day, three unanticipated, and seven controls that could not fail. **The paper's ambiguity score is limited by the pipeline, not by the physics**, and §5.6 is the right place to say so.

> **Placement.** §1.4 already separates what is borrowed from what is new, and cites Wright,
> Kanodia & Lubar (2018) for the coverage-as-merit standard rather than claiming it. A nine-axis
> table belongs beside that, because it answers the question a NASA-workshop reader will actually
> ask — *where does this sit among the alternatives* — and answers it including the two axes where
> the answer is "worse".

## 21. Where this sits against the NASA technosignature taxonomy

**For §1.1 or §1.4, and a reader will want it early.** The 2018 NASA Technosignatures Workshop
(Houston; report at arXiv:1812.08681) is the field's reference organisation of what a
technosignature *is*. A reader arriving from that report needs one paragraph telling them which
box this paper occupies, and the answer — that it occupies none of them — is worth stating
plainly rather than leaving them to infer it.

**Every category in the report assumes one of four things.** That is the useful way to present
it, because it makes the omission visible:

| NASA category | what it assumes exists | our position |
|---|---|---|
| **Radio** — narrowband or pulsed communication, the Cocconi–Morrison line onward | a **dedicated transmitter**, powered and pointed | not assumed |
| **Optical / infrared** — laser pulses, continuous beacons | a **dedicated transmitter**, higher power still | not assumed |
| **Megastructures and waste heat** — Dyson spheres, IR excess, transit anomalies | a **structure large enough to occult or re-radiate** a detectable fraction of the star | not assumed; §2.6's modulator is ~10⁻⁶ of the disc, invisible as an occulter |
| **Atmospheric technosignatures** — CFCs, NO₂ and industrial species in exoplanet atmospheres | an **industrial byproduct**, i.e. pollution as a side effect | not assumed; nothing here is a byproduct |
| **Solar-system artifacts** — lunar, asteroidal and Lagrange-point surveys | a **physical object** to be imaged or encountered | not assumed as the carrier; §3.6 rows 4, 12 and 13 note the strand and it belongs to Davies & Wagner and Benford rather than here |
| **Planetary surface** — city lights, artificial albedo | a **modified planet** | not assumed |

> **The premise of this paper is that none of those is necessary.** The carrier is emission the
> star already produces for its own reasons, **redistributed rather than generated** (§2.6), so
> there is no apparatus to build, no structure to see, no byproduct to accumulate and no object
> to encounter. That premise does not appear in the workshop taxonomy, and the searches that
> follow from it — of solar and heliospheric archives collected for space weather — do not appear
> in its inventory of what has been looked at.

**Three further differences a reader will want named, and one debt:**

- **The carrier is a relationship, not a channel.** Every NASA category is organised by *where to
  look* — this band, this line, this species. The combination space says the signal may live in a
  dimensionless relationship among observables and in none of them individually (§3.1, §3.3).
- **The gate is opened by the receiver, not the sender.** Proof-of-capability rather than proof of
  work: the difficulty is a filter on who can act, not a price the sender pays (§2.2, §2.8). We
  have found no prior proposal of a capability-gated signal.
- **The target is our own star, from archives already collected.** The report's emphasis is on new
  instruments and on piggybacking planned missions. Nothing here was observed for SETI.
- **The debt is to Wright, Kanodia & Lubar (2018).** The standard this paper is written to —
  *that a null result is worth the fraction of a defined space it excludes* — is theirs, and §1.4
  says so. The contribution here is a space in which that fraction can be **enumerated**, and
  §4.3's finding that coverage rather than sensitivity is what binds.

**And the honest sentence to end on**, which belongs beside the nine-axis table of §19: this
architecture scores better than a transmitter on cost, on duration and on information content,
and **worse on inevitability than any signature that requires no intent at all**. A reader from
that workshop will grant the first three immediately and press on the fourth, and the paper is
stronger for conceding it in the same breath.
