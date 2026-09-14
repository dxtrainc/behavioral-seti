# Inverting anticryptography: a capability-gated interstellar signal as a cipher under brute-force keysearch

**Robert Griffin**
Dxtra Inc. · rtg@dxtra.com

---

## Abstract

The SETI message-design literature is explicitly **anticryptographic**: from
Freudenthal's *Lincos* through Drake's Arecibo tests to contemporary message
design, every construction assumes a sender who maximises decodability —
maximum redundancy, minimum assumed context. We consider the inverse. Suppose a
sender, constrained by transmission cost, wishes the message to be legible *only*
to a receiver above some capability threshold, and therefore makes decoding
deliberately expensive.

We show this construction is a **symmetric cipher under brute-force keysearch**
with an unusual property: the plaintext is *self-identifying*, so the receiver
needs no crib, and the ciphertext is a *public record* the sender did not
transmit. We then observe that for such a scheme **Shannon's unicity distance and
the multiple-testing threshold of the receiver's search are the same quantity** —
the sender's difficulty parameter and the receiver's statistical burden are one
number seen from two sides.

Three consequences follow. Required record length scales as ln *K*, so keyspace
expansion is logarithmically cheap for the sender and linearly expensive for the
receiver — an unusually favourable asymmetry for a gating construction. A
receiver searching below the unicity distance for its keyspace obtains **no
result**, not a negative one, a distinction we argue must be reported separately
or accumulated coverage estimates are corrupted. And we prove that linear
multi-channel forms cannot carry a signal invisible to pairwise analysis, which
constrains what such a sender could choose.

We report an instantiation over ~10⁴ keys on 30 public time series, in which the
predicted below-unicity failure occurred and was diagnosed.

**No hardness result is claimed and no security proof is transferred.** The
contribution is the inversion, the correspondence it makes visible, and one
theorem.

**Keywords:** anticryptography · SETI · unicity distance · brute-force keysearch
· self-identifying plaintext · time-lock puzzles · multiple hypothesis testing ·
Kerckhoffs's principle

---

## 1. Introduction

### 1.1 The tradition this inverts

Sixty years of interstellar message design share one assumption: **the sender
wants to be understood.** Freudenthal's *Lincos* (1960) builds a language from
logic upward so a receiver with no shared context can bootstrap it. Drake tested
an early Arecibo message on colleagues at Green Bank in 1961 specifically to
measure decodability. The tradition has an explicit name —

> Freitas sets out a **Principle of Anticryptography**: a beacon message should
> be optimised for easy decoding by its recipients, using simple, highly
> redundant coding with abundant clues — the opposite objective to that of
> secret code-making. [7]

— and contemporary work continues it: general-purpose binary languages for
crowdsourced messages [8], controlled decipherment experiments in which one
researcher encodes and another decodes [9], and public decrypt challenges [10].
Every such construction maximises redundancy and minimises assumed context.

**This paper considers the inverse.** Suppose a sender is constrained by
transmission cost, and wishes the message to be acted upon *only* by a receiver
above some capability threshold — because a reply from a receiver who cannot act
on it is wasted energy. Such a sender has reason to make decoding **deliberately
expensive**, and to let the expense do the filtering.

To our knowledge this inversion has not been proposed. The nearest idea in the
literature is the zoo hypothesis [11], which shares the *prediction* — silence
despite presence — but supplies no mechanism, and requires many parties to
abstain indefinitely. A difficulty-gated signal requires nobody to abstain: it is
simply unrecognised until someone spends enough.

### 1.2 Why this is a cryptographic object

Once difficulty is the design goal, the construction becomes a cipher, and the
standard analysis applies. Classical cryptanalysis assumes an adversary who must
*recognise* a correct decryption; Shannon fixed the condition under which that is
possible, the **unicity distance**

  *U* = *H*(*K*) / *D*

with *H*(*K*) the keyspace entropy and *D* the plaintext redundancy per unit
length. Below *U*, spurious keys decrypt plausibly and the attack cannot
conclude.

Our setting has three unusual features:

- the **ciphertext is a public record** the sender never transmitted — an
  observational archive a third party collected for unrelated reasons;
- the **key is a selection**, not a bit string: which quantities are combined and
  under which functional form;
- the **plaintext is self-identifying** — the receiver knows the key is right
  because a statistic fires. No crib, no known-plaintext pair, no channel back to
  the sender.

The receiver's task is a brute-force keysearch with a statistical oracle. The
observation of this paper is that the stopping condition for that search is the
unicity distance in different vocabulary, and that this identity is what ties the
sender's difficulty parameter to the receiver's statistical burden.

### 1.3 Contributions

1. **The inversion (§1.1).** Sixty years of interstellar message design assume a
   sender maximising decodability. We consider a sender minimising it, and show
   the resulting object is a cipher. We find no prior proposal of a
   difficulty-gated interstellar signal; the zoo hypothesis shares the prediction
   and supplies no mechanism.
2. **The correspondence (§4).** Unicity distance and the multiple-testing
   threshold are the same quantity — the sender's difficulty parameter and the
   receiver's statistical burden seen from two sides. A Bonferroni or Šidák
   correction over *K* candidate keys is a unicity-distance calculation.
3. **The scaling (§5).** Required record length grows as ln *K*. Keyspace
   expansion is logarithmically cheap in data and linearly expensive in compute —
   the opposite of the usual intuition that a larger search space demands
   proportionally more evidence.
4. **The void/null distinction (§4.3).** A search below the unicity distance for
   its keyspace returns *no result*, not a negative one. We argue this should be
   a reported category.
5. **Payment inversion (§3).** The construction is not proof-of-work: the
   receiver performs the work. We give the consequences for what the scheme can
   and cannot enforce.
6. **A structural theorem (§6).** For linear multi-channel forms, any signal
   invisible to every pairwise test is invisible to the multi-channel form as
   well. Only genuinely multiplicative forms escape.
7. **An empirical instantiation (§7)** over ~10⁴ keys, in which the predicted
   below-unicity failure was observed and diagnosed.

---

## 2. The construction

Let *R* be a public record comprising *n* observable series
*x*₁(*t*), …, *x*ₙ(*t*), each of length *N*.

A **key** is a pair *k* = (*C*, *f*) where *C* ⊆ {1, …, *n*} with |*C*| = *m* is a
choice of channels and *f* ∈ *F* a functional form combining them. The keyspace
has size

  *K* = C(*n*, *m*) · |*F*|,  *H*(*K*) = log₂ *K*.

The **plaintext** is a modulation *s*(*t*) imposed on the channels in *C* such
that *f* applied to those channels exhibits structure absent under the null
hypothesis for *R*.

The **oracle** is a statistic *T*(*k*) computed from *R* under key *k*, together
with a null distribution for *T*. The plaintext is self-identifying in the sense
that *T* exceeds threshold if and only if *k* is correct — the analyst needs
nothing but *R* and *T* to know.

Three properties distinguish this from textbook keysearch:

**The ciphertext is not secret and not addressed.** *R* was collected by third
parties for unrelated purposes, and is public. There is no channel between the
party imposing *s* and the party searching.

**The scheme satisfies Kerckhoffs's principle in its strongest form.** Nothing is
hidden. The construction, the keyspace, the statistic and the record may all be
published without loss. Whatever difficulty exists lies entirely in the work of
enumeration.

**The keyspace is not uniform.** Whoever imposed *s* had to choose channels the
analyst plausibly measures, on records long enough to hold structure. This
concentrates the key distribution — see §5.2.

---

## 3. It is not proof of work: the payment is inverted

The construction is often described by analogy to proof-of-work. The analogy
misnames the primitive.

In hashcash and its descendants the **sender** performs work, to price a message
or prove commitment. Here the **receiver** performs the work, and the sender
performs none. The direction of payment is inverted, and with it the design goal:
the scheme does not ration access to the sender's resources, it conditions
legibility on the receiver's.

A more accurate description is a **proof-of-capability challenge** — a CAPTCHA
inverted. A CAPTCHA is a puzzle a machine is meant to fail. This is a puzzle only
a sufficiently resourced analyst can pass.

The distinction is not terminological. It determines what the scheme can enforce.
A proof-of-work scheme can impose a cost on a sender and verify it. This scheme
can impose nothing on anyone; it can only remain unrecognised until an analyst
spends enough. It is therefore self-enforcing in a weak sense — no party need
cooperate for it to hold — and coercive in no sense at all.

### 3.1 Relation to time-lock puzzles

The nearest technical antecedent is Rivest, Shamir and Wagner's **time-lock
puzzle** (1996), which requires a chosen amount of computation to open and whose
difficulty is set from *projected future hardware speed*, so that the puzzle
opens at a chosen future date. That is the same reasoning as setting a keyspace
size such that only an analyst above a given compute threshold can exhaust it.

The constructions differ in one load-bearing respect, and the difference matters
for what each can guarantee. A time-lock puzzle is **inherently sequential** —
its delay cannot be bought off with more machines, which is the property later
formalised in verifiable delay functions. **A keysearch parallelises perfectly.**
Our construction therefore cannot enforce a delay; it can only enforce a total
work requirement. An analyst with *P* processors finishes in 1/*P* the time.

This makes the scheme strictly weaker than a VDF as a timing primitive, and
strictly more practical as a capability filter: total work is a meaningful
threshold even when elapsed time is not.

---

## 4. Unicity distance is the multiple-testing threshold

### 4.1 The correspondence

An analyst enumerating *K* keys computes *K* statistics and must decide whether
the largest is evidence of a true key or the expected maximum of *K* draws from
the null. The standard remedy is a multiplicity correction: test at α/*K*, or
compare against the distribution of max *T* over *K* nulls.

This is Shannon's condition restated. Unicity distance is the ciphertext length
below which spurious keys yield plausible decryptions. The multiple-testing
threshold is the statistic level below which spurious keys fire. **They are the
same boundary, derived in two vocabularies** — one information-theoretic, one
frequentist.

The practical consequence is that the multiplicity correction in such a search
is not a conservative statistical convention that a more clever method might
avoid. It is an information-theoretic bound. No test can do better, because
below the unicity distance the information required to distinguish the true key
is not present in the record.

### 4.2 Deriving the scaling

Let *T* have null distribution with sub-Gaussian tails. The maximum of *K*
independent draws grows as

  E[max *T*] ≈ σ √(2 ln *K*).

A modulation of fractional amplitude *a* imposed on a record of length *N* yields
a statistic scaling as *a*√*N* against the same σ. Detection requires

  *a*√*N* ≳ σ√(2 ln *K*)  ⟹  **N ≳ (2σ²/a²) ln K.**

Since *H*(*K*) = log₂ *K*, this is *N* ∝ *H*(*K*): **the required record length is
proportional to key entropy**, which is exactly the form of Shannon's *U* =
*H*(*K*)/*D*, with 1/*D* playing the role of the per-sample information about the
modulation.

### 4.3 The consequence that should change reporting

If *N* < *U*(*K*), the search cannot conclude. It has not shown the plaintext
absent; it has shown that this record is too short to distinguish *K* candidate
keys, for any test whatsoever.

We propose that such an outcome be reported as **void** rather than **null**, and
that the two be counted separately:

| outcome | meaning |
|---|---|
| **negative** | *N* ≥ *U*(*K*); the key is not in the searched space at the stated amplitude |
| **void** | *N* < *U*(*K*), or no calibrated null exists; the search is uninformative |

Reporting a void as a negative overstates coverage and, in a setting where
searches accumulate, corrupts any estimate of how much of the keyspace has been
eliminated.

---

## 5. Keyspace expansion is cheap in data, expensive in compute

### 5.1 The asymmetry

From §4.2, enlarging the keyspace by a factor of ten costs a factor of
ln(10*K*)/ln(*K*) in record length — for *K* ~ 10⁴, about 25%. It costs a factor
of ten in computation.

This inverts a common intuition. Analysts frequently restrict a search space to
"preserve statistical power," reasoning that fewer tests means a laxer threshold.
The scaling says the saving is logarithmic and the loss is the risk of excluding
the true key entirely. **Where the record is fixed and compute is elastic — which
is the usual modern situation — the correct move is to search the larger space.**

Two caveats attach. The ladder below counts *nominal* keys; the **reachable**
keyspace may be smaller, since a key whose channels share no overlapping record
cannot be tested — in the instantiation of §7 the reachable pair space was 88% of
nominal and the reachable triple space 72%. The correct *K* for a threshold
calculation is the number of keys actually tested, not the number enumerable.

An illustrative ladder over *n* = 30 observables with |*F*| = 4 forms:

| tier | keys *K* | *H*(*K*) bits | record length, relative |
|---|---|---|---|
| pairs | 4.35×10² | 10.8 | 1.00× |
| triples | 4.06×10³ | 14.0 | 1.30× |
| quadruples | 2.74×10⁴ | 16.7 | 1.56× |
| triples over 60, 10 forms | 3.42×10⁴ | 18.4 | 1.71× |

A sixty-fold keyspace expansion costs 71% more record.

### 5.2 Weak keys, and why brute force is the wrong first attack

The keyspace is not uniform. Whoever chose *C* was constrained to channels the
analyst plausibly measures, at cadences and durations sufficient to carry
structure. The key distribution concentrates on the best-instrumented and most
obvious observables.

In cryptographic terms these are **weak keys**, and the standard response applies:
*no password cracker begins with brute force; it begins with a wordlist.* An
ordered search over a prior-weighted enumeration reaches a probable key far
sooner in expectation, at no cost in coverage provided the full enumeration is
eventually completed.

We stress the limit of this argument. The prior is the analyst's reasoning about
a designer's choices and is not testable. It justifies search **order**, which
costs nothing if wrong. It does not justify restricting the search **space**.

---

## 6. A structural constraint on multi-channel carriers

A natural question is whether a modulation can be hidden from all pairwise
analysis while remaining visible to a higher-order form. For linear forms the
answer is no.

**Proposition.** Let a modulation *s*(*t*) be imposed with coefficients *c*ₖ on
channel *k*. Consider the linear four-channel form with coefficients
(1, −3, 3, −1) — the third finite difference — and the pairwise forms
*x*ᵢ − *x*ⱼ. If the modulation is invisible to every pairwise form, it is
invisible to the four-channel form.

*Proof.* The pairwise form on (*i*, *j*) picks up (*c*ᵢ − *c*ⱼ)*s*. Invisibility
to every pair requires *c*ᵢ = *c*ⱼ for all *i*, *j*, hence all *c*ₖ equal to some
*c*. The four-channel form then picks up (*c* − 3*c* + 3*c* − *c*)*s* = 0. ∎

The result extends to any linear form whose coefficients sum to zero, which
includes the finite differences and cross-differences commonly used as
higher-order statistics. **Such forms are not independent carriers**; they are
differences of pairwise information.

Escaping the proposition requires a genuinely multiplicative form. A four-way
phase alignment, living in the trispectrum, leaves every pairwise
cross-correlation untouched and is therefore a real higher-order carrier.

The practical import for search design is that enumerating linear higher-order
forms enlarges *K* — and by §4.2 raises the detection threshold for every key —
while adding no reachable plaintext. It is a strictly bad trade.

---

## 7. Empirical instantiation

We instantiated the construction on 30 public observational time series,
enumerating pairs and triples under four functional forms.

The reachable keyspace is smaller than the nominal one, and the distinction
matters for §4.2: a key whose channels have no overlapping record is not a key
the analyst can test at all. Of 435 nominal pairs, 382 had the required
2,000-sample overlap, giving **1,074 completed tests** against 60,000 circular
shifts each. Of 4,060 nominal triples, 2,934 were reachable and 2,749 completed,
giving **10,996 tests** — 117 million statistic evaluations in 26 hours on one
88-core machine.

The pair search completed and returned negative at a stated amplitude. **The
triple search returned void, in the precise sense of §4.3**, and the diagnosis
is instructive: no calibrated null distribution could be constructed for the
higher-order statistics. Three candidate nulls were built and all three failed,
each because it required a false modelling assumption about the underlying
process:

| null | assumption | failure |
|---|---|---|
| circular shift | pairwise structure irrelevant to a three-body statistic | destroys it; ninefold tail excess |
| common phase screen | series are stationary | they are not; real data at median *z* = 45 |
| envelope-preserving | the nonlinearity is pure amplitude modulation | absorbs the signal; 4.25× over-firing |

Additionally, and in agreement with §6, the two linear higher-order forms were
swamped at *z* ≈ 50 while the single multiplicative form sat near its null at
median *z* = 5.1 — the predicted signature of forms that are not independent
carriers.

The full search, its code and its data are public [1].

---

## 8. What is not claimed

- **No hardness result.** The mapping to a cipher is structural. No security
  proof, reduction or hardness assumption is transferred, and none is implied.
  The keyspace arithmetic and the √(2 ln *K*) threshold stand on their own.
- **No claim of novelty for the components.** Unicity distance is Shannon's.
  Self-identifying plaintext is the standing assumption beneath every brute-force
  attack. Time-lock puzzles are RSW's. The claimed contribution is the
  correspondence of §4, its scaling consequence, the void/null distinction, and
  the proposition of §6.
- **The weak-key prior is not testable.** It orders a search; it does not bound
  one.
- **The empirical instantiation detected nothing.** It is reported as an
  existence proof for the failure mode of §4.3, not as evidence for any
  hypothesis about its subject matter.

---

## 9. Open questions

1. Does the correspondence of §4 extend to keyspaces with non-trivial structure
   — where keys are related, so that the *K* tests are dependent? The
   sub-Gaussian maximum is then an over-estimate and an effective *K* smaller
   than the nominal one may apply.
2. Is there a construction with self-identifying plaintext that is also
   inherently sequential, recovering the VDF property lost in §3.1?
3. Can the void/null distinction be given a sharper operational test than "no
   calibrated null exists" — ideally one computable before a search is run?

---

## References

[1] Griffin, R. *A combination-space search for embedded technosignatures in
solar and heliospheric archives.* 2026. Code and data:
https://github.com/dxtrainc/behavioral-seti

[2] Shannon, C. E. *Communication theory of secrecy systems.* Bell System
Technical Journal **28**(4), 656–715, 1949.

[3] Rivest, R. L., Shamir, A., Wagner, D. A. *Time-lock puzzles and
timed-release crypto.* MIT/LCS/TR-684, 1996.

[4] Boneh, D., Bonneau, J., Bünz, B., Fisch, B. *Verifiable delay functions.*
CRYPTO 2018.

[5] Back, A. *Hashcash — a denial of service counter-measure.* 2002.

[6] Kerckhoffs, A. *La cryptographie militaire.* Journal des sciences
militaires **IX**, 5–38, 1883.

[7] Freitas, R. A. Jr. *Xenology: An Introduction to the Scientific Study of
Extraterrestrial Life, Intelligence, and Civilization.* 1st edn., §24.2.4
"Alien Message Contents", 1979.

[8] Busch, M. W., Reddick, R. M. *Testing SETI message designs.*
arXiv:0911.3976, 2009. Also in Vakoch, D. A. (ed.), *Communication With
Extraterrestrial Intelligence*, 2011.

[9] Freudenthal, H. *Lincos: design of a language for cosmic intercourse.*
North-Holland, 1960.

[10] Heller, R. *Decryption of messages from extraterrestrial intelligence
using the power of social media — the SETI Decrypt Challenge.* International
Journal of Astrobiology **18**(4), 296–303, 2019. arXiv:1706.00653,
doi:10.1017/S1473550417000568.

[11] Ball, J. A. *The zoo hypothesis.* Icarus **19**, 347, 1973.

---

*Draft for IACR ePrint. All references verified against the originals.
ePrint requires PDF; see the LaTeX source alongside this file.*
