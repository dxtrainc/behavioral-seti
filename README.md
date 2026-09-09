# Behavioural SETI — analysis code

Code and intermediate results for *A combination-space search for embedded technosignatures in
solar and heliospheric archives* (R. Griffin, Dxtra Inc.) and its companion,
*Set and forget: initial conditions as a technosignature*.

Every number in the papers is produced by a script here. **The table below is the audit trail**:
paper claim → script → result file. Nothing is quoted in the papers that is not regenerable.

## Requirements

    pip install -r requirements.txt        # numpy, scipy, pyhdf

`pyhdf` is needed only for ACE SWICS, whose level-2 products are HDF4 Vdata tables.
Runtimes below are for 80 cores.

## Data

No data is committed. Everything is public; `acquire/` fetches it.

    python acquire/pull_chan.py            # ACE SWICS, sunspot area, WSO, Ca II K, LASCO
    python acquire/getpsr2.py              # ATNF psrcat -> psrcat.json

GOES-R EXIS, OMNI2, SILSO, F10.7, NMDB and LISIRD products are listed in the paper's data
availability statement and were retrieved from the archives named there.

## Audit trail

| Paper claim | Script | Result |
|---|---|---|
| §3.3, §4.8 — 13-channel sweep, 229 tests | `search/sweep13.py` | `results/sweep13.json` |
| §4.8 — 30-channel sweep, 1,074 tests, 22 survivors, **0 of 510 across the boundary** | `search/sweep30.py` | `results/sweep30.json` |
| §5.2 — tail-fitted null, pair space | `search/sweepT.py --mode pairs --shifts 10000` | `results/sweepT_pairs_zeus.json` |
| §5.2 — tail fit validated against counting | `validate/validate.py` | printed |
| §4.9 — triple sweep, 10,996 tests, 198 screened | `search/sweepT.py --mode triples` | `results/sweepT_triples.json` |
| §4.9 — direct-count confirmation, 93 of 198 | `search/confirm.py --shifts 500000` | `results/confirm_triples.json` |
| §4.9 — **form-level control that voids the sweep** | `validate/fptest.py` | printed |
| §4.2 — fast-band limits | `search/fastsearch.py`, `search/euvsfast2.py`, `search/nmfast3.py` | printed |
| §4.2, Fig 2 — **injection recovery curves** | `validate/inject2.py` | `results/injection.json` |
| §4.6, Fig 3 — p-mode comb, 135.1 / 135.0 µHz | `search/pmode.py`, `search/pmode_cutoff.py` | — |
| Fig 3 — comb curves for plotting | `figures/figdata3.py` | `results/figdata_comb.json` |
| §5.7 — 53 unreachable pairs, all permanently blocked | `validate/span30.py` | `results/reach30.json` |
| Companion §3.1 — population back-extrapolation | `search/backex2.py` | printed |
| Companion §3.2 — subset ladder search | `search/subset.py` | printed |
| Companion §3.3 — Shklovskii + Galactic correction | `validate/shk2.py` | `results/psrcat_shk.json` |
| Companion §3.3 — subset search on corrected Ṗ | `search/subset_c.py` | printed |

## Two things a reader should check first

**The p-floor against the significance threshold.** Counting exceedances among *M* circular shifts
floors *p* at 1/(*M*+1). If that floor sits above α/*N*, no test can survive and a clean-looking null
means nothing. This happened here: the 30-channel sweep was first run at 6,000 shifts against a
threshold of 4.66×10⁻⁵ and returned "0 of 1,074 survive", which is reported in §5.2 as void. `confirm.py`
refuses to run when the floor is inadequate.

**Whether a detector can fail.** Three statistics in this work could not, and were caught only by
injection or by a form-level control: a log-ratio autocorrelation that sat at 0.954 for data and
surrogates alike; a residual-product form sitting 4.5σ below its own null; and a phase-randomised
surrogate used against a statistic computed from the power spectrum it preserves. `validate/fptest.py`
is the template for testing a new form before trusting it.

## Known corrections

Recorded because they are load-bearing, and because the papers report them rather than hiding them.

- The §4.2 limits were originally described as 99%-confidence figures established by injection. They
  are analytic **threshold-crossing** amplitudes, recovered about half the time, and no injection had
  been run. `validate/inject2.py` supplies the real curves; the cosmic-ray limit holds up at 95%
  recovery, Lyman-α and soft X-ray recover 1% and 16% at their quoted values.
- The continuum-normalised power is **not** exponentially distributed (mean 1.44, tail ~4×10³ heavier
  than Exp(1)), so the analytic threshold is itself optimistic.
- A re-implementation of the mode-comb detector returned 127 µHz where the original gives 135, from a
  4001-minute median detrend and a linear continuum in place of a 60-minute mean and a log-space one.
  Use `search/pmode_cutoff.py`; `figures/figdata3.py` reproduces it with a strided median for speed.
- Mesh numerics in an early draft did not cohere; they are now derived from the link budget.
- A cost estimate was wrong by eight orders of magnitude by costing probe buses rather than modulators.

## Licence

Code (`acquire/`, `search/`, `validate/`, `figures/`) — **MIT**, see `LICENSE`.
Derived results (`results/`) — **CC BY 4.0**, see `LICENSE-DATA`.

The underlying observations are not relicensed: they belong to NOAA, NASA, LASP,
SILSO, NMDB, ATNF and the other archives listed in `LICENSE-DATA`, and should be
cited to them. No observational data is committed to this repository.

## Software

Written and run with Claude Code (Anthropic). Scope and documented failures are set out in the paper's
Software section. Responsibility for every number rests with the author.

## Manuscript changelog

Moved here from the manuscript at the reviewer's request: a changelog belongs with the
repository, not in a paper submitted to a journal.

### Round 4

Changed since the last review. §3.6 is new and is the largest addition: a designer-priority 
table ranking fifteen channels by what a sender would choose, against what we can resolve and 
whether we have searched there. It sits before the Results deliberately, so its ordering is a 
prediction rather than a description of what was found, and §1.4 now carries its uncomfortable 
consequence — that rows 2&ndash;7 are capable and unsearched, so this paper has not yet 
tested the framework's own best guess. §4.12 gained the cross-viewpoint test as the 24th 
search and has since been extended to the fast band on 1,879 days of MAVEN L2: null, with the 
binding limit set by MAVEN rather than GOES and 580× weaker, and with the band capped at 6 h 
because Earth and Mars see solar rotation at different synodic periods. §4.13 reports lunar 
laser ranging outside the tally, with the reason — that no dynamical method can reach a 
2.1×109 kg modulator — stated before the nulls. §2.4 gained the neutrino gate as content 
behind a later gate and §2.5 the coupling calculation behind §4.13. §5.5 proposes an untried 
neutrino-line search. Five references added, all verified against ADS and the publishers; the 
existing Learned et al. 2008 entry is a different paper from Learned, Pakvasa & Zee 2009. 
Earlier in this round: §4.10 replaced the note that dismissed the self-keyed excess on the 
sign of its correlation, §4.11 answered the annual-geometric-term objection, the stale 3.7% 
false alarm became the measured 6.3%, and the outcome table gained the false-alarm-and-power 
row it had been missing.