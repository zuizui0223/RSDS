# Source selection is part of the ecological estimand: first 16-source original Movebank GPS wave

**2026-11-15 first biological event run** — source-preallocated before new GPS outcomes.
Original Actions run [#43690283579](https://github.com/zuizui0223/RSDS/actions/runs/43690283579):
all **16/16 original-source jobs SUCCESS**, merged first source outcome
**SUCCESS**. Important: job success includes explicit scientific
HOLDs, not 16 usable biological effects.

Original first-run atlas artifact:
`rsds-437-wave1-native-q-multi-source-meta-atlas-v1`
(artifact **14151799059**, raw JSON SHA256
`8e023272683b2f46ab30c35ec62632eab17b5e0aef89c49f774d03d98f114c6c`).

This manuscript note accompanies the main
[437-package sampling and source methods](RSDS_SYSTEMATIC_437_SOURCE_FRAME_WAVE1_ORIGINAL_GPS_20261115.md).
Results are **descriptive and post-historical-RSDS**, not a globally
representative or causal estimate. The wave's selected source roster,
original provider file identities, sample weights and admission rules
were SHA-frozen before any new wave animal GPS event was opened.

## Causal/measurement identification problem

The first observable outcome is the **ability to identify a comparable
ecological effect**, before the effect itself. Let a registered source
study `s` have complete native-taxon effect `Y_s`, and let `R_s=1`
when every registered native has sufficient, schema-matched original
GPS day/night contrasts under the fixed q rules.

When `R_s=0`, `Y_s` is not zero; it is **not identified** from this
source under the original method. Selection from source package
metadata, original provider-file budget, GPS schema, native species
support, temporal cadence and calendar overlap define several
successive measurement filters. Meta-analytical synthesis conditional
only on `R=1` targets a *different selected sample* from the
original frozen 120 candidate-package frame.

## Full accounting: no positive-source cherry-picking

| Source screening/qualification | Number | Denominator / limits |
|---|---:|---|
| Frozen public packages in repository | 437 | Metadata packages, not independent ecologies |
| GPS/GNSS sensor source packages | 259 | File schema not validated |
| CC0, prior-excluded, source-ID dedup, metadata threshold candidates | 120 | 109 single-native + 11 multi-native |
| Deterministically sampled first wave | 16 | 8 single + 8 multi |
| Frozen provider-file UUID+MD5 under 200MB | 13 | 3 exact source-file HOLDS |
| Original GPS source schema loaded without source-wide HOLD | 11 | 2 schema/parse source HOLDS |
| >=1 source native taxon q estimable | 6 | 2 of these are partial multi-native |
| **Every registered source native q estimable** | **4** | 2 single + 2 multi; actual eligible meta DOI |
| Every registered source native daylength slope estimable | 3 | 2 single + 1 multi |
| All originally registered multi-native taxa with common-month covariance | **0** | Do not impute or claim shared temporal niches |

The original 16-source status ledger is:

- 3 source-file metadata/byte-budget HOLD
- 2 original GPS-wide schema/HOLD
- 5 original GPS files read but *no* supported q effects
- 6 source studies with at least one native q, of which 2
  retain unestimable registered native taxa

This is a source-level **11/16** original GPS parsing admission
and **4/16** all-native q endpoint, not a 4/13 result
obtained by deleting earlier file/HOLD sources. The **27 original
registered native study labels** stay in the denominator.

## New observational ecological effects (not a common global sign)

Four complete-source examples, each one original repository study:

| Source / registered native taxa | Complete source equal-native q | Support |
|---|---:|---|
| Feral goats (Capra hircus) | **−0.16936** | 12 q-tagged animals, one source |
| Grey seals (Halichoerus grypus) | **+0.04221** | 15 q-tagged animals, one source |
| Four grazing-goose taxa | **−0.17446** | 4/4 native q estimable, one source |
| Red and Arctic foxes (Vulpes vulpes, V. lagopus) | **+0.39162** | 2/2 native q estimable, one source |

Geese provide **four** estimates within the same field source, *not*
four independent regions. The two fox species also comprise **one**
source. Co-listed taxa and source-day records do not prove physical
sympatry, competition or shared experienced solar environment.
Original source q is observed GPS straight-line displacement speed
relative to recorded day/night, not latent activity duration. Other
registered multi-native sources have partial q or no q, all kept as
unavailable rather than zero.

All first-wave source-native `q~astronomical daylength` estimates are
also observational within animal-year slopes; a time/calendar
association cannot be attributed causally to photons, temperature
or seasonal environment without stronger comparisons.

## Why the source-selecting probability matters

The first-wave sampling rate deliberately differs by native metadata
stratum:

- single-native: 8 selected among 109 eligible public packages
- multi-native: 8 selected among 11 eligible public packages

An equal 16-source mean would give multi-native studies a vastly
different representation from the fixed 120-package metadata frame.

Among *complete respondents only*, `q` source mean was
**−0.06357** for the two single-native studies and **+0.10858**
for the two multi-native studies. The descriptive 109/120 vs 11/120
poststratified **respondent-only** contrast is **−0.04779**.

This is NOT a design-unbiased estimate of all 120 candidates: 12
selected sources lack complete native `q` and the source's ecological
observation process is likely informative. The all-16 selected wave
would have a poststratified q within approximately **[−0.14280,+0.28177]**
*only if* the 12 missing source-level q values were forced
hypothetically into the range observed among the four complete
source effects **[−0.17446,+0.39162]**.
That is a **conditional tipping sensitivity exercise**, not an
identified outcome bound, formal confidence interval, or world
estimate. Without a bound on missing effects, even its sign is
unconstrained.

Daylength complete-source responders are 3/16. The same nonresponse
caution applies. Do not combine these responders or the snapshot's
prior 12 source group coefficients as if they represented an
unbiased, independently sampled global study population.

## Joint native covariance: important methodological negative result

Although the new source wave deliberately oversampled multi-native
metadata candidates, **0/8** selected multi-native study slots
passed the fixed common native-calendar month support requirement.
Hence no original GPS joint multi-native covariance matrix was
certified **in this wave**. Failure may arise from taxon-uneven
GPS schedules, native support, calendar overlap, or schema; absence
of joint covariance does not prove biological synchrony is absent.

The earlier selected *four-source* same-calendar
[PR #24](https://github.com/zuizui0223/RSDS/pull/24) remains
separately valid under its different data/selection architecture.
The zero result for this wave **must not reclassify** those already
frozen positive original covariance results.

## Main next method instead of more speculative moderator regressions

Model and report **source measurability as a selection process**
rather than a nuisance drop:

`P(R_s = 1 | original publication metadata, taxa, reported animals,
reported positions, source-file budget, raw schema, observation cadence)`.

Only an outcome-blind, source-frozen full frame and explicit
nonresponse diagnostics permit principled extrapolation. With just
16 selected source slots, it is premature to fit high-dimensional
ecological moderators or claim predictor-level inference for
nonrespondents. First improve availability of observed contrasts
without selectively relaxing calendar or event thresholds, e.g.:
new pre-specified alternative source-first estimands explicitly
separate from the earlier q endpoint; use all 120 metadata candidates
and expand source recruitment with prespecified wave quotas.

The scientific result is both positive and limiting:

**Positive:** A robust source-first, DOI/study-ID-aware analysis
can retrieve native GPS time-niche contrasts and preserve exact
source-level provenance over an outcome-blind new-source wave.

**Limit:** Public GPS metadata-level eligibility substantially
overpredicts the number of comparable source-complete biological
effects. Source-selection and observation processes must be treated
as explicit parts of cross-study meta-analytic generalization.

All older 48-source irradiance-direction negative external replication
remains **NOT REPLICATED**; old 12 ERA5 raw-version/HOLD is **unchanged**.
No new evidence for causal solar response or solar compass is claimed.
