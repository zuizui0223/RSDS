# First original GPS observation-process audit: study measurability is not interval throughput

**2026-11-15 | RSDS PR #27 | Post-outcome methodological diagnostic | Original source data replayed, not reclassified.**

Original code/read-only first-result:
[35-source source-native observation funnel](../research/RSDS_35_SOURCE_GPS_OBSERVATION_FUNNEL_FIRST_RESULT_20261115.json).
Original Git blob `49dcf5c798f8eee788c8001906a1fe0e2c195086`,
original JSON SHA-256 `e7c923ad331c8b7c5a8cdfb1dc11d3f510961ce2d145d79c3c0b4c0d16ae8acf`.
The first full [35 original source GPS Action run](https://github.com/zuizui0223/RSDS/actions/runs/43697404274)
has **37/37 jobs SUCCESS**: original numeric and synthetic-integrity
tests, 35 original-source private GPS or immutable HOLD jobs, final
35-source source-equal aggregate. No new ecological q outcome or original
species/source roster was substituted.

A companion source-group
[q-measurability rank sensitivity](../research/RSDS_35_SOURCE_Q_MEASURABILITY_RANK_DIAGNOSTIC_FIRST_RESULT_20261115.json)
is computed ONLY from the frozen aggregate; no GPS reruns.
[Source-level rank CI](https://github.com/zuizui0223/RSDS/actions/runs/43698864972)
also succeeded.

## What was measured

**One original repository study-ID = one source replication unit**.
Among the original 35 source studies / 50 registered native source-taxon
rows, only **25 previously parsed source studies / 34 registered native
rows** could be reprocessed under the same original MD5 files. Nine
nonhuman original file/sourcewide HOLDS and the one human-location privacy
HOLD were preserved without GPS access. Original source q remained
fully estimable in **10/35**, daylight in **8/35**, and q available
for at least one native in **13/35**; these are unchanged old endpoints.

The new audit **only** partitions the *identical original GPS
movement step admission* into nested counts. Segment records below
are **not independent study replicates**.

| Gate | Original GPS intervals left | Conditional admission rate |
|---|---:|---:|
| Adjacent original hour-selected fixes | 1,502,093 | — |
| Actual elapsed 0.5–2.0h | 1,218,362 | 81.11% |
| Same astronomical solar date and daypart across endpoints and midpoint | 1,068,433 | 87.69% |
| Finite and <=60 km/h great-circle GPS step speed | 1,066,632 | 99.83% |

The final accepted segment share is 71.01% of original adjacent
source steps. **Time opportunity** removes 283,731 adjacent pairs
(18.89%) and astronomy/daypart continuity removes a further 149,929
(12.31% of those time-eligible). The original 60-km/h geometry
filter removes only 1,801 (0.17% of solar-continuous pairs).

These data identify **where the source GPS processing rule rejects
observations**; they do not identify how each animal would behave
under counterfactual GPS sampling, or prove that source habitats or
biological velocity do not differ.

## Source units, not GPS events: distinguish complete from incomplete q

The following figures average the original native fractions *within each
source* first, then take **source-equal** medians in original
q-complete studies vs previously q-incomplete but GPS-parsed studies.

| Source-level score median | q-complete, 10 source studies | q-incomplete, 15 parsed source studies |
|---|---:|---:|
| Original time 0.5–2h gate fraction | 0.9049 | 0.8121 |
| Original astronomical continuity among time-passing intervals | 0.8686 | 0.8869 |
| Original speed validity among solar-continuous intervals | 1.0000 | 1.0000 |
| **All interval gates passed / adjacent** | **0.7081** | **0.7077** |
| UTC even-hour fix share | 0.5040 | 0.4999 |

One previously parsed q-incomplete source had **zero original adjacent
intervals**, so the rate scores have an undefined denominator. The
q-incomplete group medians above use **14/15** sources, not zero-fill.
There are no inferential p-values or source-matched causal comparisons.
The source GPS observation gates do not independently explain q
availability, because original animal-year/month and native tag
support are stricter *non-collapsible* conditions.

### Direct counterexamples from original tagged animals

- *Capra hircus*: 30,233 adjacent pairs; 26,562 time-eligible,
  20,539 final qualified intervals, original 12/13 tagged animals
  with >=10 qualified DAY/NIGHT dates in at least one year; source
  native q estimable. Its hour-selected source GPS fixes were 100%
  on even UTC hours.
- *Aquila nipalensis*: 137,934 adjacent pairs; 123,888 time-eligible,
  112,847 final qualified intervals (81.8% overall), **but only 2/21
  original tagged animals** crossed the original annual paired-date
  eligibility threshold. Its native q remains UNAVAILABLE, not zero.
  Hour-selected fixes were approximately 49.7% even UTC.
- *Chelonia mydas* and *Eretmochelys imbricata* in one original source:
  the former has only 2 total original tags, while the latter has 21
  qualified q tags. Dropping the smaller-sample registered species
  would improperly turn a source-incomplete comparison into a
  false independent source-complete result.

The observed **source-level thresholded independent animal support**
therefore cannot be replaced with number of fixes, number of
valid segments, raw data-package animal metadata or UTC phase balance.

## Post-outcome rank diagnostic: measurement fraction is not a validated predictor

For each technical source score, compare the score distribution between
the **10 original source-complete q studies** and the **15 parsed but
source-incomplete studies** by the empirical source-pair AUC.
This is only a descriptive statistic of the already-selected 25 source
studies; it has **no out-of-source training/validation** or calibrated
world-level interpretation.

| Score, larger = more likely source-complete | Source-pair AUC | Eligible original source studies |
|---|---:|---:|
| Original adjacent GPS interval count | 0.680 | 25 |
| 0.5–2.0h time-interval acceptance fraction | 0.707 | 24 |
| **Final accepted interval fraction** | **0.543** | 24 |
| Original UTC even-hour share | 0.650 | 24 |
| Number of q-qualified original tags, outcome-definition component | 0.913 | 25 |

The high tagged-animal AUC **does not count as an independent
prediction**: >=5 eligible animals for every registered native taxon
is literally part of the source q-complete endpoint definition.
Its apparent discriminability is tautologically influenced by
how q availability is defined. In contrast, interval acceptance
alone barely distinguishes source-complete q in the first sample.

The time-only fraction gives descriptive AUC 0.750 in the original
first wave and 0.667 in the later second wave. The combined accepted
interval fraction gives AUC 0.500 and 0.524 within the same waves.
Those source waves were not a fully prospective randomized
generalization experiment (Wave2 quotas chosen after Wave1 outcomes).

## Implications for meta-analysis

1. **Do not count GPS fixes as source replicates.** One original
   field-study DOI/repository ID is one source. Native species and
   original tagged animals are nested inside it. There are 35
   preselected source groups, only 25 originally parsed and only
   10 original source-complete q. A source with 112k qualified
   intervals can be q-incomplete.
2. **Technical observation cadence is not the only gate.** Elapsed
   interval availability, solar date/daypart pairing, geometric
   plausibility and minimum independently observed animal-year
   support are separate, and q is measured only after all.
   The geometry filter and solar continuity themselves reflect
   actual source geography/animal behavior; they are not purely
   equipment-generated missingness.
3. **Source selection is informative.** The prior 10-respondent
   composition-adjusted q mean +0.07568, and its earlier 4-respondent
   value -0.04779, cannot be extrapolated as a global night/day
   biological direction. Event acceptance/UTC duty-cycle adjustment
   alone will not restore missing original independent-tag samples.
4. **Better future inference** requires prespecifying response-blind
   source tracking design metadata, including independent tag count,
   original sampling schedule and valid daypart opportunity. Then
   evaluate source-complete q and the sampling process on entirely
   untouched new original studies, with honest study-level
   uncertainty and no retrospective threshold tuning.

The original 48-source externally signed sunlight-direction
**negative replication remains NOT REPLICATED**; the 12 legacy
ERA5 provider raw-JSON-source identity issues remain **HOLD**. Neither
the observational GPS time niche nor the astronomical daylength
association identifies causal photoperiodic selection.

**Methodological contribution:** an auditable *source measurement
selection decomposition*, validated against old independent-tag
q, that shows why source GPS observation throughput is not sufficient
to produce source-comparable ecological effect sizes. This is
more defensible than a broad, unjustified worldwide biological
mean estimated from the ten responding studies.
