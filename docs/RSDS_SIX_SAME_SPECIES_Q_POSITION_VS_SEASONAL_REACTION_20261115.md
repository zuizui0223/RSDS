# RSDS biological result: day–night movement niche position can persist while daylength association changes sign

**15 November 2026 | Source-native original Movebank GPS; post-outcome descriptive synthesis.**

## Ecological question

> Can the **direction** of an animal species' recorded night-minus-day displacement persist across archived studies while its **seasonal daylength association** differs?

This is more informative than a universal "nocturnal versus diurnal" average. It explicitly separates two biological descriptors: (i) a source-tag-equal mean night-minus-day GPS straight-line step displacement contrast (q), and (ii) the association of within-tag/year (q) with astronomical daylength, per hour. Neither is equivalent to the true duration of animal activity or a controlled photoperiod response.

## Verified results

The previous archive-only 63-DOI/88 source-native atlas included six identical species with at least two original DOI-specific (q) estimates. After the prior outcome-guided species choice, the **six new same-species study IDs were fixed before the new (q) outcomes**. Four passed the original provider bitstream UUID, bytes ≤200 MB and MD5 gates; two were permanent file/source-identity HOLD. Four original files were examined under the unchanged original GPS estimator, **three supported (q)**, and the other original source (*Branta leucopsis*, Barents Sea) had insufficient supported source-native (q), not a measured zero.

All **three** supported new (q) estimates reproduced the previously specified signed direction with their *original per-source tagged-animal bootstrap nominal 95% percentile intervals excluding zero*. This is a **conditionally prospective source-outcome comparison**, **not** a new globally preregistered species test: the species were selected after the earlier results; the new fox original GPS was previously analyzed for solar direction; and absence of earlier source matches in eight archived files is not proof of wholly unexposed GPS.

| Species | Prior q source-native points (independent DOIs only) | New original study q | New original tagged-animal 95% interval | Nominal q direction |
| --- | --- | ---: | --- | --- |
| *Anser albifrons* | −0.12099, −0.12730 | **−0.24516** | [−0.29789, −0.19018] | negative repeated |
| *Odocoileus hemionus* | −0.02185, −0.00446 | **−0.02882** | [−0.03348, −0.02422] | negative repeated |
| *Vulpes vulpes* | +0.18975, +0.29043 | **+0.24182** | [+0.17628, +0.31417] | positive repeated (GPS previously examined for another endpoint) |

**The more important newly highlighted ecological contrast is the lack of a common seasonal association.**

### 1. *Anser albifrons*: q direction remains negative; daylength slope is nominally opposite

Earlier agricultural grassland goose source (DOI `10.5441/001/1.fk899541`): source-native (q) negative, and (q)–astronomical-daylength slope **−0.02653/h** with original tag bootstrap 95% **[−0.04778, −0.00477]**. The separate firework-response goose source had negative q but unavailable daylength slope. New Pannonic migration source `10.5441/001/1.46b0mq21` (only five q-eligible tagged animals) has (q=-0.24516) with negative CI and daylength slope **+0.02393/h [ +0.00871, +0.04020]**. These are opposite *within-source* nominally supported associations; no cross-source slope difference CI was estimated.

### 2. *Odocoileus hemionus*: q direction remains negative in three source archives; seasonal slope changes sign

Earlier Utah source `10.5441/001/1.711`: q **−0.02185 [−0.02321,−0.02047]**, slope **+0.005619/h [+0.004388,+0.006921]**. Earlier Wyoming source `10.5441/001/1.703`: q **−0.00446 [−0.00578,−0.00309]**, slope **+0.005538/h [+0.004102,+0.006973]**. New Wyoming Sublette source `10.5441/001/1.708`: q **−0.02882 [−0.03348,−0.02422]**, slope **−0.003920/h [−0.008944,−0.000331]**. Its 133 independent q-eligible tags give a reasonably large *within-source* estimator, but even here different study ID is not proof of disjoint animals or location relative to the earlier Wyoming tracking study.

This is a particularly sharp *observational* example that an approximate conservation of the mean day–night GPS movement direction need not imply conservation of how that difference covaries with season.

### 3. *Vulpes vulpes*: strong q direction consistency; seasonal slope not resolved in the new source

Earlier UK wet-grassland source q **+0.18975** and coexistence source **+0.29043** have positive original tagged-animal CI support. New Australian fox source `10.5441/001/1.72hh609t` has q **+0.24182 [+0.17628,+0.31417]**. However its daylength coefficient +0.05417/h has interval **[−0.02504,+0.13625]**, including zero. One earlier source slope was positive, the other included zero. No reliably conserved seasonal direction is identified. Crucially, this new fox source's original GPS was **already opened for the solar-heading endpoint**, so it cannot be called untouched-GPS external validation.

## Complete six-source accounting

- Source q supported: **3/6**, all nominal same-sign source bootstrap intervals.
- Source file HOLD: **2/6** (original source size >200 MB and earlier source-file ambiguity, respectively).
- Original parsed/MD5-confirmed but source q unavailable: **1/6** (*Branta leucopsis*, Barents Sea).
- Three studies have source-native daylength slope estimates; **two species** (*Anser albifrons*, *O. hemionus*) have at least one earlier originally supported daylength slope of the opposite sign.
- *Ovis canadensis* has zero eligible NEW study in the original fixed catalog, not a negative replication; *Anas platyrhynchos* has new original source-file HOLD, not a q sign.
- The old 48-original-source held-out sunlight/radiation→direction sign verdict remains **NOT REPLICATED**, and old 12 ERA5 provider raw-SHA source cases remain **HOLD**.

## Scientific interpretation and ceiling

**Provisional ecological hypothesis**: A species' night/day movement niche **position** can be relatively persistent in sign across source contexts, while its **seasonal modulation** is more contingent on local conditions and source-specific ecology. Ecological interpretation should ask which species maintain a relative day/night position and which environmental/seasonal contexts modulate its response.

**Not established**: opposite source-level slopes as a formally significant between-population interaction; phenotypic plasticity, adaptation, causality of light or daylength; independence of individual tags/populations across two Wyoming studies; cross-study sampling-clock harmonization; global prevalence of directional stability. The q and seasonal slopes were derived using the same declared *form* of original GPS estimator, but different actual source seasons, latitude/solar ranges, animals and recording duty cycles can account for part of the contrast. The new interpretation is post-result and multiple species/contexts were examined. Nominal source-specific bootstrap intervals are not corrected for this search.

Before asserting generality or a true "conserved position, labile response" mechanism, quantify original study-season/daylength support and source exposure/animal overlap. A next **separate** original-data sensitivity route must predefine comparable daylength ranges and within-source observation cadence restrictions; it cannot turn the post-hoc findings here into a preregistered confirmation.

## Immutable sources and execution

- [Previously frozen 63-source exact source-native atlas](../research/RSDS_ORIGINAL_SOURCE_NATIVE_SAME_SPECIES_CROSS_STUDY_FIRST_RESULT_20261115.json)
- [Six follow-up source IDs selected and frozen pre-file/pre-q](../research/RSDS_SIX_SAME_SPECIES_FOLLOWUP_ORIGINAL_SOURCE_ONLY_FIRST_RESULT_20261115.json)
- [Old solar GPS exposure audit, limited to eight original output surfaces](../research/RSDS_SIX_FOLLOWUP_PRIOR_ORIGINAL_SOURCE_EXPOSURE_FIRST_RESULT_20261115.json)
- [Provider file selection, size and MD5/HOLD receipt](../research/RSDS_SIX_FOLLOWUP_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1.json)
- [First new native q and daylength original GPS results](../research/RSDS_SIX_FOLLOWUP_ORIGINAL_GPS_Q_FIRST_COMBINED_RESULT_20261115.json)
- [Strict post-outcome q-position vs seasonal-slope derivative audit](../research/RSDS_SIX_SAME_SPECIES_Q_POSITION_VS_DAYLENGTH_SLOPE_POSTOUTCOME_FIRST_20261115.json)
- [Verifying source-only GitHub Actions run](https://github.com/zuizui0223/RSDS/actions/runs/43712910795)
