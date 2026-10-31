# Divergent diel displacement and astronomical-daylength associations across twelve original animal-tracking studies

**RSDS original-data ecological manuscript v1 — 31 October 2026**
**Stage:** first cross-source synthesis of immutable public GPS results; observational, source-clustered, post-outcome integration.

## Abstract

Astronomical daylength organizes environmental opportunity, but whether animals share a general daily allocation of movement and a common signed seasonal response remains uncertain. Inference can be distorted when multiple taxa from one tracking program are treated as independent ecosystems. We combined two previously source-validated Movebank original-GPS cohorts (eight field studies, twelve study–taxon groups) with four separately original-file-verified multi-native tracking studies (ten additional registered study–taxon groups). All input original-data results were pinned to their previously committed Git blob identifiers, with no telemetry downloaded or newly selected based on estimated effects. The final inventory contains **12 published study DOIs, 22 study–taxon records and 20 distinct registered taxon labels**, of which **nine studies and fifteen study–taxon groups** support recorded night-minus-day straight-line GPS displacement q and within-individual-year q–photoperiod associations.

Nominal original-animal 95% bootstrap intervals for q lay below zero in eleven groups, above zero in four, and were unavailable in seven. For the q–daylength slope, five groups had nominal positive intervals, two negative intervals, eight zero-spanning intervals, and seven were unavailable. The **study-DOI-equal mean** q changed from +0.03464 across seven earlier eligible studies to **+0.01026 across nine**; equal weighting of the fifteen native study–taxon rows instead gave **−0.02251**. The conditional nine-study DOI bootstrap interval [−0.04923,+0.07674] spans zero and leave-one-study-out means change sign. In one new, separately verified field study, **three goose species observed on the same 120 eligible calendar days** all showed relatively greater recorded daytime displacement, yet two had opposing nominally supported slopes against a shared study-reference daylength. A within-month study-level contrast between two geese persisted descriptively, but used only five year–month blocks.

These original-data results show signed heterogeneity in observed diel movement and astronomical-daylength associations at multiple levels, rather than a universal photoperiod response. They do not establish latent activity budgets, solar navigation, climatic causation, or representative worldwide population parameters. Missing-source and earlier radiation-direction external nonreplication outcomes remain intact.

## Introduction

Astronomical daylength varies deterministically with solar date and latitude, while irradiance, temperature, wind, seasonal resource availability and predator encounters follow overlapping but distinct ecological trajectories. Biological use of movement opportunities can respond through the time of activity, the frequency and intensity of displacement, or the direction of travel, none of which is identical to the others. GPS surveys can measure travel chords and their timing, but not true active-state duration without an explicit observation model.

Two inferential errors threaten macroecological generality: comparing night/day displacement records gathered under different collar schedules as though they were identical samples; and treating four taxa tracked in one field study as four independent populations or landscapes. We synthesize originally archived telemetry cohorts at the **published study DOI** level, preserve the underlying native taxon results, and ask whether observed displacement and photoperiod trends converge on common signs. We use a new co-observed goose source to distinguish a taxon's overall diel displacement position from the trend along a **common astronomical calendar**.

This manuscript integrates already observed first outcomes. The addition of four source-selected field studies is not a new blind confirmatory test of any universal effect. Earlier RSDS investigations found GPS acquisition-phase aliasing and month/season support problems; their evidence and original hold decisions remain visible rather than being filtered out.

## Methods

### Original published study frames and evidence integrity

The first eight field-study DOIs include Utah (four native ungulate study–taxon cohorts), Wyoming (two ungulates), and one native taxon per source in wild turkeys, bobcats, red foxes, fishing cats, great egrets and koalas. Four separately source-selected multi-native Movebank publications cover Galápagos tortoises (three registered taxa), migratory ducks (two), African vultures (two), and grazing geese (three). The four new original study DOI/study-ID selections were frozen from source metadata before new study GPS response rows were opened. Original file provider UUID, bytes and MD5 were separately verified. The present synthesis verifies the **immutable first original source-result Git blob SHA** for both strata and the independent four-study validation/metadata snapshot, without re-reading original GPS events or altering prior source eligibility criteria.

The 12 published study DOIs include **22 study–taxon records**, but only **20 distinct registered species labels**, because native taxa recur across earlier field studies. New source schema or tracking-cadence failures are retained as unavailable original-data records. Their absence is not an ecological zero. No observed source was replaced to improve an outcome.

### Endpoints and ecological units

For each tagged animal's accepted original GPS records, the diel displacement contrast was
\[
q=\log(1+\widetilde v_{\mathrm{NIGHT}})
  -\log(1+\widetilde v_{\mathrm{DAY}}),
\]
where \(v\) is the straight-line GPS step displacement speed (km h\(^{-1}\)). A positive q means larger recorded relative nighttime displacement, not a greater proportion of hours actually active. Original source-level q is the equal available-animal mean. The daylength coefficient is a within-original-animal-calendar-year association of q with local astronomical photoperiod (hours), with its original eligibility and bootstrap uncertainty retained.

Previously archived nominal original-tag bootstrap 95% intervals classify signed **native study–taxon** outcomes as positive, negative, zero-spanning, or unavailable. These intervals are neither simultaneous nor multiple-comparison adjusted.

For **DOI-equal** summaries we first average each published study's available native source-taxon point estimates equally, then give each published study one equal weight. A missing source study contributes **no zero**. The alternative naïve study–taxon equal-mean is reported to expose differential weighting of multi-taxon studies; it is not a rival globally identified estimator. To examine source selection uncertainty descriptively we resample the available DOI-level means 5,000 times at a deterministic seed, and calculate leave-one-DOI-out means. This procedure is conditional on a small **purposively selected** source set, not a confidence interval for any world superpopulation.

### Separate within-source co-observed calendar analysis

The new goose original source DOI 10.5441/001/1.fk899541 registered *Branta leucopsis*, *Anser albifrons*, and *Anser anser*. A prior frozen analysis intersected their day/night GPS-qualified **120 common calendar days** over **five year–month blocks**. It calculated source tag q and q–**study-reference** astronomical daylength relationships as a nested study-level observational comparison. These shared-date values must **never replace** the native full-source q/year/daylight values used in the twelve-study synthesis, because the sample calendar and daylight reference differ. The three goose taxa represent a single field-study source. Within-month partial contrasts are treated as exploratory with only five time blocks, not independent cross-population replication.

### Prior scientific nonreplication and telemetry validity

All old source-specific GPS observation-cadence and clock-phase findings remain documented. A prior 48-source negative-direction shortwave radiation association failed its independently held-out directional replication test. Older ERA5 raw JSON response hashes failed a strict 12-context replay while newer canonical physical hourly weather values were reproducible; those historical source identity holds were not reclassified. Neither finding is retested with the present q/photoperiod-only analysis.

## Results

### Source completeness, support and signs

The complete source frame contains **twelve original published field studies**; **nine** have sufficient eligible native source GPS data for q and daylength estimation. Among **fifteen estimable native study–taxon records**, q had nominal original-animal intervals below zero in **eleven** and above zero in **four**. Seven other native source rows were unavailable: one oversized egret source file, three Galápagos tortoise registered-taxon label/schema mismatches, two African vulture study-step cadence support failures and one unsupported northern pintail within a duck source. Source holds and unsupported individual responses did not enter any numerical mean.

Photoperiod associations were more heterogeneous. Of fifteen estimable source–taxon records, **five** original-tag bootstrap intervals were nominally positive, **two** negative and **eight** included zero. The signed trends were not a common response to increasing astronomical photoperiod.

### Source-grain weighting alters the apparent overall movement direction

| Equal-mean method | First 8 study frame | Expanded 12 study frame |
| --- | ---: | ---: |
| Published studies with identifiable q | 7 | 9 |
| Mean of all supported native study–taxon rows | +0.00480 | **−0.02251** |
| Mean of DOI-level equal available native taxon means | **+0.03464** | **+0.01026** |

In the twelve-source frame, the conditional DOI bootstrap 95% interval for q was **[−0.04923,+0.07674]**. The leave-one-source-out DOI-equal mean ranged from **−0.01217 to +0.02656**. The nine study-level q means included four positive and five negative point estimates. Hence neither the apparent positive DOI-equal mean nor the negative naïve native-row mean supports a universal diel movement sign.

The DOI-equal native photoperiod slope was **−0.00638 per daylight hour**, with a conditional descriptive DOI bootstrap interval **[−0.01880,+0.00340]**. As with q, the sample is not representative of world taxa, climatic environments or biome histories.

### Geese: similar relative daily displacement, divergent calendar trends

In the independently published goose source, all three taxa had negative q on the same **120 calendar days**, but with different signed study-reference daylight slopes:

| Native goose taxon | Matched-source calendar q | Shared-study reference daylength coefficient / h |
| --- | ---: | ---: |
| *Branta leucopsis* | **−0.18403** | **−0.02171** |
| *Anser albifrons* | **−0.08391** | −0.01581 (source-tag interval included zero) |
| *Anser anser* | **−0.06918** | **+0.02066** |

The opposite nominal slopes for *B. leucopsis* and *A. anser* within one calendar support frame show that the *position* of a recorded diel displacement niche and its *change along the sampled seasonal clock* are distinct quantities. In an exploratory within-month analysis, the *A. albifrons* minus *A. anser* shared-calendar contrast had a within-month study-reference daylength coefficient of **−0.03594 / h**, with a year–month block bootstrap 95% interval **[−0.11046,−0.00392]**. Only five blocks and three pairwise comparisons were available, and the study-reference daylength is not the actual received photoperiod at every migrating individual location. These results do not identify interspecific competition, photoperiodic physiology or causal solar cues.

## Discussion

Our strongest supported empirical result is that **the signed observational organization of daily displacement and its astronomical seasonal association is non-universal within the available original field-study sample**. The direction of any macroecological q average changes with the sampling unit and with which independent original studies are included, while the magnitude remains source dependent and uncertainty spans zero. This justifies considering the published tracking study and original observation schedule as distinct levels in any hierarchical generalization.

Co-observed species provide an additional axis of evidence: animals using the same source calendar can all record relatively greater daytime displacement yet show opposing q–daylength trajectories. This difference between a niche *location* and an astronomical calendar *slope* is an ecological observational finding; its causal explanation remains open. Equal calendar dates do not imply equal locations, resources, temperature or true solar exposure.

The present evidence does not measure biological activity-state budgets, climatic adaptation, solar compass navigation, wind-assisted orientation or direct radiative stress. We cannot infer a worldwide vertebrate mean, phylogenetic tendency, or longitudinal change from twelve purposively sampled source studies. An independent next-stage test would require original source-level observation-process harmonization, location-appropriate photoperiod/solar geometry, version-verified climate, and externally reserved study DOIs before inspecting their biological effects.

## Source files, code and audit

- [Twelve original published GPS studies / 22 native source–taxon first macro dataset](../research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json)
- [Four immutably hashed original source result contracts](../research/RSDS_12_STUDY_ORIGINAL_GPS_TEMPORAL_MACRO_V1_CONTRACT.json)
- [Source-grain and native first-result audit implementation](../scripts/build_rsds_twelve_original_study_temporal_macro_v1.py)
- [Original data native-taxonomy completeness and source unit regression tests](../tests/test_rsds_twelve_original_study_temporal_macro_v1.py)
- [Detailed Japanese original ecology interpretation](../docs/RSDS_TWELVE_ORIGINAL_STUDY_TEMPORAL_MACROECOLOGY_20261031.md)
- [Original three goose native co-observed calendar and source-HOLD evidence](../research/RSDS_FOUR_MULTINATIVE_ORIGINAL_GPS_ECOLOGICAL_VALIDATION_V1_20261029.json)

**Privacy:** no original tagged animal identifier, coordinate, individual effects, precise event dates/timestamps or telemetry raw bytes are republished.
