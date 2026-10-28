# Astronomical daylength, diel displacement, and solar-referenced travel across eight original animal-tracking studies

**RSDS original-data ecology manuscript — October 26, 2026 — exploratory v1**

## Abstract

Astronomical daylength and solar azimuth structure environmental opportunity, but animal activity timing and actual travel direction are different aspects of realized movement niches. Whether these axes show recurrent ecological patterns across species and regions remains unresolved. We extended an original-data study of Utah and Wyoming ungulates (two Movebank source studies, six source–taxon outcomes) with six newly fixed original Movebank surveys selected on source metadata before opening their animal GPS events. The new studies span wild turkeys, bobcats, red foxes, fishing cats, great egrets and koalas. Five new tracking sources passed provider checksum and fixed file-size qualification; the egret original GPS exceeded the predeclared 200 MB cap and was retained as a source metadata hold. The combined archive represents **eight original study DOIs and twelve study–taxon records**; eleven records support an original animal-level night-minus-day observed GPS minimum displacement contrast, ten support solar-relative heading, and nine support repeated-original-tagged-animal comparisons between two coverage-selected observational contexts (288 tag-pair occurrences, not independent ecosystems).

The observed contrast \(q=\log(1+\widetilde v_{\rm night})-\log(1+\widetilde v_{\rm day})\) had source-tagged-animal bootstrap 95% intervals below zero in seven records and above zero in four, contrasting day-biased and night-biased **recorded displacement**, not activity duration. Its within-individual-year association with astronomical daylength had nominal intervals above zero in five source–taxon records, below zero in one and crossing zero in five. Most strikingly, newly observed wild turkeys had an original-study daylength coefficient of **−0.02630** (source-animal descriptive interval [−0.03719, −0.01533]), while new bobcats had **+0.01282** ([+0.00155, +0.02320]), and previously observed mule deer had positive within-source coefficients in two distinct field studies. Source-mean cosine of GPS travel bearing minus Sun azimuth showed two positive, two negative and six unresolved nominal bootstrap interval signs among ten estimable records. Joint circular covariance was available with source-animal bootstrap uncertainty for only four new sources, all intervals crossing zero; the older six estimates were point-only and were not silently granted confidence intervals.

The results establish *observational nonuniversality* of signed diel displacement and solar-relative movement in this purposively selected sample. They do not yet identify daylength memory, sunlight navigation, weather mediation, causal environmental plasticity, or a globally representative phylogenetic parameter. Historical failure to replicate an independent 48-source radiation-heading hypothesis and twelve original ERA5 response-version holds are preserved.

## Introduction

Solar elevation, photoperiod, solar azimuth and apparent solar time are astronomical properties; local air temperature, wind and shortwave irradiance covary but differ in their spatial and seasonal trajectories. Organisms encounter these dimensions through species- and habitat-dependent opportunities for movement, predation avoidance, foraging, soaring, rest and migration. Aggregate time-of-day records and Sun-relative GPS travel bearings therefore need not obey the same ecological generalization rule.

We evaluate three separate observational propositions:
1. **Temporal niche diversity:** sampled native source taxa do not all exhibit the same relative night-versus-day GPS displacement.
2. **Astronomical seasonal heterogeneity:** daylight duration may be associated with relative recorded displacement, but with heterogeneous slopes across field-study source contexts rather than one universal sign.
3. **Direction decomposition and limited transfer:** a source-average travel-to-Sun cosine has no clear universal sign and its geographical-heading marginal must be distinguished from the residual joint circular covariance of geographical bearing and Sun azimuth.

These are not predictions that a preselected new sample must share the sign of any previous source; all eligible, null, unsupported and source-HOLD results remain included.

## Data, pre-event source freeze and field-study units

The input universe was the existing fixed catalog of 437 original Movebank packages. The first two published tracking studies were Utah ungulates (DOI 10.5441/001/1.711; four native taxa) and Wyoming mule deer/bighorn sheep (DOI 10.5441/001/1.703; two taxa). Independently of the new animals' movement results, six additional source UUIDs/DOIs and study IDs were fixed on October 26, 2026 to broaden the source sample across movement guilds and geography, after excluding previously used original study IDs and DOIs. This selection is **purposeful and post-prior-result**, not a probability sample or a blind preregistration of sign-direction hypotheses.

Five original new source tracking files were identified unambiguously and verified against source-publisher bitstream UUID, bytes, and MD5 checksums. The other (great egret) tracking file exceeded the prespecified 200 MB size cap; it remained a full denominator record with no behavioral effect, rather than being replaced. Accelerometer-only streams were excluded by original file NAME alone before the new GPS outcome was accessed. All animal tracks, names, exact dates/times, locations and per-tagged-animal coefficients stayed ephemeral and were not republished. Each DOI represents one original study, regardless of the number of tags, fixes, calendar years or co-observed taxa.

## Observational endpoints and uncertainty

We reused original timestamp-only one-fix-per-UTC-hour sampling and the previously specified 0.5–2-hour GPS displacement interval constraints, 60 km/h maximum reported step displacement speed, and solar-day/night consistency checks. Source animal days with at least three admitted daytime and three nighttime steps, and original animal-years with at least ten such days, permit per-tag \(q_{it}\) estimates. Each native original study's q is a group mean of per-animal eligible night/day \(\log(1+\mathrm{km\,h^{-1}})\) displacement differences and is not an activity-hours or time-budget outcome.

The temporal-daylength endpoint is the within-animal-calendar-year OLS association of daily q with astronomical daylight hours, requiring sufficient original days and daylight range, then averaging slopes equally within original source animals. Source-tagged-animal nonparametric bootstrap 95% percentiles are descriptive, and source studies remain the cross-population units.

For a moving daylight GPS segment with bearing \(B\) and contemporaneous solar azimuth \(A\), Sun-relative direction was measured by \(Y=\cos(B-A)\). Distinguishing the geographical bearing distribution is essential:
\[
E\{\cos(B-A)\}=
E(\cos B)E(\cos A)+E(\sin B)E(\sin A)
+\mathrm{Cov}(\cos B,\cos A)+\mathrm{Cov}(\sin B,\sin A).
\]
The first two terms constitute a geographic/Sun marginal composition and the last two joint circular covariance. The original Utah/Wyoming cross-axis matrix provides point estimates of joint covariance **without comparable per-source-animal confidence intervals**; the four supported new sources provide bootstrap intervals. We do not fabricate old intervals, count them as significant or declare their true coefficients zero.

We separately evaluated whether the exact same tagged individuals were observed in two 1° grid×calendar-year contexts chosen for maximal q-qualified source animal overlap (not weather, photoperiod branch, or direction sign) and could contribute both q and daylight solar-heading measurements. Pairing conditions were met in nine of twelve source–taxon rows, with 269 occurrences in the earlier two studies and 19 across wild turkey, bobcat and koala, totaling 288 source-pair occurrences; overlap does **not** create 288 independent ecological populations.

## Results

### Native source outcomes: time allocation proxy

Eleven source–taxon rows supported q. **Seven** had nominal source-tagged-animal bootstrap intervals entirely below zero (including six older ungulate rows and the new wild turkey survey), whereas **four** had intervals entirely above zero (bobcat, red fox, fishing cat and koala). This is not proof of a guild-specific day/night adaptation, because taxon, region, season, telemetry and ecology are confounded.

The within-animal-year q-versus-daylength association was positive with an interval excluding zero in five rows (four older ungulate source–taxon rows and new bobcats), negative in one (wild turkeys), and unresolved in five other estimated rows. The egret remained a metadata HOLD. The opposing signed daylength associations in newly sampled wild turkeys and bobcats falsify a simplistic *source-sample universal common signed response*; they do not identify a controlled species or life-history explanation.

### Solar-referenced geographical movement

Ten source–taxon rows supported the observed raw solar-relative cosine. Two had positive, two negative and six zero-spanning individual-bootstrap nominal intervals; fishing-cat heading had only three qualifying tags and egret GPS was not available under the file-size cap. The four newly supported joint circular covariance interval estimates all overlapped zero. Six older joint-component values were point estimates without source-animal intervals, preventing pooled sign-confidence claims. Thus there is no demonstrated global directional law and certainly no solar-compass or light-avoidance effect established by this source-level decomposition.

### Repeated original tagged individuals

Across all eight field studies, nine source–taxon records had two coverage-selected contexts with repeated original tagged individuals. Three of these were from new source studies: wild turkeys (eight paired animals), bobcats (six) and koalas (five). In the remaining new sources, the identical original-tag pairing was insufficient or unavailable; the lack of pairing cannot be interpreted as no biological temporal change. The corresponding two contexts per source were not sorted by solar irradiation, temperature or common astronomical season. No common exposure response can be identified from their signed contrast alone.

### Post-result observation-process audit: native UTC GPS collection phase

A subsequent source-checksummed telemetry quality analysis revisited exactly the same eight original field-study DOIs and twelve native study–taxon cohorts, retaining the original Great Egret file-size HOLD. The original baseline native day/night q and q–daylength estimates were reproduced before alternate temporal subsampling. It became apparent that a *one-fix-per-UTC-hour-bin source-selection rule* does **not** imply a true hourly collar acquisition schedule. Among eleven original GPS-estimable study–taxon cohorts, all four Utah ungulate groups had only 4.38–12.92% of their chosen source fixes in EVEN UTC-hour bins, while the other seven source groups were substantially closer to 50%.

To distinguish a genuine collection-clock artifact from changes in which animals and dates happened to qualify, we compared the original source record, EVEN-only 1.5–2.5-hour moving segments and ODD-only 1.5–2.5-hour segments using **exactly the same tagged animals, calendar years and astronomy-qualified observation dates eligible in all three schedules**. Eight native source–taxon groups had at least five triple-matched original tagged animals, three were not estimable on this stricter denominator and one remained the original Great Egret file-size HOLD.

In the independent Wyoming ungulate tracking source (DOI 10.5441/001/1.703), the observed q EVEN−ODD measurement difference was **+0.005393** (113 identical mule-deer tags, nominal animal-bootstrap 95% [+0.003574,+0.007132]) and **+0.013339** (41 identical bighorn tags, [+0.011617,+0.014897]). Relative to the same ORIGINAL record and same days, both species had a nominal **positive EVEN−ORIGINAL** q difference and a nominal **negative ODD−ORIGINAL** difference. Therefore, the phase-dependent measurement differences cannot be attributed solely to differing animals or recorded dates across two separate comparisons. They may still arise from the phase at which the same animal’s movement chord is sampled, travel tortuosity, missing fixes, collar schedules and changes in the segments eligible under the day/night rule.

Across all eight original field studies, only **three** study–taxon cohorts had an EVEN−ODD difference whose *unadjusted descriptive original-tag* bootstrap interval excluded zero on the positive side: both Wyoming taxa and Utah elk (the latter only six original triple-matched tags in a strongly ODD-acquisition-biased study). The remaining five estimable cohorts had zero-spanning intervals. **These are observation-process outcomes, not direct phenological, sunlight or behavioral adaptation tests.** Study sources, not study–taxon rows or individual GPS fixes, determine the level of ecological replication. Consequently, signed q or q–daylength contrasts from the eight original studies must remain explicitly conditional on the observation schedule until the source-specific cadence and measured step length are modeled.

The observation audit does not refit any ERA5 effects. The original radiation–direction held-out negative-direction replication failure remains unchanged; the original ERA5 raw provider JSON version holds remain historical, despite later successful numeric reproduction of the fitted coefficients on an independently versioned current ERA5 physical-data snapshot.

**Original-data source receipts:** [eight study original 1h-bin vs EVEN/ODD q audit](../research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json), [native GPS clock-phase alias](../research/RSDS_EIGHT_STUDY_GPS_CADENCE_PHASE_ALIAS_DIAGNOSTIC_20261028.json), [strict same-tag/year/date triple-source first result](../research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json), and [SHA-pinned biological interpretation and support denominators](../research/RSDS_GPS_TRIPLE_DATE_CADENCE_ECOLOGICAL_INTERPRETATION_V2_20261028.json).

## Ecological implications and inference limits

The strongest supported finding is that the **observational organization of daytime/nighttime displacement is heterogeneous and does not follow an invariant signed daylight association** across the broader set of original field studies. The comparison extends the prior Utah/Wyoming tests to additional birds, carnivores and a marsupial, but it is not an eight-ecosystem random-sample estimate of a phylogenetic process. Observed night/day GPS speed can differ from true biological activity duration because of duty cycles, behavioral state, habitat, step length and device accuracy. Source study and taxon effects cannot yet be separated in the new single-taxon field archives.

Apparent Sun-directed movement depends on the geographical bearing marginal and sampling distribution of Sun azimuth. Raw orientation alone is an especially poor mechanistic proxy. The all-positive/all-negative hypothesis for joint covariance is not confirmed; insufficient old source CIs and new-source uncertainty preclude a global signed estimate. Demonstrating an actual photic mechanism requires a new, independent set of tagged-original movement intervals with fully version-verified hourly radiation, direct geometric solar altitude/azimuth, local time, wind and heat controls, and prespecified original-study/animal random effects. Existing ERA5 source SHA mismatches should be investigated by canonical numerical-physical arrays rather than silently interpreting provider serialization changes as animal biology.

Only source-derived aggregate data, fixed source SHA checks and public DOI-level attribution are released. Earlier externally nonreplicated radiation direction and all original GPS/weather supports or HOLDS are unchanged.

## Reproduction and related RSDS study artifacts

- [Eight original study macro dataset](../research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json)
- [Pinned source-equal macro synthesis implementation](../scripts/build_rsds_eight_study_temporal_solar_macro_atlas_v1.py)
- [Source manifest originally sealed before new GPS event access](../research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_V0_20261026.json)
- [Five new source observations and one fixed HOLD](../research/RSDS_NEW_SIX_ORIGINAL_GPS_ECOLOGY_FIRST_RESULT_20261026.json)
- [Previous Utah/Wyoming first joint time and Sun-direction cohort results](../research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json)
- [Japanese ecological synthesis and limitations](../docs/RSDS_EIGHT_STUDY_TEMPORAL_DIRECTION_MACROECOLOGY_20261026.md)
