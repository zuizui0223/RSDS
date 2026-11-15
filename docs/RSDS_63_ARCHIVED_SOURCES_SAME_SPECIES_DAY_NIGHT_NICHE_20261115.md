# RSDS: does a species retain a day–night GPS movement signature across source studies?

**15 November 2026 | Existing archived source-native results only, no new GPS re-reading.** The first cross-source native atlas was generated under the immutable 2026-10-19 Movebank catalog and exact archived artifact digests, in [original source-native Action #43706508911](https://github.com/zuizui0223/RSDS/actions/runs/43706508911). The analysis source is [the frozen source-native atlas](../research/RSDS_ORIGINAL_SOURCE_NATIVE_SAME_SPECIES_CROSS_STUDY_FIRST_RESULT_20261115.json) and the full [artifact and original DOI contract](../research/RSDS_ORIGINAL_SOURCE_NATIVE_CROSS_STUDY_REPEATABILITY_AUDIT_V1_CONTRACT.json).

## Why this is different from the earlier 51-source sign vote

The original 51 new-source frame contains 35 earlier and 16 Wave3 original studies, while the earlier 12-DOI legacy source frame contains 22 registered native rows. **Together there are 63 distinct original DOI / repository study IDs and 88 registered source-native rows**, but not 63 random, mutually independent ecological populations. The original Wave2 19 source **first** aggregate ZIP artifacts were separately digest-verified and reopened for species-level **aggregate** q, not GPS. This avoids a fundamental mistake: an average of four goose species in one original source is **not** the q estimate of any named one of them.

Of the 88 archived source-native registrations, **38 have a supported q and 33 have a supported q–daylength slope** (the latter is a different endpoint). **17 species names** occur in at least two distinct archived source DOIs; only **6 of those 17** have source-native q supported in two or more. The 11 other repeat names with missing original data cannot be reclassified as a failed biological replication.

## Six identical-species comparison candidates

The original estimator is source-native, animal-equal **NIGHT minus DAY log1p observed GPS straight-line step speed**. All individual-bootstrap 95% intervals below are from the original source dataset. They are *nominal, unadjusted, source-specific* intervals, not meta-analytic confidence intervals.

| Species | Distinct original DOIs with q | Original q sequence, one value per DOI | What original nominal 95% intervals establish |
| --- | ---: | --- | --- |
| *Branta leucopsis* | 3 | −0.16547, −0.17464, −0.02611 | First two negative; third [−0.06660,+0.01469] includes zero |
| *Anser albifrons* | 2 | −0.12099, −0.12730 | Both negative intervals |
| *Vulpes vulpes* | 2 | +0.18975, +0.29043 | Both positive intervals (second source only five original tagged animals) |
| *Odocoileus hemionus* | 2 | −0.02185, −0.00446 | Both negative intervals, Utah versus Wyoming |
| *Ovis canadensis* | 2 | −0.06576, −0.01997 | Both negative intervals, Utah versus Wyoming |
| *Anas platyrhynchos* | 2 | −0.02998, +0.00340 | Old source negative interval; newer near-zero source [−0.00025,+0.00695] includes zero |

**Result accounting:** in these six post-outcome comparison candidates, **five of six species have the same nonzero q point-estimate sign across all supported original DOIs**. Four species have *all their available original source-bootstrap intervals* excluding zero in the same direction: *Anser albifrons*, *Vulpes vulpes*, *Odocoileus hemionus*, and *Ovis canadensis*. *Branta leucopsis* has three negative point signs but only two negative original intervals. *Anas platyrhynchos* has a negative vs positive point sign, but one of the two intervals includes zero: it **does not establish a supported sign reversal**.

Neither 5/6 nor 4/6 estimates a world species fraction: it conditions on repeated species names, available source-specific q, mixed old purposively selected studies, metadata-priority Wave1–3, and multiple after-outcome checks. No uncorrected source CI qualifies as a new confirmatory p-value.

## Biological interpretation: direction, magnitude and seasonal modulation are not the same

The two repeatedly archived ungulates provide a useful internal comparison under separate original field studies:

- *Odocoileus hemionus* has negative full-source q in both Utah (−0.02185, 95% [−0.02321,−0.02047]) and Wyoming (−0.00446, [−0.00578,−0.00309]); yet its q–astronomical-daylength slope is **positive in both** (Utah +0.005619/h [0.004388,0.006921]; Wyoming +0.005538/h [0.004102,0.006973]).
- *Ovis canadensis* also has negative full-source q in both original studies (Utah −0.06576 and Wyoming −0.01997, both nominal CIs negative), but the daylength slopes in these original full-source samples both **include zero**. Thus two species can share an observed negative night-relative displacement level while having different strengths of seasonal photoperiod association.
- *Branta leucopsis* has negative q point estimates in three different source DOIs. The Svalbard migration dataset has a **positive** original q–daylength slope (+0.08000/h [0.01391,0.16333]) while the other two source-specific daylength slope intervals include zero. This is *not* proof of a reversed seasonal mechanism.

Strongly unequal q point magnitudes across source studies cannot automatically be interpreted as evolved intraspecific ecological change. Recording schedules, original tag/year observation support, year and location, migration versus local movement, environmental conditions and source populations differ. In particular, some q bootstrap intervals do not overlap across regions, but those original source tag-bootstrap intervals are **not a calibrated independent-source difference test**. Comparisons of full-source tag-equal q and earlier study-shared-date/q are different estimands: never pool them as the same effect.

## Important source identity and independence checks

The source catalog provides distinct DOI, study ID and original *titles* for the repeated names, including goose agricultural grazing, firework escape and Svalbard spring migration; British wetland red fox versus a two-fox coexistence study; migratory mallards versus landscape-linked mallard movement; and Utah versus Wyoming ungulates. These are different original archive records. However **distinct source DOI and study IDs do not prove distinct animal tags, populations, habitats or observation schedules**. We did not reopen raw original GPS or tag identities to establish disjoint tagged animals, and source titles alone cannot certify biological replication independence. No original GPS coordinates, dates or animal tag identifiers were re-exported.

The original 51-source selection remains strongly selected by source-observability (16/51 q estimable). The older 12 sources were earlier purposive selections. **Do not call this a randomized meta-analysis, species-invariant life-history rule or a world nocturnality estimate.** An additional issue is that domestic/managed source domains occur in the original candidate frame and require source-metadata-only adjudication before any ecological inference.

## Most defensible ecological hypothesis after this audit

> **The relative day/night displacement niche of a species can retain a consistent direction across archived field-study contexts, even when its magnitude and its association with seasonal daylight differ. The degree of consistency is conditional on actual repeated, independent source population support and on a comparable telemetry observation process.**

This is a new **research hypothesis supported by source-specific observational examples**, not yet a validated general law. More informative than a naive all-species night/day average is asking which ecological guilds and contexts maintain same-species direction versus alter it at fixed observation support. A particularly important contrast is that **q position and q–photoperiod association are different biological axes**, as Utah/Wyoming *Odocoileus* illustrate.

## Next data action: test transfer only on predeclared truly independent populations

First compile overlap/continuity documentation at study source level (original provider metadata, no raw person/home trajectories) and define what counts as separate original tagging populations. Then select entirely untouched original external same-species study pairs under fixed same-time sampling windows. Keep the *source-native* rather than across-species source mean. If calendar-aligned original tagged-animal samples are unavailable, report HOLD rather than treating shared species identity as proof of replicated behavior.

The original solar-irradiance → direction held-out negative sign still **NOT REPLICATED**, and the 12 archived ERA5 raw provider source SHA cases still **HOLD**.
