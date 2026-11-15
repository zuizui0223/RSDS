# RSDS: 35 original Movebank sources expose nonrandom measurability and q meta-sign fragility

**2026-11-15. First original-GPS result frozen; source-unit sensitivity separately tested.** No additional GPS was opened after the first successful 19-source pipeline.
[Source SHA-locked 35-study result](../research/RSDS_437_WAVE1_PLUS_WAVE2_35_ORIGINAL_GPS_FIRST_META_RESULT_20261115.json),
[source-unit sensitivity](../research/RSDS_437_35_SOURCE_Q_SIGN_FRAGILITY_FIRST_RESULT_20261115.json).

## Two distinct scientific endpoints, not one ecological global mean

**Endpoint A (methodological):** Can a public source selected *before opening its new original GPS movement outcomes* supply the minimum q evidence for **every registered native species**, and if not where does the measurement support fail?

**Endpoint B (ecological, descriptive only):** Among the source studies passing A, what is the direction/heterogeneity of observed night-minus-day **GPS straight-line step displacement** `q`? What happens when sources, sample waves or unknown effects are varied?

The source population is **437 repository packages, 259 GPS/GNSS packages, 120 metadata-priority packages** (109 one-native, 11 multi-native). It is neither a globally randomized survey of wildlife populations nor a random sample of animal chronotypes. The original SHA ranking selected wave 1 (8 single + 8 multi) before its original GPS outcomes. Wave-2 *quota* was chosen after looking at wave 1, but the *19 new identities* were still determined by ranks already frozen in wave 1: 16 further single-native sources and all 3 remaining multi-native sources. **Total 35 source study-ID units, with 50 original registered native study×taxon rows.**

## What ran on original GPS (not hypothetical or synthetic)

[Wave-2 original first execution](https://github.com/zuizui0223/RSDS/actions/runs/43692435152) validated original source file MD5, exact source/native identities, 1h UTC source clock eligibility and original 0.5–2h segment q, and stored all first source outcomes. **19/19 source-specific CI jobs and the combined first 35-source meta job succeeded.** A CI job succeeding includes a scientifically genuine HOLD or absence of supported q.

The first combined result was copied byte-for-byte with SHA256
`96a34f572ad8109e2dab5ac28202faceb6b7750724d33b96b9ab94dbdd16f48d`;
no raw animal GPS/time/position or tag data were committed.
The [independent source-unit sensitivity CI](https://github.com/zuizui0223/RSDS/actions/runs/43693286399)
and its first artifact are also source-sha reproducible.

### Gate-by-gate evidence (not success after deleting failures)

| Gate | Source studies | Interpretation |
|---|---:|---|
| Source SHA selected | **35** | One original source ID, not 35 independent ecosystems |
| Original file eligible under MD5/200MB and not human privacy | **28** | Seven original identities retained as preaccess HOLD |
| Original GPS parser/schema ran without source-wide HOLD | **25** | Three additional source-specific HOLD |
| At least one native species has supported q | **13** | 12 parsed sources lacked a supported native q |
| All registered native taxa supply supported q | **10** | Three partial multi-native sources not admitted as complete |
| All registered native taxa supply q~daylength | **8** | Separate association endpoint, not causality |
| New multi-native same-calendar/month covariance | **0/11** | Failing data support is not evidence of absent biological partition |

Of 35 selected, 7 were original-file/privacy HOLD, 3 source-wide original GPS HOLD, 12 parsed-but-no-native-q, 3 partial-native-q, and 10 source-complete-q. The human-tracking source remained private and unaccessed; its place in the **metadata** denominator is retained but it is **not wildlife evidence**.

### Source q effects (complete groups only)

Wave 1 contributed 4 complete studies; wave 2 added 6 single-native studies. Wave-2 source q (not latent activity proportions):

| Native study source | Observed source q |
|---|---:|
| *Sus scrofa* | +0.06594 |
| *Pelecanus occidentalis* | +0.48955 |
| *Puma concolor* | +0.09236 |
| *Anas platyrhynchos* | +0.00340 |
| *Felis catus* (pet-source) | −0.00293 |
| *Caracal caracal* | +0.05776 |

The 10 source-complete q signs are **7 positive, 3 negative**. This does not establish that the world is nocturnal. The original GPS records have different movement ecology, tracking schedules, geography, study design, season and management contexts. *Felis catus* and domesticated animals are not automatically wild populations.

## The interesting positive methodological result: sign instability

The original source metadata strata are **109 single-native and 11 multi-native**. Wave 1 sampled 8 of each, and wave 2 completed the full 11-source multi-native metadata frame but only 24/109 single-native sources. The respondent-only poststratified source mean—weighting observed source-complete stratum means by the original metadata frame composition **109:11**, NOT by number of animal GPS records—changed from:

**Wave 1 only:** q = **−0.0477935** (4 complete respondents)

**Waves 1 + 2:** q = **+0.0756846** (10 complete respondents)

This is a source-sampling *composition* sign reversal, **not** evidence that animal activity changed with time. The Wave-2 sample was not a wholly blinded replication of an ecological directional effect; the between-wave quotas were decided after Wave-1 outcomes. Pseudorandom source SHA rank is not a truly randomized global population sampling design.

Recalculating across the frozen 10 completed sources:

- Leave **one original source study** out: q remains positive, between **+0.02155 and +0.10705**.
- Leave **two original source studies** out, while retaining at least one complete source in both metadata strata: q ranges **−0.00440 to +0.13300**. The negative end removes *Pelecanus occidentalis* and the earlier two-fox study (*Vulpes vulpes* + *V. lagopus*) together. This was selected AFTER examining the first source effects, and is only a diagnostic.
- Of the **35** fixed original source studies, **25** still lack a full registered-native q. **Only under the hypothetical additional assumption** that the missing selected-source q values lie inside the range of the ten observed full-source q values does the respondent-source frame mean lie between **−0.09501 and +0.35688**. This is neither a valid confidence interval nor a partial-identification bound without that assumption.

The actual inference problem is not rescued by applying a conventional random-effects estimator to the ten survivors. A source's original GPS sampling architecture influences whether it survives, and thus the distribution of q among survivors.

## Finite metadata-frame source measurability, not ecological mean

In the now fully examined **11/11** originally selected multi-native metadata frame, source-complete q was available in **2/11 (18.2%)**. Among the **24/109** single-native sources selected through two SHA waves, **8/24 (33.3%)** passed.

A purely descriptive source-measurability projection using the 109:11 metadata composition is **31.94%**; it assumes the original SHA ordering resembles an outcome-independent selection process, and cannot be treated as a formally randomized design because Wave-2 quotas were adapted. The strict *logical* complete-source count among the finite 120 metadata candidates lies between the **10 observed** and **95 maximum** (all 85 uninspected single-native candidates potentially pass), or 8.3–79.2%. This is deliberately too broad to support a global ecological inference.

A separate domain problem remains: one source is a protected human-location source and some animal sources are managed/domestic. A wildlife-only target requires independent metadata-based source-domain adjudication, not retrospective exclusion based on q.

## What to improve next (do not overwrite previous first results)

1. **Data-architecture harmonization:** classify cadence/sensor duty cycle from source-only technical metadata and source-event sampling support. Preserve the exact original 1h/0.5–2h result as a reference. Any more permissive observation protocol requires a separately frozen new estimand; it may recover source q but cannot be reclassified as earlier independent confirmation.
2. **More single-native original study slots:** the remaining 85 metadata-eligible SHA-rank sources remain untouched. If expanded, fix quotas/stop rules before event access, preserve native/nonresponse denominators and source identity overlaps.
3. **A valid multi-level meta-analysis** needs comparable full-source q definitions, full joint within-source covariance where possible and enough genuinely independent source DOI/study IDs to estimate heterogeneity without pseudo-replicating taxa or tags. The zero source-month covariance in this wave is a data-support failure, not biological absence.
4. **Missingness sensitivity** should be treated as an estimand; source effects among complete respondents alone are selected conditional on measurement support and cannot be extrapolated as global night/day movement.

Old **48-source held-out directional negative shortwave effect remains NOT REPLICATED**, previous **12 original ERA5 raw-SHA versions remain HOLD**, no causal photoperiod, solar compass or activity-clock claim is made.

**Potential paper framing:**
*When movement datasets become meta-analytically comparable: observation-dependent filtering and apparent temporal niche generality in public animal GPS archives.*

Methods novelty is the **fully auditable observation-stage source admission + registered-native joint measurement + downstream transfer/sign sensitivity**, not SHA ranking, random-effects meta-analysis, or leave-one-study-out itself.
