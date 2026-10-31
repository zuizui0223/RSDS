# Daily movement niche rank can recur across sites while astronomical seasonal trajectories diverge

**RSDS original-GPS comparative ecology manuscript v1 | 31 October 2026**  
**Status:** exploratory post-result synthesis of four originally source-verified field studies, not a prospective independent validation of causal photoperiod effects.

## Abstract

Animals may differ in their relative daytime versus nighttime movement and in how that balance changes along the annual astronomical clock. These are distinct dimensions of temporal niche organization, but most comparisons confound species with recording calendar, study region and telemetry schedule. We assembled previously source-verified original Movebank GPS results from four independently published field studies with fixed, within-study co-observed qualified calendar dates: Illinois coyotes, bobcats and white-tailed deer (169 days); Wyoming mule deer and bighorn sheep (193 days); Utah mule deer and bighorn sheep (199 days); and three goose taxa (120 days). The dataset represents **four source-study DOIs, ten study–native-taxon groups and eight unique taxon labels**, not ten independent ecosystems.

For each study, we distinguished the position of the observed night-minus-day GPS straight-line displacement contrast \(q\) from its association with a common study-reference astronomical daylength. All four study sources had at least one native pair with a nominal original-tag bootstrap interval excluding zero for their relative q position. In the *same repeated native pair*, mule deer had a higher q than bighorn sheep in both Wyoming (+0.01702, original tagged 95% interval [+0.01416,+0.01987]) and Utah (+0.03767, [+0.02158,+0.05716]), across the original and both 2-hour UTC fix-phase sampling schedules. Yet Wyoming showed positive mule-deer and negative bighorn daylength slopes, while Utah's native slopes were both positive and their between-species slope contrast was inconclusive. Oppositely signed nominal native daylength slopes were also seen among Illinois and goose taxa on their respective shared source calendars. Thus **three of four** selected sources contained opposite-signed photoperiod-associated native trends, even when native ranks in diel GPS displacement were resolved. Exploratory within-month contrasts retained one signed difference in each of Illinois, Wyoming and the goose source, but used only five to nine year–month blocks.

These patterns suggest that the *relative position* of a recorded temporal displacement niche need not predict its *astronomical seasonal trajectory*. They do not establish physiological daylength sensing, sympatry, interspecific competition, direct solar radiation responses, real activity duration, or worldwide generality.

## Introduction

The same animal can redistribute observed GPS movement between astronomical day and night as seasonal conditions change. A source group's mean night-minus-day displacement describes its recorded daily movement position, while a daylength coefficient describes change along the sampled annual calendar. Neither is a direct measurement of time actually active.

For multi-native studies, a stronger comparison holds the source observational calendar fixed: registered taxa appear in the original study on exactly the same qualified solar dates, reducing some date-support differences without equating individual environments or detection schedules. A particularly useful diagnostic is whether the ranking of the same two native taxa recurs in another source region while their daylength-associated trajectories also recur.

This study synthesizes such immutable original source outcomes, separating within-source native comparisons from between-source evidence. It is explicitly exploratory: the Utah mule-deer–bighorn pair was selected after inspecting Wyoming, and the four original sources were not a random sample of the Earth's ecological communities.

## Original GPS evidence and methods

We source-locked five published source-result receipts to their exact Git blob SHA: Illinois three-native first source, Utah/Wyoming two-source same-native calendar, new three-goose validated calendar, Illinois/Wyoming month-conditional source results and a separate 12-study full-source temporal macro summary. Each original source had previously passed its original Movebank source-file checksum and signed native q/daylength support conditions. We did not access animal GPS positions, track IDs, precise dates, weather arrays or restricted raw files.

The original response was
\[
q = \log(1+\widetilde v_{\rm NIGHT})
       -\log(1+\widetilde v_{\rm DAY}),
\]
where \(v\) is straight-line GPS displacement speed in km/h from eligible consecutive fix segments, not true time spent in activity states. Signed native study-reference photoperiod slopes were measured on qualified within-source shared dates. Original source-tag bootstrap percentile intervals of 95% were used to classify native source means as positive, negative or unresolved. These are nominal, multiple-comparison-unadjusted intervals, not simultaneous population inference.

Four study calendars were independently preserved: Illinois 169 days in nine year–month blocks, Wyoming 193 days in eight, Utah 199 days in ten, and geese 120 days in five. The astronomical daylength reference was 43°N for Illinois and the two ungulate studies, but an independent study-specific reference latitude for geese; magnitudes cannot be interpreted as matched direct daylight exposures between source geographies.

For the repeated mule-deer–bighorn pair, original q rank and daylength trends were checked in all three archived source GPS schedules: original UTC hourly-bin selected source fixes, even-hour-only 2-hour intervals, and odd-hour-only 2-hour intervals. This is an **instrument sampling** sensitivity, not an experimental light treatment. Year–month-centered analyses of q species contrasts from Illinois, Wyoming and the goose source were included only as prior source-specific exploratory checks. Comparable within-month analysis was **not** performed for Utah.

## Results

### Study-scale summary

| Original published field study | Qualified common source dates | Native taxa | Native q pairs with a resolved nominal q rank | Two opposing signed native daylight slopes? |
|---|---:|---:|---:|---|
| Illinois | 169 | 3 | 3/3 | Yes |
| Wyoming | 193 | 2 | 1/1 | Yes |
| Utah | 199 | 2 | 1/1 | No |
| Goose source | 120 | 3 | 2/3 | Yes |

This totals **four field studies and ten native source rows** (eight unique taxon labels). The sum of 681 source-calendar dates is not 681 independent ecological units, and overlap in source calendar date does not establish spatial co-occurrence.

### Recurrence in the mule deer–bighorn relative movement position

Mule deer had nominally larger q than bighorn sheep on shared original source calendars in both Wyoming (**+0.017018**, original tagged interval **[+0.014157,+0.019873]**) and Utah (**+0.037671**, **[+0.021581,+0.057158]**). The positive rank was retained in both even/odd GPS sampling-phase subsets as previously archived. This is source-conditioned repeatability of a **relative observed displacement measure**, not evidence of identical time budgets.

Daylength associations contrasted with the q rank. In Wyoming, source-native mule-deer and bighorn daylength slopes were **+0.007827/h** and **−0.004520/h**; the nominal between-native contrast was **+0.012347/h** with a positive interval. In Utah the corresponding source-native slopes were **+0.013736/h** and **+0.040937/h**, both positive, and the native contrast **−0.027201/h** had an interval spanning zero. Thus a repeated q ranking did **not** entail a repeated supported species difference in seasonal astronomical response.

### Other source calendars show heterogeneous temporal slopes

Illinois three-native original tracking records on 169 dates showed q: coyote **+0.26852**, bobcat **+0.09339**, white-tailed deer **+0.00245**, with all three native-pair q ranks supported by conditional original-tag intervals. Their source-reference daylength slopes were **+0.03744/h** (coyote), **−0.02343/h** (bobcat) and unresolved near zero (deer). The goose source on 120 common dates had negative q in all three natives (*Branta leucopsis* **−0.18403**, *Anser albifrons* **−0.08391**, *Anser anser* **−0.06918**), yet the first and third taxa had supported opposing signed source-reference daylight trends (**−0.02171/h** vs **+0.02066/h**). Native source ranking and daylight slope were therefore demonstrably different observational outcomes in both mammal and avian field studies.

In previous year–month-adjusted analyses, Illinois bobcat minus white-tailed deer had a nominal negative within-month calendar daylight contrast **−0.08466/h**; Wyoming mule deer minus bighorn had a nominal positive within-month original-fix contrast **+0.01027/h**; the goose source's *A. albifrons* minus *A. anser* had a nominal negative contrast **−0.03594/h**. These are not independent replicated effects and have only 5–9 month-level blocks per study. Utah did not undergo the same within-month check.

## Discussion and inference ceiling

Two conclusions are supported at the observational source level. First, common-calendar native taxa show distinct daily movement positions, and the same mule-deer–bighorn q rank can recur in two distinct tracking studies under different GPS timestamp sampling phases. Second, those native positions do not determine a universal direction of change along the annual astronomical clock: supported opposite native slopes appear in three of the four selected field studies, and the same ungulate pair's seasonal contrast does not transfer automatically from Wyoming to Utah.

The appropriate unit of between-study replication is the **field-study DOI**, not native pair comparisons, tags, fix records or a sum of common calendar dates. Each source contains correlated sampling and biological histories, and only two distinct studies observe the repeated mule-deer–bighorn pair. Native study-reference daylength is deterministic with date; it can covary with weather, resource availability, migration, breeding condition, observer/GPS duty cycle and spatial habitat. Therefore no causal daylength sensing, community competition, sun compass or thermal mediation has been identified.

This four-source co-calendar frame is **not the same data population** as the separate 12-study full-source q/daylength macro. Illinois adds one co-calendar source DOI outside that twelve-study fixed macro universe, and within-source shared calendar q does not use the same dates/tags as native full-source annual q. Treating both as one homogeneously weighted 13-study numerical mean would be invalid. The earlier externally nonreplicated 48-source solar-radiation directional result and the original 12 ERA5 raw provider SHA version holds remain unchanged.

A prospective ecological generalization should preselect original multi-native field studies independent of q or daylight slope, require comparable GPS observation effort and same-site/day resource context, and then estimate source-clustered contrasts with explicit solar geometry, temperature, wind and geography.

## Reproducibility and source evidence

- [Four independently published same-calendar original source study outcomes](../research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_FIRST_RESULT_20261031.json)
- [Frozen original source DOI, source-calendar/latitude and SHA contract](../research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_CONTRACT.json)
- [No-GPS-retrieval original source comparison script](../scripts/build_rsds_four_coobserved_calendar_ecological_niches_v1.py)
- [Original study native rank and seasonal trends tests](../tests/test_rsds_four_coobserved_calendar_niches_v1.py)
- [Detailed Japanese ecological synthesis](../docs/RSDS_FOUR_COOOBSERVED_SOURCE_STUDY_DIEL_POSITION_VS_DAYLIGHT_SLOPE_20261031.md)
- [Separate 12-study full-source temporal macro](../research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json)

No raw original GPS records, tagged animal identifiers, precise fix dates/times, coordinates, private source grid locations or animal-level coefficients were published.
