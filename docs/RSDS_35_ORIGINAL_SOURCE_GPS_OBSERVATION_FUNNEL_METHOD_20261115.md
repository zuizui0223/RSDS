# RSDS original GPS observation funnel — separating telemetry measurability from temporal ecology

**15 November 2026. Post-outcome method diagnostic on the fixed Wave1+Wave2 source frame.**
Historical first source GPS event results and their scientific support classifications
are immutable. The new route does **not** try to rescue unsupported time-niche q
by relaxing observation thresholds, substitute convenient taxa, or publish raw tracks.

## Main question

**Which parts of the original animal-GPS observation process prevent a
public dataset from contributing a comparable night–day movement effect?**

This is necessary before interpreting heterogeneity across studies as ecological
niche variation. A positive/negative source effect and the probability that a
source can *produce* that effect are separate quantities.

The [original 35-source source/biological status](../research/RSDS_437_WAVE1_PLUS_WAVE2_35_ORIGINAL_GPS_FIRST_META_RESULT_20261115.json)
is treated as frozen. It contains 35 registered study-ID/DOI source units
and 50 registered source-native taxon rows, including the following
original stages:

| Fixed original stage | Sources |
|---|---:|
| Selected source package and registered native labels | 35 |
| Original provider MD5/bytes admitted (human source privacy restriction and file budget intact) | 28 |
| Original GPS parsing and native denominator available | 25 |
| At least one original registered native q estimable | 13 |
| Every originally registered native q estimable | 10 |
| Every originally registered native q~daylength estimable | 8 |

Of the original 35 sources, seven were provider-file/human privacy HOLD;
three had source-wide GPS/schema HOLD; 12 parsed but lacked any native q;
three had only partial native q; ten had complete native q. These are
**distinct** measurement gates. A previously failed source remains failed
even if its repository file has since changed.

## New source-native *interval* gates (original thresholds unchanged)

A private original GPS source is read under the same registered bitstream
UUID, original file size, provider MD5 and original source-native CSV schema.
The original parser retains **at most one original fix per UTC-hour bin**
per original tagged animal. No interpolated new observations.

For **each original registered native taxon** in a parsed original source,
audit the following nested stages, before averaging native effects:

```text
original UTC-hour-selected GPS fixes
      |
adjacent original tagged-individual fix pairs
      |
0.5 <= original elapsed hours <= 2.0                 time-only opportunity
      |
same astronomical solar date and DAY/NIGHT state    local solar continuity
at pair start, midpoint and endpoint
      |
finite great-circle GPS speed <= 60 km/h             geometric plausibility
      |
>=3 qualified DAY + >=3 NIGHT steps                  original tag x solar date
      |
>=10 such dates in an original tag x calendar year   native q-tag opportunity
      |
>=5 independently tagged animals per native source   original native q identifiable
      |
ALL registered source-native taxa identifiable       source-complete q
```

Separately, a native daylength association requires adequate original
within-animal-year astronomical daylength span (>=0.5 hours) and
>=5 eligible original tagged individuals. **Each node has a different
denominator.** A large timestamp dataset cannot substitute for missing
hourly intervals, matching DAY/NIGHT source observations, valid geometry
or original tag-year support.

The technical audit records, by original source-native taxon: the number
of original hourly fixes, adjacent pairs, time-eligible pairs,
solar-period-eligible pairs, geometry-eligible pairs, qualified
tag-date records and qualifying tag-years/animals. It also reports the
fraction of original hour-selected fixes at even vs odd UTC hours.
This parity statistic is a source sampling-schedule diagnostic,
**not** an inferred biological circadian rhythm.

At each native source, compare the new nested counts with the **unchanged**
`one_taxon` original total adjacent steps, qualified displacement
steps, date counts, q-tag counts, daylight-tag counts and native q
mean. If any differ, the route fails closed; do not silently
reinterpret the original endpoint.

## Exact accounting identity and the non-collapsible animal denominator

For any source-native taxon with positive interval counts, define
`N0` adjacent original selected-fix pairs, `N1` elapsed-time-eligible
pairs, `N2` pairs additionally continuous in solar date/daypart, and
`N3` pairs additionally passing the fixed speed gate. They are
**nested integer counts**:

```text
N3 / N0 = (N1 / N0) × (N2 / N1) × (N3 / N2)
```

The right side is an exact **decomposition of observation admission**,
not a fitted model, probability estimate, or causal path analysis.
The denominators are original GPS intervals **within** an animal,
and thus cannot serve as independent ecological replication units.

Crucially, the subsequent source-native group indicator

```text
R_source,native(q) = 1{ >= 5 distinct tagged animals each supplying
                          >= 10 paired DAY/NIGHT-qualified solar days
                          within at least one original calendar year }
```

**cannot** be determined by multiplying source-wide interval
acceptance fractions. Two original sources with similar `N3/N0`
can have radically different counts of eligible *independent tags*.
This is the key distinction the diagnostic quantifies and verifies.
No independence, equal duty cycle, or Missing At Random assumption
is created by large `N3`.

## What the diagnostic can and cannot identify

- The first *time-interval gate* isolates an observation-cadence constraint
  on the selected original timestamp series.
- The solar-date/daypart gate reflects astronomical state classification
  and the animal's latitude, longitude and recording time; it is
  **not** pure technological duty-cycle or actual measured sunshine.
- The speed gate combines data quality and true movement geometry.
  Removing it may alter the biological response distribution.
- The >=3/3 DAY/NIGHT daily and >=10-days/year tag gates depend on
  availability of observations in both dayparts and cannot be treated
  as independent of the animal's phenology, habitat or activity state.

The aggregate comparison will give *source-equal* medians and ranges
separately for already q-complete versus q-incomplete original sources.
These groups are selected **after** the first effects are known, and
any difference is **descriptive**. Source studies, not original fixes,
are the replication units; no event-level p-value or causal
unobservability model will be claimed.

## Source independence and source domain

The original 437 packages and 120 stage-one metadata candidates were
not a global probability sample of world wildlife populations. The second
wave quota was chosen after observing the first-wave result, although
the source identities followed the earlier frozen SHA ranks.
One `Homo sapiens` source is a permanent **human-location privacy HOLD**:
no provider file queries or GPS event access. Managed/domestic sources
remain source-accounted rather than removed for a favorable sign.
All 35 original source IDs and 50 registered taxa remain in the output
even if their original data never passed screening.

## Relation to the three existing RSDS methodological routes

- [PR #23](https://github.com/zuizui0223/RSDS/pull/23) identified
  between-study effect heterogeneity and missingness limits in nine
  supported older DOI, relying on approximate native percentile CIs.
- [PR #24](https://github.com/zuizui0223/RSDS/pull/24) estimated
  same-source native date/month joint covariance in four purposively
  selected studies, with a **date-equal** rather than tag-equal effect.
- [PR #25](https://github.com/zuizui0223/RSDS/pull/25) and
  [PR #26](https://github.com/zuizui0223/RSDS/pull/26) exposed the
  source selection process and native observation nonresponse among
  35 newly sampled sources, and showed a descriptive response-only
  `q` sign change from -0.04779 to +0.07568 when the response
  source composition changed.
- This route decomposes those sources' original UTC cadence and
  astronomic/geometric q eligibility *without changing any q*.

The old 48-source externally held-out **negative irradiance direction
was NOT REPLICATED**, and the older 12 source raw ERA5 checks remain
**HOLD**. None of those original scientific classifications may be
promoted by this source-technical diagnostic.

## Next interpretable source model

Only after the source-native technical funnel is externally verified
can one test a specific, prespecified association such as whether the
source-level fraction of 0.5–2h adjacent steps predicts source q
measurability **within the original 25 parsed studies**. Such a model
requires a distinct source-level response, source-held-out checking,
small-n uncertainty and documentation of partially observed
source-native measurement. It would still NOT identify missing
biological q for sources with unsupported sampling processes.

**Current claim ceiling:** a reproducible, original-GPS-validated
measurement-stage source selection audit, not a novel estimator of
animal chronotype, environmental causality, or global nighttime activity.
