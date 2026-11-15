# RSDS: 437-package systematic source sampling and pre-event original GPS wave 1

**Frozen on 2026-11-15, after historical RSDS biological results but BEFORE
opening wave-1 original GPS event responses. Research development route,
not a prospective worldwide confirmatory meta-analysis.**

## Complete source accounting

The immutable public Movebank data-package catalog from **2026-10-19**
contains **437** packages. Exactly **259** had GPS/GNSS as an original
sensor label. A GPS sensor label does **not** establish an original
CSV schema, usable day/night intervals, independent field site, or even
the presence of compatible native animal records.

The full source-only metadata classifier preserves all 437 units with
explicit reason codes. Before original-file access:

| Metadata reason | Packages |
|---|---:|
| Not a GPS/GNSS package | 178 |
| Previously examined DOI/source-study ID including aliases | 14 |
| Repeated study-ID nonprimary package version | 2 |
| Missing/ambiguous study ID or DOI | 2 |
| Non-CC0 license; review not biological exclusion | 52 |
| Under first-wave reported animals/locations priority threshold | 66 |
| Missing/ambiguous taxon label | 3 |
| **Eligible for SHA-ranked first-wave priority** | **120** |
| **Total** | **437** |

The strict 120 are metadata-only candidates **after** the historical
13 DOI/source-ID gate, single source ID/DOI, CC0 first-wave license,
reported >=10 animals and >=10,000 locations, unique registered names.
The 120 consist of **109 single-native**, **11 multi-native** original
packages. Retain all metadata and its original SHA, including each
excluded/held candidate. Metadata records are not all 120 independent
field populations or equally usable GPS studies.

[Original source frame script](../scripts/rsds_437_systematic_source_frame_v1.py);
[full catalog metadata contract](../research/RSDS_437_SOURCE_FRAME_WAVE1_V1_CONTRACT.json);
[byte-identified first 16-source selection](../research/RSDS_437_SOURCE_FRAME_WAVE1_PRE_EVENT_SELECTION_20261115.json).

## Why a fixed 8+8 wave?

A deterministic SHA256 rank of a source-only salt, frozen original
catalog blob, original repository study ID and package UUID (not species
name or effect) chose **8/109 single-native** and **8/11 multi-native**.
The first original sample frame source/rank was committed before original
file metadata or new GPS events were requested. Under uniform hash rank
assumptions, metadata-**frame** inclusion fractions are 8/109 and 8/11,
respectively. They are NOT inclusion probabilities for all world
wildlife populations.

Multi-native sources were intentionally oversampled because their
**within-source cross-taxon covariance** is unknown from historical
percentile intervals. Consequently, equal weighting of the 16 rows
would misrepresent even the fixed 120-source metadata frame. A full
response conditional finite-frame contrast would poststratify source
means by 109/120 and 11/120. **With original file, schema, sampling and
ecological nonresponse, a respondent-only poststratified mean is not
design-unbiased.** The new aggregator reports it only descriptively
if both source strata respond, with explicit nonresponse and hypothetical
missing-value sensitivity.

The wave includes managed/domestic-associated taxa and genuinely wild
organisms. This is a source-metadata stratified *animal telemetry*
frame, not a claim of random wildlife taxa or ecological guild
representation. Do not silently remove managed populations after
seeing their q coefficients; predeclare wildlife/management ecological
moderators before future confirmatory source comparisons.

## Provider original file screening (completed before new event read)

The first source-only metadata result was checked byte-for-byte with
the original GitHub Actions artifact and committed as
[RSDS_437_WAVE1_ORIGINAL_FILE_METADATA_FIRST_RESULT_20261115.json](../research/RSDS_437_WAVE1_ORIGINAL_FILE_METADATA_FIRST_RESULT_20261115.json).
All 16 preselected source slots were retained:

- **13** had a unique metadata-selected original file with provider MD5,
  original bitstream UUID and byte count, under the fixed 200MB budget
- **3** remained original-file size or missing-size **HOLD**
- **0** animal GPS response rows opened at this stage

No replacement file, study, version or species was selected because
another original file exceeded a byte budget or showed an undesirable
ecological result. The first screening full JSON has SHA256
`a0f30a3318a2e52597201c4d0698efe3211401669178669ab3248a1fb3c0a826`.

## Precommitted original GPS effects

[Original source-first q and solar daylength contract](../research/RSDS_437_WAVE1_ORIGINAL_GPS_FIRST_EFFECTS_V1_CONTRACT.json).

The 13 MD5-qualified source files are next checked for provider
bitstream UUID, size, checksum and event schema. No outcome event is
accessed before those checks. Each original registered native label is
processed without effect-based species filtering.

- q: equal eligible original tagged animals, within animal equal
  year contributions, nightly median GPS straight-line displacement
  speed relative to daily median GPS displacement speed on log1p
  scale; NOT activity-duration fraction
- source native q~astronomical daylight per hour: separate
  within-animal/year observational association; not causal daylight
  stimulus
- exact original UTC hour-bin acceptance and 0.5–2h qualifying
  source-adjacent step constraints are reused from existing RSDS
- native means require >=5 eligible original tagged individuals;
  missing means stay **null**, not zero
- any multi-native shared-calendar covariance requires support for
  **every originally registered taxon**, matching source dates,
  >=5 native tagged animals, and existing >=5 month blocks/40 days
  eligibility; otherwise covariance HOLD without deleting other
  native taxa
- when available, 1000 common calendar-month block joint resamples
  plus optional original tag-positive reweighting generate a full
  source-native **2K x 2K covariance matrix**. The date-equal
  vector is an additional distinct estimand, **not** the older
  full-year tag-equal q effect

Any other original data schema, native-taxonomy, cadence, animal support,
or provider failure returns a source-specific HOLD with the original
source identity unchanged.

## Statistical synthesis and interpretation

One original source (DOI/repository study ID) is the study-level unit.
Up to 16 study groups may be available; the first 3 source-file HOLDs
remain among the 16 and do not count as zero. A source-level effect is
admitted only if **all registered native taxa in that source** have
estimable group effects. Partial native support is recorded separately
and never silently redefines the fixed original study estimand.

We distinguish:
1. Source metadata selection and file-screening probabilities.
2. Animal/GPS-event measurement eligibility and native-taxon
   support, possibly highly nonrandom.
3. Within-source native-covariance under shared-month resampling.
4. Between-source heterogeneity and transfer to entirely new original
   field studies.

Only after all four are defensible can one fit a multilevel
covariance-aware source/study/species synthesis and give an
out-of-source prediction range. Earlier [PR #23](https://github.com/zuizui0223/RSDS/pull/23)
uses imperfect source-tag percentile SE proxies; [PR #24](https://github.com/zuizui0223/RSDS/pull/24)
establishes original same-source native covariance for **four
purposefully selected** prior studies. The new source wave is an
explicitly selected different data frame and must not be pooled
with those effects until measurement estimands are aligned.

Every known previous 48-source negative sunlight–direction external
replication remains **NOT REPLICATED**, and the old twelve ERA5 raw-JSON
version checks remain **HOLD**. GPS activity, true ecological
chronotype, sunlight causal selection and global animal prevalence
are not identified by a repository selection protocol.

## Reproduction

Python 3.11+ from repository root:

```bash
python -m unittest discover -s tests -p 'test_rsds_437_systematic_source_frame_v1.py' -v
python -m unittest discover -s tests -p 'test_rsds_437_wave1_original_file_metadata_v1.py' -v
python -m unittest discover -s tests -p 'test_rsds_437_wave1_original_first_effects_v1.py' -v
python -m unittest discover -s tests -p 'test_rsds_437_wave1_original_gps_meta_v1.py' -v
```

The original GPS analysis runs in isolated temporary CI source jobs,
one per SHA-frozen package, with group-level JSON only uploaded as
an Actions artifact. Never publish original protected GPS events,
animal IDs, precise dates, coordinates or raw individual effects.
