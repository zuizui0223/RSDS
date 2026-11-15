# RSDS Wave 2: moving from selected GPS studies to source-measurability meta-inference

**2026-11-15 — method route, separate from finalized first-wave science.**
Earlier [PR #25](https://github.com/zuizui0223/RSDS/pull/25) accounts for 437 Movebank public data packages, 259 with GPS/GNSS sensors, 120 stage-1 metadata-eligible source packages (109 with one native label, 11 with ≥2) and a SHA-ranked 8+8 first wave. The source-first first wave produced **4/16 source-complete q** and **3/16 source-complete q~daylength**, with zero new multi-native same-month covariance.

## Stage 1: locked second wave prior to any wave-2 GPS

A **new original selection contract** was committed after interpreting wave 1 but **before** opening any wave-2 provider-file metadata or event rows. The original 2026-10-19 437-package source frame and original SHA ranks are immutable. From those same ranks we took the **next 16 single-native source IDs (ranks 9–24 of 109)** and **the final three multi-native source IDs (ranks 9–11 of 11)**. This adds **19** untouched sources, for a combined **35/120** metadata candidates (24 single and all 11 multi-native). Remaining metadata candidates are 85 single-native.

The **quota** was chosen after viewing wave-1 results. The particular **identities** still follow the old source-only SHA rank and are outcome-blind with respect to wave 2. We must NOT describe the whole two-wave selection as prospectively fixed before any outcomes. No source replacement after newly observed q.

[SHA-pinned 19-source selection](../research/RSDS_437_WAVE2_19_SOURCE_SELECTION_PRE_EVENT_FIRST_V1.json).

## Stage 2: original file admission and privacy

The original 19 source package rows were first screened for provider bitstream UUID, original filename, bytes and MD5 without opening any original event. Exactly **15** had an original file <=200 MB with an MD5/bitstream ID and retained primary-source filename. **3** exceeded the frozen file cap or lacked qualifying bytes; they remain source-identity HOLDs.

A fourth source has the native metadata label *Homo sapiens*. It is assigned **HUMAN_LOCATION_PRIVACY_HOLD** before even retrieving provider file listings. No individual human trajectories, coordinates, timestamps or GPS source rows may be downloaded. Retain that package in the metadata screening denominator but it is not an animal-ecology measurement or a wildlife source. This flags a future domain-definition problem with the originally selected 120-package catalog, not a reason to swap in a more favorable animal study.

[Fixed 19 original-file metadata result](../research/RSDS_437_WAVE2_19_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1.json).

## Stage 3: unchanged original animal-GPS measurement

Fifteen original GPS file streams are processed under the established prior RSDS 1h source-timestamp selection and original adjacent 0.5–2h eligible step rules. The response is
`q = log(1+median observed NIGHT GPS km/h) - log(1+median observed DAY GPS km/h)`, on qualifying same-solar-date/daypart steps, with minimum original tagged animal support. This is **observed movement displacement**, not a latent daily activity budget or true nocturnality.

All **23 originally registered native labels** in the 19 wave-2 study packages remain represented, even if q, daylength associations, source common-date joint covariance or whole original GPS sources cannot be evaluated. The four provider/private-source HOLDs remain among the 19. For multi-native study-source covariance, all registered taxa must share the required original dates/month support; the paired calendar + animal-weight estimator is separate from full-source tag-equal q. Never claim unsupported covariance is zero.

The joint meta outcome of wave 1 + 2 must preserve **35 original source IDs** and the originally registered **50 native study–taxon records**. Species inside a common source DOI are not new independent populations.

## The estimands to report

**A. Source measurability**, not a biological effect: for each original source `s`, define `R_s^{(q)}=1` when all registered native taxa have supported original-GPS q; separately `R_s^{(L)}` and a stricter common-calendar covariance gate. Partition source failures into original provider file, original schema/rows, native sampling support, and joint calendar support. The first-wave gate was especially strict in co-registered multi-native sources.

**B. Finite 120-package metadata-frame measurability:** after wave 2, all **11/11** originally metadata-eligible multi-native source packages have an observed *gate outcome* (even if the gate failed), versus **24/109** single-native. If using SHA-prefix sampling as an approximate metadata-frame sampling process, report `(109 × observed single pass rate + 11 × observed multi pass rate)/120` only as a finite-archive exploratory projection. The SHA rank and the post-wave1 source-quota decision are **not a genuinely randomized sampling experiment**. Also report strict finite-frame logical bounds: at least the observed complete sources, at most observed complete sources plus 85 uninspected single-native sources. None of this estimates global animal populations.

**C. Source-complete biological q:** only admit a study q if **all pre-registered native taxa within that study** provide supported tag q. A partial multi-native source still has useful within-species descriptive q but does NOT become an independent single-species field-study record in the source-complete analysis.

**D. Response-missingness sensitivity:** list mean and variability for respondents in each source stratum, then a descriptive stratum-composition-adjusted mean (109:11) **only when each stratum has ≥1 observed complete q**. This is *not* unbiased for the 120 metadata frame while missingness depends on original GPS cadence/schema/native support. Display explicit hypothetical effects for incomplete sources rather than treating their q as zero, and do not present a standard random-effects pooled interval as global.

## Relation to earlier method streams

- [PR #23](https://github.com/zuizui0223/RSDS/pull/23): nine DOI source-equal heterogeneous q, source-native percentile SE proxies and conditional missingness sensitivity. Historical purposive selection, not the newly selected source frame.
- [PR #24](https://github.com/zuizui0223/RSDS/pull/24): original GPS joint shared-month covariance for four previously selected, co-calendar field studies. DATE-equal estimates cannot be pooled with original full-year TAG-equal q without explicit harmonization.
- [PR #25](https://github.com/zuizui0223/RSDS/pull/25): frozen source/candidate sampling and wave-1 original GPS first scientific outcomes.
- This wave: prospective **with respect to new wave-2 source events only**, source-accounted stage-specific measurability and native q estimates.

**Unchanged:** older 48-original-source signed sunlight-direction held-out negative external test **NOT REPLICATED**, older 12-source ERA5 raw-JSON **HOLD**, no general claim about causal sunlight, solar compass, biological competition or world prevalence.

## Route status

New source metadata was already frozen before any wave-2 original-GPS download; source/event job execution must validate provider MD5 and preserve private original telemetry within ephemeral runners. First combined 35-source meta is a newly declared endpoint, not a reclassification of archived first-wave receipts.
