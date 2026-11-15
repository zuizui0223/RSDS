# Five original studies: source-native q–daylength associations within tag × year × month

**15 November 2026 | Result from original Movebank GPS in GitHub Actions | Post-outcome sensitivity, not causal photoperiod evidence**

## Scientific status

Previous RSDS PR #37 found that *Anser albifrons* and *Odocoileus hemionus* could retain negative original-source **NIGHT−DAY GPS step displacement q** while their **full-year q–astronomical-daylength slopes** had nominal original-tag bootstrap 95% intervals with opposite signs across previously selected source studies.

**The new test is directly about that seasonal interpretation.** We reprocessed five preexisting, original provider-file UUID/MD5/size-frozen studies, used the exact old `one_taxon` original source/native estimator to replay full-source q and daylength slopes *before* any new month conditioning, then estimated an additional *within-original-tag × calendar-year × calendar-month* q/daylength slope.

All five original q/daylength first effects and supported original tag counts were reproduced to the fixed (10^{-10}) absolute tolerance. The five source studies, including the three original multi-native study source lists, were retained. No source failed at final analysis. First GitHub Actions full original-data result: [run #43722240065](https://github.com/zuizui0223/RSDS/actions/runs/43722240065). Frozen results are [source-only JSON](../research/RSDS_FIVE_SOURCE_NATIVE_TAG_YEAR_MONTH_FIRST_RESULT_20261115.json), [pre-new-month-effect scientific contract](../research/RSDS_FIVE_SOURCE_NATIVE_YEAR_MONTH_ADJUSTMENT_CONTRACT_20261115.json), [code](../scripts/rsds_five_source_native_tag_month_conditioned_v1.py) and [known-truth tests](../tests/test_rsds_five_source_native_month_v1.py).

## Central result — the proposed seasonal sign inversion is NOT retained as two opposite supported slopes

The source-month slopes below are **original tagged-animal-equal** means of each qualifying tag-year's within-month centered q~astronomical-daylength slope, averaged equally across qualifying years and tags, not raw GPS-event-weighted regression. Their original-tag bootstrap 95% intervals are *nominal* and **not** cross-study difference confidence intervals.

| Taxon / original DOI | Earlier full-year q–daylength slope (per hour) | New source-native within-month slope | Original tag-bootstrap nominal 95% interval | Qualified original tags |
|---|---:|---:|---|---:|
| *Anser albifrons* old grazing `10.5441/001/1.fk899541` | −0.026534 | −0.018166 | [−0.064963,+0.025732] | 37 |
| *Anser albifrons* new Pannonic `10.5441/001/1.46b0mq21` | +0.023928 | +0.011594 | [−0.041146,+0.047780] | 5 |
| *Odocoileus hemionus* old Utah `10.5441/001/1.711` | +0.005619 | +0.003152 | [+0.001777,+0.004565] | 1,647 |
| *Odocoileus hemionus* old Wyoming `10.5441/001/1.703` | +0.005538 | +0.008806 | [+0.005575,+0.012102] | 115 |
| *Odocoileus hemionus* new Wyoming Sublette `10.5441/001/1.708` | −0.003920 | −0.002149 | [−0.006737,+0.002646] | 131 |

**Main interpretation:** The apparent opposed full-year daylength associations in *A. albifrons* are no longer a pair of nominally opposite supported signs once within-month variation is considered; both intervals include zero. In *O. hemionus*, the two prior positive associations remain nominally positive within month; the new Sublette negative point slope still has a bootstrap interval overlapping zero. **There is no supported negative within-month slope opposing the two positive ones**.

The point signs themselves need not vanish when an interval includes zero. Neither post hoc month adjustment nor a widened interval can prove no month effect, and the month-eligible source tags/date support is not the same as the old full-year support. Current results therefore **limit**, rather than conclusively explain away, the full-year seasonal sign contrast.

## Fixed-latitude sensitivity

The exact same source-native month/tag/date support was also evaluated with astronomical daylength computed at the **original source-median tagged-animal latitude**, held constant across dates; no latitude or original coordinates are exported. This prevents changes in a migratory animal's actual latitude from directly entering this version of the daylength covariate.

The resulting five nominal source interval classes exactly agree with the local-daylength primary: **0/2 goose source intervals exclude zero; 2/3 mule-deer source intervals positive, 1/3 includes zero**. Representative means: old Wyoming mule +0.008741/h (fixed-lat bootstrap [+0.005492,+0.011823]); new Sublette mule −0.002825/h [−0.007364,+0.001971]. Similar nominal signs under fixed latitude do **not** control weather, tracking-duty-cycle confounding, breeding, migration stage or cross-source cohort overlap.

## Methods and original source integrity

- Five studies were selected **after** viewing old positive full-year patterns. Source identities, exact native taxonomy/DOI and provider original bitstream UUID/bytes/MD5 were frozen from older Git source blobs **before the new month results**. No source or tag outcome was substituted to improve signs.
- Original event definition unaltered: timestamp-only UTC-hour representative per original tag, adjacent 0.5–2 h valid steps, same local solar date and daypart at endpoints/midpoint, speed at most 60 km/h, at least three accepted DAY and NIGHT segments per tag/date; at least ten qualified dates per original tag-year.
- A tag-year-month must have ≥5 original qualified dates and ≥0.05 h daylength range; a tag-year needs ≥2 such months, ≥10 eligible dates and ≥0.05 h² centered daylength variance. At least five independent original tagged animals are required for a source mean. Resample source tags 2,000 times, not individual GPS events. Fixed thresholds were set before new month slope outcomes and unchanged after.
- The old Utah ZIP source initially failed **generic source CSV parsing**; the read path was fixed to the **previously used Utah native ZIP parser** without altering the original MD5, target species, observation criteria or sampling threshold. Old Utah source q and daylength means and tagged-animal count replayed exactly after this repair. This is an input-reader repair, not outcome-based biological reclassification.
- Raw GPS events, individuals/IDs, coordinates, exact dates, original tag-level slopes, bootstrap draws and source reference latitude remained confined to ephemeral Actions jobs. Only source/native means, intervals, independent eligible tag counts and aggregate technical support leave the runner.
- Source original field-study DOI/study ID may not equal independent tagged animal populations. Two Wyoming studies in particular require tag/provenance overlap review. Original tag-bootstrap intervals are **within-source**, without study-cluster, year-month, region or multiple-hypothesis correction.

## Scientific conclusion and next unresolved source-data question

A reasonable **observational** conclusion is that sign stability of source-recorded NIGHT−DAY movement (q) does not automatically imply identical seasonal full-year q–daylength relationships. However, the **stronger** hypothesis that two species show real, opposite, month-robust seasonal reaction norms **is not supported** in these five source-specific data conditions.

Next, separate two reasons for the apparent disappearance: (1) conditioning on year/month vs (2) the change in eligible original animal/date support. On the **same month-qualified original tag-year dates**, calculate both the date-unadjusted and month-FE tag slopes, resampling identical source tags jointly. Only then can a within-source paired difference be assessed descriptively; no ecological daylight causality will follow.

The original **48-source held-out sunlight/heading direction sign NOT REPLICATED**, **12 original ERA5 source raw SHA HOLD**, and earlier six-source frozen q outcomes remain unmodified.
