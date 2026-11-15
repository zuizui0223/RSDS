# What changes after adjusting for month on *identical tagged animals and dates*?

**RSDS, 15 November 2026 | Post-outcome observational sensitivity | Fixed five original Movebank source studies.**

## Biological question and result

The prior five-original-study analysis showed that original species-native NIGHT−DAY GPS displacement `q` mean could retain a similar sign across studies while the **full-year** `q`–astronomical-daylength slope could have opposite nominal within-source individual-bootstrap 95% intervals. Year×month conditioning removed the apparent nominal opposite-signed seasonal source CIs. However, month admission changed source-tag/date support.

We therefore reran **the same five original source studies, with the original same-provider UUID/bytes/MD5 and source study/native lists**, replayed each old full-year `q` and photoperiod coefficient and the already frozen first within-month source coefficient, then calculated the following *on precisely identical original qualified tag-year-month dates and tags*:

- **A** = original-tag-equal source average of *unadjusted* original tagged year slopes `q ~ L` on the month-eligible sample;
- **B** = source average of original tagged year×month-demeaned `q ~ L` slopes on exactly the same sample;
- **Δ = B−A** = a within-source paired original-tag difference, with **the same 2,000 tag bootstrap draws** simultaneously used for A, B and Δ.

This controls the changing source-tag/date sample in the *paired* comparison. It does **not** treat time of year as an intervention, and it does not address source differences in weather, season, movement, duty cycle, location or independence of original GPS tagged populations.

**All five source first scientific q/daylength baselines and first within-month coefficients were replayed; all five same-tag/date analyses were estimable.** CI [first matched-support original GPS Action #43723252524](https://github.com/zuizui0223/RSDS/actions/runs/43723252524), [Git-frozen source-level result](../research/RSDS_FIVE_SOURCE_MATCHED_TAG_DATE_MONTH_DECOMPOSITION_FIRST_RESULT_20261115.json), [original prerequisite first month study](../research/RSDS_FIVE_SOURCE_NATIVE_TAG_YEAR_MONTH_FIRST_RESULT_20261115.json).

## All five original studies: paired source-native results

All source effects below are q/log1p GPS chord-speed change per 1 h astronomical daylength, not actual locomotor activity time. Intervals are **original independent tagged-animal bootstrap nominal 95%** source intervals, unadjusted for study selection or multiple comparisons.

| Species and original DOI | Eligible tags | A: unadjusted, same sample [95%] | B: month-FE, same sample [95%] | Paired Δ=B−A [95%] |
|---|---:|---|---|---|
| *Anser albifrons* agricultural goose `1.fk899541` | 37 | −0.01705 [−0.03960,+0.00538] | −0.01817 [−0.06567,+0.02588] | −0.00112 [−0.04358,+0.04642] |
| *Anser albifrons* Pannonic `1.46b0mq21` | 5 | **+0.02393 [+0.00988,+0.03930]** | +0.01159 [−0.04115,+0.04778] | −0.01233 [−0.05829,+0.02514] |
| *Odocoileus hemionus* Utah `1.711` | 1,647 | **+0.005538 [+0.004779,+0.006258]** | **+0.003152 [+0.001733,+0.004588]** | **−0.002386 [−0.003699,−0.001069]** |
| *O. hemionus* Wyoming `1.703` | 115 | **+0.005391 [+0.004176,+0.006750]** | **+0.008806 [+0.005764,+0.012241]** | **+0.003415 [+0.000750,+0.006288]** |
| *O. hemionus* Sublette `1.708` | 131 | **−0.004098 [−0.008620,−0.000995]** | −0.002149 [−0.006785,+0.002448] | +0.001948 [−0.002465,+0.006218] |

**First ecological reading:** Month conditioning does not support opposite seasonal `q` slopes as a general same-species pattern. Old Utah/Wyoming mule-deer source-month slopes stay nominal positive; new Sublette's negative source point is not distinguishable from zero after month FE. Both goose sources' month slopes include zero.

## The key new distinction: sample selection versus month conditioning

1. *Anser albifrons* agricultural goose: original full-year −0.02653 with negative source CI, but even the **same-month-qualified sample before adjustment (A)** has a CI including zero (−0.01705 [−0.03960,+0.00538]). Here **changing admitted tag/date support already removes the earlier nominal negative sign support**; month conditioning produces no independently supported paired Δ (interval includes zero).
2. *Anser albifrons* Pannonic: on the **same selected source tags/dates**, A remains nominal positive [+0.00988,+0.03930] while B includes zero; **Δ is not separately supported** (95% [−0.05829,+0.02514]). Thus a significance-classification change occurs with month adjustment, but its paired source-tag effect uncertainty is broad, especially given only five original eligible tags.
3. *O. hemionus* Sublette: even on identical month-qualified tags/dates, A has a nominal negative CI [−0.00862,−0.000995], B includes zero. But **paired Δ also includes zero**, so absence of B significance does **not** demonstrate a supported negative-to-zero ecological mechanism.
4. *O. hemionus* prior Utah and Wyoming: both retain **positive A and B original source intervals**. Their paired changes have **opposite** signs: Utah −0.002386 [−0.003699,−0.001069], Wyoming +0.003415 [+0.000750,+0.006288]. Month conditioning affects the numerical slopes in context-dependent directions even though nominal positive source orientation persists.

This is the most important scientific message: **a stable direction of recorded day–night displacement can coexist with context-dependent seasonal associations; however the original full-year “reaction norm reversal” is not convincingly replicated within months, and the loss of nominal sign support is not uniformly explained by month effects rather than observation support.**

## Fixed latitude and source overlap caveats

We repeated **the identical tagged-date sample** using astronomical daylength at each source's fixed median tagged-mean latitude to avoid within-month shifts in an animal's current latitude contributing to the time covariate. Paired difference signs and interval classifications remain the same for the five sources: the **two old mule-deer signed Δ results remain opposite**, while the other three remain intervals including zero. This does not independently control changing environmental/weather opportunity or dispersal location.

Distinct original study DOI ≠ distinct telemetry animal, population, site or independent seasonal experiment. The Utah and two Wyoming source files may share movements or sampling design. A tag-bootstrap interval is not a properly independent study-cluster or source-to-source interaction confidence interval; opposite source-specific interval signs are not the same as a formally tested slope interaction. All comparisons were selected after seeing prior source effects. **Nothing about daylight causal selection or adaptation has been established.**

## Reproducibility and fail-closed record

- [Original source-file and q/month first-result locks](../research/RSDS_FIVE_SOURCE_MATCHED_DATE_MONTH_DECOMPOSITION_CONTRACT_20261115.json)
- [Implementation: exact source replay + paired tag bootstrap](../scripts/rsds_five_source_same_tag_date_month_decomposition_v1.py)
- [Synthetic known-truth and source-HOLD tests](../tests/test_rsds_five_source_matched_month_v1.py)
- [Five original source first-result aggregate](../research/RSDS_FIVE_SOURCE_MATCHED_TAG_DATE_MONTH_DECOMPOSITION_FIRST_RESULT_20261115.json)
- [Execution, source-private telemetry, all CI success](https://github.com/zuizui0223/RSDS/actions/runs/43723252524)

The original GPS/astronomical q calculator's scientific thresholds, original five source DOIs, previous first output and provider original file UUID/MD5/bytes remain pinned. The original first source/native and month effects are *replayed before new data*. No original raw tracked animal IDs, exact GPS positions, date labels, original individual slopes or bootstrap samples leave ephemeral CI runners. Source groups remain five even if an original provider file fails.

The previous **48-source held-out solar-irradiance→heading signed effect stays NOT REPLICATED**, and the original **12 ERA5 source raw-SHA cases remain HOLD**.

## Ecological interpretation and next question

The strongest positive RSDS result is now **conditional, source-specific structure of time use** — recurring night/day q directions in some species and distinct seasonal association responses across original study contexts — rather than evidence of a universal q direction or a universal species-level astronomical photoperiod reaction.

A genuinely causal or broadly transferable ecological explanation would require predeclared original cohort independence, comparable observation schedules and location/monthly climate/phenology support in independent field sources. No claim about a globally adaptive seasonal response or solar-compass use should be made from these estimates.
