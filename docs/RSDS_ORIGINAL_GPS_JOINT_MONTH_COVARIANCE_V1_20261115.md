# RSDS joint native covariance, original GPS — distinct post-outcome method route

**2026-11-15. Source-original field-study covariance extension.** This
route does not modify the 12-study source-equal first result or PR #23's
nine-DOI random-effects *summary-CI* sensitivity. It replaces the
uninformative guessed within-DOI correlation endpoints **only for four
historically preselected co-calendar source studies**, using joint
resampling of **original GPS observations** under strict MD5 and
historical-result replay.

## Source frame and comparisons

| Field study | Published DOI | Registered/native taxa | Existing identical calendar days | Original months |
|---|---|---:|---:|---:|
| Illinois | 10.5441/001/1.649 | 3 | 169 | 9 |
| Wyoming | 10.5441/001/1.703 | 2 | 193 | 8 |
| Utah | 10.5441/001/1.711 | 2 selected of 4 registered | 199 | 10 |
| Geese | 10.5441/001/1.fk899541 | 3 | 120 | 5 |

**Important:** Utah is a *selected pair*, not a complete four-native
source census. Those selections were outcome-aware and cannot be called
prospective. The original four-studies eight native pair contrasts are
replayed from SHA-pinned source results before interpretation. Each
source DOI is one study, not one per species or per monthly contrast.

## Estimands: share both the calendar and within-month blocks

For every source, derive daily
`q_{s,d} = mean_{tag observed on day d}(log1p(NIGHT GPS median km/h)
 - log1p(DAY GPS median km/h))` for all native taxa on the same
historically fixed calendar days. A datum is an original animal×date
qualified under the existing GPS cadence/solar-date/daypart rules.

Select month blocks with at least five original common-calendar days,
within-block daylight range of at least 0.05 h, at least five eligible
blocks, at least 40 qualified total days, and combined within-block
daylength sum-of-squares of at least 0.25 h². **Eligibility depends only
on date support/astronomy, never q effects.** For the *same* eligible
days, form a source-native joint vector:

```
theta_s = (
  mean_date(q_native1), ... mean_date(q_nativeK),
  within-month daylight slope_native1, ... daylight slope_nativeK
)
```

The slope is sum_blocks(sum_d (L_d - mean_block L) q_s,d) divided by
sum_blocks(sum_d (L_d - mean_block L)^2). All taxa use the exact same
sampling calendar and reference astronomical daylight (43°N for
Illinois/Utah/Wyoming, source-native position-blind study reference
for geese; the geese reference is never exported).

This is **DATE-balanced, not animal-balanced**; it is intentionally
distinct from the older full-native equal-animal-year estimates.
There is still no standardized observation-cadence interval across
source studies, so no worldwide pooled biological effect is allowed.

## Paired covariance — two transparent sensitivities

1. **Calendar-only:** 2000 draws of whole study year×calendar-month
   blocks with replacement. The same sampled blocks apply to every
   species and to both q and slope coordinates. This estimates
   conditional variation of fixed-tag source daily means and preserves
   shared seasonal covariance.
2. **Calendar + tagged-animal positive reweighting:** for each same
   month draw, draw one `Exp(1)` positive weight per observed original
   tagged animal, separately per taxon, and apply it to all dates for
   that animal. This is Bayesian-bootstrap-*style* uncertainty
   sensitivity, not a nominally calibrated hierarchical frequentist
   bootstrap. Positive weights ensure no previously supported calendar
   day disappears merely because its sole tagged animal was resampled
   zero times.

Estimate the full `2K × 2K` matrix of joint replicate covariance.
The resulting matrix is positive semidefinite by construction, and
species contrasts use `c'Vc`; no synthetic independence assumption
for species that share source months. Do NOT export exact dates, tag
identities, reference latitude for geese, or the private replicate
draw vectors. Only source-level matrices and summary support leave CI.

**Claim boundaries:** four purposively selected DOI, animal
GPS *displacement* not activity duration, joint same-date support
does not prove sympatry, calendar daylength cannot be experimentally
separated from season, strong measurement/schedule heterogeneity
remains. The estimated covariance is not an identified covariance of
independent future populations and should not be attached to the
older full-year equal-tag q estimates. Any unavailable original
provider stops with its original DOI and HOLD rather than source
substitution.

## Validation and gating

- Recheck original publisher source bitstream UUID / byte length / MD5.
- Replay previous original native full-source q and q~daylength
  coefficients (equal-tag, historic) to numeric tolerance `1e-12`.
- Replay identical original shared source calendar size, tags, and
  the 8 previously fixed pair within-month slopes to `2e-10`.
- Unit tests of joint covariance PSD, `c'Vc` pair contrasts,
  shared month resampling and fail-closed missing native source/day.
- Keep prior 48-source sunlight-direction **NOT REPLICATED**, and
  12-context original raw ERA5 version **HOLD** unchanged.

Run unit checks:
`python -m unittest discover -s tests -p 'test_rsds_original_joint_calendar_covariance_v1.py' -v`

Run source-only original GPS in an ephemeral Action runner:

`python -m scripts.run_rsds_original_joint_covariance_v1 --study Illinois --out outputs/illinois_joint.json`

Allowed source names: Illinois, Wyoming, Utah, Geese. Source
matrices will be artifacts, **not automatically reclassified as
biological external replication or a 12-DOI meta-analytic model**.

## Next quantitative synthesis

After source matrices pass original GPS reruns, extract a
covariance-aware data schema for `metafor::rma.mv` using full
source-native effect and `V` blocks. For now, the q/date and
month-daylight endpoints from **four selected studies do not comprise
a representative global ecological sample**. Further source
acquisition must be source-list prespecified, admit all unavailable
studies, and collect true jointly resampled animal *or* month blocks
under comparable GPS observation designs. Do not fit a high-dimensional
species×geography meta-regression with four studies.
