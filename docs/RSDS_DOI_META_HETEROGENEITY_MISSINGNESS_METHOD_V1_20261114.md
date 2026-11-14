# RSDS methods v1 — study-level meta-analysis sensitivity, not a global pooled effect

**14 November 2026; branch-only post-outcome extension.** Parent:
[12 original Movebank studies, 22 native taxon rows](../research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json),
original Git blob `af0574c6bb6383402a97b5bebadb8cf0408e9740`.
The parent, its old external 48-source solar-radiation **nonreplication** and
the original 12 ERA5 raw-JSON **HOLD** are immutable.

## Main methodological result

The target is **heterogeneity / sensitivity to ecological and observational
study architecture**, *not* a claim that all wild animals prefer night or day.
There are 12 source-field-study DOI (not 22 independent studies); 22 native
study-taxon records; 15 eligible q and photoperiod records; **9 eligible DOI**.
A DOI may contain more than one native taxon, and a source can be entirely
unavailable for capacity, taxon-schema or cadence reasons. Missing is not zero.

The frozen older DOI-equal descriptive q mean was **+0.0102623**, while the
unbalanced native-taxon-row mean was **-0.0225097**. This sign reversal
depends on the unit of analysis, motivating an explicit dependence analysis.

## Diagnostics implemented

`python -m scripts.rsds_doi_meta_heterogeneity_missingness_v1 --out /tmp/rsds-meta.json`

1. **One source = one observation.** Average only the estimable native-taxon
   point effects inside a source DOI. Preserve 12 DOIs and 22 rows in the
   denominator and output; fail closed on duplicate or modified identities.
2. **Sampling precision is only approximate.** Each source animal-percentile
   bootstrap interval `[L,U]` supplies a *proxy* `se=(U-L)/3.92`, not a
   published normal-likelihood standard error or joint bootstrap. For
   studies with multiple taxa, vary the assumed within-DOI equicorrelation
   between **rho=0** and **rho=1**; no naïve taxon-level independence.
3. **Exploratory random-effects diagnostics.** Calculate Paule–Mandel
   between-DOI `tau²`, modified Hartung–Knapp approximate interval with
   `t(df=8)`, and approximate prediction interval from `sqrt(tau²+SE(mu)²)`.
   Report prediction spread *alongside*, not beneath, the mean effect.
   These are **model-based sensitivity intervals** conditional on
   purposively selected studies; they do not have calibrated world-level
   coverage. Their high heterogeneity can reflect GPS cadence, tracking
   support, habitat, taxonomic composition and measurement architecture,
   not only true ecological variation.
4. **Two missingness tipping analyses.** Compute the common effect required
   across the three wholly missing DOI for the twelve-DOI equal mean to cross
   zero. Separately fill the seven missing native rows with either the
   minimum or maximum of the fifteen observed native point effects, then
   recalculate equal *registered* taxon means within each of the 12 DOI.
   The observed-range condition is an **illustrative assumption**, not a
   constraint that missing biological effects must satisfy.

## First calculations (conditional, not confirmatory)

| Parameter | rho 0 | rho 1 |
|---|---:|---:|
| q random-effects pooled mean | +0.00936 | +0.00941 |
| q tau (between DOI SD) | +0.10198 | +0.10194 |
| q model 95% interval | [-0.06935, +0.08807] | [-0.06929, +0.08811] |
| q model 95% prediction interval | [-0.23862, +0.25734] | [-0.23849, +0.25731] |
| q~daylength random-effects mean, per hour | -0.00230 | -0.00214 |
| q~daylength model 95% interval | [-0.01270, +0.00810] | [-0.01269, +0.00840] |
| q~daylength model 95% prediction interval | [-0.02970, +0.02509] | [-0.02949, +0.02520] |

The intervals are driven by only nine (nonrandom) original study sources.
Point-estimate heterogeneity is large; **no universal q sign or a universal
photoperiod slope** is supported. Both outcomes are different estimands
(`q`: observed GPS night-minus-day log1p straight-line speed;
`q~daylength`: within-tag/year association per hour of daylength). They
must not be combined as if they were two measures of one latent response.

**Missingness:** three wholly unestimable DOI with a common q of
**-0.030787** would already set the *all-12-study* equal q mean to zero
(holding the 9 observed DOI means fixed). If all seven missing native
effects were hypothetically within the observed native q range
[-0.16547, +0.18975], the equal full-frame DOI q mean would lie between
approximately **-0.03932** and **+0.06429**. Neither is a formal
partial-identification bound without the range assumption. The daylight
tipping value for the three missing DOI is **+0.01915/h**.

## What makes an actual meta-analysis stronger next?

**First priority: common estimands and within-source uncertainty.** Reopen
licensed original GPS only under source permissions and predeclared inclusion
rules, apply one uniform observation schedule or explicit cadence-bias model,
and calculate (i) source-level native taxon effects on the same clock/solar
time definition, (ii) *joint animal/bootstrap replicate vectors* and the full
within-source covariance matrix, and (iii) matched calendar-day contrasts.
Do not mix full-record annual slopes with shared-day / within-month slopes.

**Second: a prespecified systematic source frame.** Inventory *all* eligible
public tracking studies with full screening reason codes and nonresponse
mechanisms; lock inclusion before obtaining outcomes. The currently 12
purposefully sampled DOI do not represent a global ecological study sample.
No funnel plot, Egger test, or environment/guild meta-regression is
interpretable as confirmatory with this selected set and 9 supported DOI.

**Third: genuinely nested meta-analysis.** Once joint covariance and more
independent studies exist, use `metafor::rma.mv` with structured V, random
source and species terms, and `clubSandwich` CR2/small-sample df,
alongside robust sensitivity to unknown covariance. This is **planned**, not
claimed as an executed CR2 result; the existing metadata cannot identify
the covariance matrix.

**Fourth: prospective out-of-study replication.** Freeze geographic,
phylogenetic, seasonal/cadence moderator hypotheses and test them on
entirely untouched *new* published GPS sources. Do not promote the
old exploratory sign tests or the failed 48-study solar negative effect
to preregistered external successes.

Methods references:
[Cochrane Handbook chapter 10](https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-10);
[metafor dependent-effects and CR2 guidance](https://cloud.r-project.org/web/packages/metafor/metafor.pdf).

## Strict evidence boundary

This is a **post-outcome sensitivity analysis**. No novel independent
animal observations were analyzed. No field mechanism, photoperiod
causality, solar compass, competition, or cross-world prevalence is
identified. The 48-original-source negative sunlight-direction external
replication remains **NOT REPLICATED**. The source's old ERA5 raw SHA holds
remain **HOLD**. Preserve all published nulls, HOLDs and preexisting
scientific contracts.
