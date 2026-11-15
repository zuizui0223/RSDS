# Method: when observed animal GPS-source eligibility changes a meta-analytic sign

**2026-11-15 — simulation study, separate from original-source ecological evidence.**

The previous [RSDS original GPS observation funnel](../docs/RSDS_35_SOURCE_GPS_OBSERVATION_FIRST_FINDINGS_20261115.md) quantified technical attrition across 35 source-study IDs. A 0.5–2 h adjacent interval constraint removed 18.89% of 1,502,093 original accepted-hour adjacency opportunities; same astronomical-day/date/daypart continuity removed a further 12.31% of the remaining intervals; the 60 km/h step-speed filter removed only 0.17% of those solar-continuous intervals. **Source-complete q was available in 10/35 studies**, not all 25 parsed GPS studies.

Nevertheless, 25 observed parsed studies revealed original time-gate ratio **source AUC 0.707** for recognizing q-complete vs q-incomplete studies, while the final accepted interval fraction had **AUC 0.543**. These are *post-outcome within-sample source-level ranking diagnostics*, not validated external prediction. The number of eligible native original tagged animals is part of the q admission rule and cannot be used as a novel independent predictor.

## Core identification question

A summary-data meta-analysis computes an average `q` only from sources where `R_s=1`, meaning that a full registered-native time-niche effect could be estimated. The finite metadata-source-frame target `mu = E_frame(q_s)` is distinct from `E(q_s | R_s=1)`. The difference may arise from informative original observation design, unbalanced native support and source eligibility.

The **observed original GPS data do not identify `q_s` for missing sources**. Accordingly this route creates a **known-truth synthetic source frame** where latent `q_s` exists for every source, and explicitly varies the unknown link between true q and whether an original full-native source measurement could be admitted. This mechanism is hypothetical, not inferred from the missing real `q_s`.

## Frozen finite-frame synthetic design

- Exactly **109 single-native** and **11 multi-native** synthetic source populations, matching the *metadata* counts of the old fixed Movebank catalog; each is a source unit, not a GPS fix or an independent tagged animal.
- Known true source `q_s ~ Normal(+0.04, 0.12²)`, independent of synthetic technical timestamp cadence `c_s ~ Uniform(.35,.95)`. No animal data are simulated or used to estimate this distribution.
- The hypothetical selected samples are sequential `8/8` then `24/11` from independently permuted source IDs. Crucially, unlike the *actual* Wave2 quota (chosen after Wave1 results), the simulation fixes these quotas in advance. This is an **idealized mechanism illustration** only.
- Generate `R_s ~ Bernoulli(logistic(alpha_stratum + .9(c_s-.65)/.20 + beta_q(q_s-.04)/.12))`. The original source strata have different intercepts; `beta_q` is **0, -1.5, +1.5** as three frozen cases.
- Each of the three scenarios reuses exactly the same source true q, technical cadence, randomly selected source IDs and uniform response draws in each world. This paired design isolates the observation-selection mechanism.
- **1,000 independent finite 120-source worlds per scenario**. All worlds are retained; where a respondent-only estimator has no eligible source in one stratum it is **UNAVAILABLE**, never zero-filled or replaced. The known-propensity design estimator can still be evaluated using a zero observed contribution because, unlike biological q imputation, that zero is a lawful Horvitz–Thompson realized-sample term under a known inclusion law.

## Estimators, with strictly separated interpretation

1. **Truth:** finite frame all-120 source average of known q. The target is recalculated in each synthetic world.
2. **All 35 selected q observable (oracle):** poststratified sample mean using 109:11 source-frame weights; quantifies ordinary source selection variability absent measurement losses.
3. **Respondent-only q mean (naive):** poststratify the sources with `R=1`. This is not design-unbiased for all 120 when the response probability depends on unobserved q.
4. **Oracle Horvitz–Thompson:** use true conditional `P(R|q,c,stratum)` in each selected source's design weight, with all 1,000 worlds retained (even if a stratum has zero respondents). This is **impossible to implement in actual GPS archives**, because `q` is precisely missing where `R=0`.
5. **Oracle cadence-only Horvitz–Thompson:** integrate source q from the *synthetic data-generating distribution* to know `P(R|c,stratum)`. With informative missingness depending on true q, knowing this marginal propensity does not identify `E(q|R=0,c)`, even if cadence is perfectly measured.

Compare estimator bias and RMSE against the true finite source frame mean over **independent simulation worlds**. Report Monte Carlo standard error and intervals of the *simulated mean bias only*, not CIs for any wildlife biological effect. Show the fraction of simulated worlds where respondent-only q has the wrong sign, and the fraction where Wave1 respondent-only q is negative and combined-wave q positive. The latter is a **possible sampling consequence**, not the probability of the actual observed Wave1→Wave2 sign reversal.

## Interpretation / manuscript position

The contribution is an inspectable **failure-mode demonstration connecting original GPS observation architecture → source measurability → heterogeneous meta-estimands**, not new theory for inverse-probability weighting, random-effects meta-analysis, imputation or causal inference. Simulation parameters and data-availability dependencies are illustrative and must never be described as estimated from the original Movebank studies.

The empirical 35-source run establishes *where* original data were filtered and how little total interval acceptance alone discriminates source-complete q. The known-truth simulation shows *why* even apparently reasonable source metadata or technical-cadence adjustments cannot generally restore missing biological q if missingness depends on the unseen effect itself.

### Required next empirical validation

To estimate real source comparability, specify an observation design that can recover q without changing the historical estimand, e.g. an outcome-blind subset of original source-cadence protocols, predeclare source and species inclusion, and validate on **independent new repository DOI/study IDs**. Do not promote the synthetic oracle to a practical adjustment for existing 35 sources.

**Invariants:** previous 48-source sunlight-direction held-out sign remains **NOT REPLICATED**, original 12 source ERA5 raw SHA versions remain **HOLD**, and historical 35 source/time q labels, original native IDs and all protected GPS human-source HOLDs remain unaltered.
