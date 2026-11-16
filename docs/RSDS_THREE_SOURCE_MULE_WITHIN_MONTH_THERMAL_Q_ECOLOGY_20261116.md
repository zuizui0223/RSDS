# RSDS: Same-species thermal associations persist within original calendar months across three mule deer tracking studies

**16 November 2026 | Original Movebank GPS and current ERA5 physical snapshots, three original source DOIs, post-discovery exploratory analysis.**

## Main biological discovery

Across the three previously fixed public original *Odocoileus hemionus* field-study archives, **the association between daily ERA5 maximum air temperature and observed night-relative GPS displacement was positive in every original source, even after controlling for local astronomical daylength, antecedent 30-day precipitation minus FAO reference ET0, and original tag × year × calendar-month**.

The astronomical-daylength coefficients, meanwhile, were not consistently signed in the very same model. This suggests that **a temperature-linked component of the observed diel displacement niche is more directionally repeatable across the three source contexts than the daylength-linked component**, not that temperature was experimentally shown to alter animal behavior.

Here `q=log1p(median observed NIGHT original GPS straight-line step km/h) - log1p(median observed DAY original GPS straight-line step km/h)`. A positive slope with `Tmax` means q *increases* with modeled daily heat on the analyzed source dates; it does **not** mean all animals are genuinely more active at night, or that average q is positive. The frozen mean q was negative in all three original sources.

### Original source-native, tag-year-month controlled full-joint model

Original individual/tag bootstrap nominal 95% percentile intervals, **equal original tagged animals**, from source-level original data with fully checked q and daily climate snapshots:

| Original source (same species) | Qualified original tags | **Tmax coefficient (per °C)** | **Daylength coefficient (per h)** | **P−ET0 coefficient (per 100 mm)** |
|---|---:|---|---|---|
| Utah `10.5441/001/1.711` | **1,596** | **+0.001087 [+0.000974,+0.001199]** | −0.000951 [−0.002672,+0.000812] | +0.008928 [+0.004351,+0.013852] |
| Wyoming `10.5441/001/1.703` | **109** | **+0.000571 [+0.000240,+0.000914]** | +0.008663 [+0.005045,+0.012251] | +0.009566 [−0.003167,+0.020904] |
| Sublette `10.5441/001/1.708` | **129** | **+0.002046 [+0.001700,+0.002398]** | −0.020726 [−0.026716,−0.014470] | −0.039143 [−0.057250,−0.021673] |

All 3 original source studies were estimable; no individual provider weather cell was dropped and every source's **2026-11-16 original daily ERA5 physical SHA256** was exactly reproduced. Before these new month models, each original provider GPS UUID/bytes/MD5, the old native full-source q/daylength mean+number of tags, and all 6 prior original 2026 daily q/weather source effect means+tag counts were independently replayed.

**Do not calculate a pooled animal-level 3-source p-value**: old Utah and two Wyoming sources are not demonstrated to be three statistically independent ecological populations. The three original source-level bootstrap intervals are not a formally calibrated meta-analysis or a test of between-source effect sizes. They are three observational source-specific nominal 95% interval descriptions.

### Why same-tag/date comparison matters

The previous 3-source daily model compared q~L|Tmax+water and q~Tmax|L+water across all day-qualified original tagged animals, but residual temperature/solar effects could represent changes in study-month mix.

The new model was fitted within each individual-tag × calendar-year × calendar-month. In every original source, the **raw calendar-month-unadjusted joint regression** was also calculated on the **identical accepted tag/year/month/date/weather records** and resampled with the exact same 2,000 original-tag bootstrap draws. A source cannot report one of the three joint coefficients unless the full 3-variable model is identifiable on the same original tagged-year sample.

**Paired Tmax change: within-month adjusted minus unadjusted on identical source tags and dates**

| Original source | Paired Tmax Δ per °C | Nominal original-tag bootstrap 95% interval | Interpretation |
|---|---:|---|---|
| Utah | −0.0000166 | [−0.0000478,+0.0000158] | no supported change |
| Wyoming | −0.0000164 | [−0.0000987,+0.0000654] | no supported change |
| Sublette | **+0.0003356** | **[+0.0001630,+0.0005059]** | positive change, not attenuation |

Thus **the positive original source temperature association did not vanish when removing between-month contrasts**. For two sources its absolute coefficient was effectively unchanged within paired source tag/date sampling; in Sublette it increased nominally.

The within-month photoperiod source effects still differ. Utah includes zero; Wyoming is nominal positive; Sublette is nominal negative. However the paired month-minus-raw *daylength* contrasts themselves all have 95% intervals that include zero, so one must not claim a fully established between-source photoperiod reaction norm. In the earlier q~L-only month adjustment (without weather controls), Sublette's negative was not nominally supported; the new joint model is a different estimand, and such conditional coefficients can change under correlated climate/season covariates.

The within-month climate **water balance** association is heterogeneous: nominal positive for Utah, unresolved for old Wyoming, nominal negative for Sublette. This is a lagged **modeled climatic water balance**, *not* measured vegetation biomass, body water, available forage or a hydrological treatment. No general positive/negative water rule follows.

## Data and strict reproducibility gates

- [Pre-new-outcome source/model/version contract](../research/RSDS_THREE_MULE_WITHIN_MONTH_HEAT_CONTRACT_20261116.json)
- [Previous exact three source native original GPS+ERA5 result (unchanged)](../research/RSDS_THREE_MULE_ORIGINAL_DAILY_HYDROTHERMAL_FIRST_RESULT_20261116.json)
- [First completed original month-weather analysis JSON](../research/RSDS_THREE_MULE_WITHIN_TAG_MONTH_HEAT_FIRST_RESULT_20261116.json)
- [First successful whole Actions run #43874362941](https://github.com/zuizui0223/RSDS/actions/runs/43874362941)
- [Code](../scripts/rsds_three_mule_within_month_joint_weather_v1.py), [all-source summary](../scripts/aggregate_rsds_three_mule_within_month_weather_v1.py), [known-truth synthetic and source-HOLD tests](../tests/test_rsds_three_mule_within_month_weather_v1.py)

Original independently tracked GPS source data were loaded only on ephemeral GitHub CI runners. All registered original source DOIs and species were retained, without selecting favorable source dates on q signs. On each original tag-year, months needed ≥5 original source-qualified **weather-complete** q dates and ≥0.05 h observed astronomical daylength span. Eligible tag-year needed ≥2 such calendar months, ≥15 original tag dates and ≥0.05 h² of within-month centered daylength variation. FWL partial OLS needed independent residual variation in Tmax ≥0.25°C, daylight ≥0.05 h and 30-day P−ET0 ≥0.05 (100-mm units), with at least 5 qualified original tags per source. Multiple source-years per tag and source tags were averaged equally; 2,000 source **tag**, not GPS-event, bootstrap resamples.

New source weather physical content had to match the immediately preceding new 2026-11-16 hourly/day-level snapshot fingerprint **exactly**. A mismatch would have been a new-version HOLD rather than an explanation for changed signs. No raw GPS point, original individual ID, exact date, private source 1° cell identifier, individual coefficient or bootstrap draw was exported from the CI runner.

## Scope of the claim, and next ecological mechanism test

**Supported observational description**: This selected species shows consistently positive *within-original-tag-month* modeled-temperature slopes for night-relative GPS displacement across all three existing field-source studies, even with daylight and 30-day climatic water balance included. Daylength and water slopes are less directionally repeatable in this sampled design.

**Still not supported**: Temperature is a causal *thermal avoidance* mechanism; animals actually become nocturnal; changes in vegetation productivity, forage availability, reproductive stage or predator risk are mediators; original tags or study populations are definitely independent; three studies establish a worldwide mammal law or prospective generalization. All target species, studies and month models were selected after earlier positive observational patterns. GPS interval/duty-cycle variation may also affect q.

**Next mechanism-oriented distinction**: Is the positive q–Tmax association driven by *decreasing observed daytime movement*, *increasing observed nighttime movement*, or both? Estimate the separate DAY and NIGHT original `log1p(median GPS step-speed)` components on the identical source tag/month/date/weather opportunity, not just their difference. A stronger ecological explanation requires independent source populations, schedule-robustness and direct heat/movement corroboration, ideally external source studies selected without seeing their q or thermal slopes.

Previously negative held-out 48-source irradiance→heading directional sign replication remains **NOT REPLICATED**. Earlier 12 ERA5 **raw JSON SHA HOLD** contexts remain unchanged and are not retroactively validated by these newly frozen physical weather versions.
