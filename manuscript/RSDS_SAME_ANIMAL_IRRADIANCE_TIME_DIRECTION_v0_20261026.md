# Do solar-relative movement directions and diel displacement covary with shortwave weather anomalies? Matched original ungulate GPS across two regional studies

**RSDS original-source ecological manuscript v0, 2026-10-26**  
**Exploratory AFTER earlier source effects. Not a preregistered causal solar response or world-random species generality test.**

## Abstract

Astronomical daylight, incident solar radiation and the Sun's bearing are correlated but distinct environmental dimensions. Original ungulate GPS analyses in Utah and Wyoming had shown recurrence of the daylength–night/day displacement association in mule deer despite regionally different unconditioned bearings relative to the Sun. These separate measurements did not establish whether *the same tagged animal* changes its time-dependent movement and travel direction with radiative weather anomalies. We therefore revisited both original public Movebank GPS datasets, reverified original provider MD5 bitstreams, and exactly reproduced all six previously qualified native taxon/group night–day and daylength estimates **before** linking new ERA5 hourly weather. Within the same tagged animals and date-eligible source gridyear, observed original GPS NIGHT−DAY log1p displacement rate and daytime solar-relative travel vector were reconstructed. For direction, modeled ERA5 shortwave and cos/sin(GPS bearing−solar azimuth) were residualized within original animal×calendar month×one-hour apparent solar-clock cells against Sun altitude, altitude squared, Sun bearing, daylength, temperature, wind speed and wind–Sun geometry. For temporal displacement, original daily q was related to the daylight-step radiation residual of the same animal and date, conditional on daylight, temperature and wind.

Both original field studies replayed perfectly. All **six study×native taxon** cohorts yielded a new coupled observing frame, with **434 source-tagged animals** supporting both source weather–direction and daily night/day-weather partial associations across the groups, 10 distinct weather gridyear requests and 12 selected original taxon×gridyear context analyses. In all six source groups, the animal-bootstrap 95% interval for the q/daypart shortwave slope included zero. Among solar-heading vector components, Utah bighorn sheep had positive source-conditional coefficients for both cosine (+0.01686 per +100 W/m², bootstrap interval [+0.00636,+0.02647]) and sine (+0.02917, [+0.01765,+0.04126]); Wyoming bighorn original-source means did not show intervals excluding zero. Within Wyoming bighorn, two original gridyear contexts yielded opposite conditional sine-associated signs even though pooled study effect was unresolved. These are *study-specific weather associations in selected GPS supports*, not familywise-corrected effects, causal solar-compass mechanisms, or evidence of between-region coefficients differing. Our findings expose a now-measurable but context-sensitive coupling problem: time-displacement and solar-relative direction need not present the same signed conditional shortwave signature in source observations. General inference requires new original independent field-study cohorts and sensitivity to environmental response versions and geographic routes.

## 1. Background and ecological hypotheses

Time allocation, movement speed and movement angle need not share one universal Sun-driven biological coefficient. Astronomical Sun elevation and azimuth are deterministic functions of original event time and location, whereas surface shortwave varies with clouds and atmospheric state. Movement direction is constrained by home ranges, foraging trajectories and migratory corridors. Two original public ungulate field studies showed that mule deer had similar positive daylength–observed night/day displacement gradients but dissimilar *raw* solar-relative movement orientation. The earlier 48-source original GPS radiation association did not replicate its directional sign in a separately frozen source cohort. Neither result alone yields a jointly measured individual radiation–time–direction relationship.

The present study asks (i) whether shortwave variation beyond modeled solar geometry and meteorological controls is associated with observed daily night/day speed displacement; (ii) whether that variation is associated with daytime cos/sin bearing relative to the actual Sun; (iii) whether original source animals with one association also exhibit the other; and (iv) what, if anything, recurs across two regional datasets in the same native taxa. These are *post-prior-result exploratory* questions, not blinded tests.

## 2. Original tagged-individual field material

The original files are public Movebank Data Repository provider originals, not effects extracted from prior published analyses. Two previously identified studies were held fixed before new weather:
- **Utah**, DOI **10.5441/001/1.711**, four native ungulates (pronghorn, wapiti, mule deer, bighorn), provider original **MD5 08200c9dbdf7989fb5a35467213b7e5a**.
- **Wyoming**, DOI **10.5441/001/1.703**, mule deer and bighorn, original **MD5 bd1c02734b6de3ddc9207fc1251edd66**, source original license **CC BY-NC 4.0**; no original animal points redistributed.

For all six original native taxon rows, the previously published q/night-day group means, year-specific q/daylength group slope means and admitted original tagged-animal counts reproduced to 1e−12 **before** any new ERA5 response. Original output files preserve provider data-integrity and the full source taxon frame.

After original one-fix-per-UTC-hour time-only thinning, native GPS intervals remain 0.5–2 h and at or below 60 km/h. A tagged-animal source day requires ≥3 valid DAY and ≥3 NIGHT steps and the original year needs ≥10 qualified dates. The response q is \(\log(1+\mathrm{median\ NIGHT\ GPS\ km/h})-\log(1+\mathrm{median\ DAY\ GPS\ km/h})\); this is **not** the fraction of latent animal activity conducted at night.

Within each original native taxon, at most the two most-supported 1° grid×original year contexts are selected using only distinct original qualified animals and their original valid q-days, *blind to new weather and heading-effect sign*, before any hourly ERA5 request. Across six original native group rows this yielded 12 context analyses and 10 unique verified original ERA5 grid×year requests. Each study×taxon remains one taxonomic original field-study context; context replicates are **not independent biomes**.

### Identical original animal/date support, but different movement responses

The direction response uses the **same individual and source-observation dates** as q, then only that day's original moving daytime GPS intervals with displacement ≥0.1 km, 0.5–2 h interval and Sun altitude ≥10° at original beginning/end/midpoint. Let GPS travel bearing be B, instant Sun azimuth A, altitude h, apparent-solar hour τ, original astronomical daylight L, ERA5 shortwave R and wind flow heading W. The two orthogonal direction components are cos(B−A) and sin(B−A). Neither is body-axis orientation.

Each original moving daytime interval is joined to the published ERA5 UTC-hourly surface shortwave (W/m²), 2m air temperature, 100m wind speed and meteorological wind FROM direction (converted to flow TO). Provider response SHA256 hashes are retained but raw weather requests and individual GPS remain ephemeral. All selected original moving daytime intervals in these source contexts obtained verified hourly ERA5 matches. Gridded weather is a surrogate for the animal's exposure, not microhabitat shade or flight-altitude wind.

### Prespecified exploratory conditional coefficients

For each animal, within calendar month × **one-hour apparent solar-clock** blocks, retain ≥6 independent timestamped original qualified daytime moving intervals from ≥3 solar dates, with total ≥40 steps and at least two months and fine-clock blocks. Within those blocks, the source covariates were demeaned. Both modeled shortwave and the cos/sin direction responses were residualized by linear least squares on solar altitude (linear/quadratic), cos/sin instantaneous Sun azimuth, daylength, temperature, wind speed and cos(wind flow TO−Sun azimuth). The result is a source-animal observational slope for cos/sin relative direction per +100 W/m² **residual shortwave**. Adequate conditional shortwave residual variance (SD≥20 W/m²), independent solar dates and rank support were required.

For the same animal and dates, after at least three eligible daytime directional-weather steps supported the shortwave residual, original daily q was regressed on mean residual radiation, with within-animal calendar-month centering and further daylength/temperature/wind controls. Minimum temporal support was 20 such days in at least two months and residual R SD≥10 W/m². Thus q and heading use **same source animal/date/exposure frame but different day/night and daylight moving sub-intervals**. Individual q and two heading coefficients are averaged across selected original contexts within the animal before equal animal-weighted source group summaries. Animal resampling 2,000 times gives **descriptive source-conditional 95%** bootstrap intervals; not independent probability sample population intervals or familywise testing.

## 3. Results

### 3.1 Original source reproducibility and actual coupled sampling support

Provider original bitstream checks passed for both studies. All four Utah taxa and both Wyoming taxa retained the **exact** previous q/daylength point results and original individual counts before new weather. The new fine-clock plus same-date/weather gates qualified **434 tagged-source animal groups** jointly for q×shortwave and cos/sin×shortwave slopes (Utah 99 pronghorn, 89 wapiti, 84 mule deer, 85 bighorn; Wyoming 54 mule deer and 23 bighorn). This is a count of tagged individuals across six original source×taxon cohorts, **not six independent ecosystems or 434 independent recruited populations**. More original source animals qualified for daytime directional weather effects alone than for the joint temporal outcome, notably Utah mule deer 180 directional animals but 84 paired.

### 3.2 Source conditional irradiation effects on time and direction differ

| Original source / native taxon | Same-tag paired N | Daily q × residual R (+100W/m²) | Daytime cos(B−A) × residual R | Daytime sin(B−A) × residual R |
|---|---:|---:|---:|---:|
| Utah pronghorn | 99 | −0.00173 | +0.00807 | +0.00315 |
| Utah wapiti | 89 | +0.00309 | −0.00239 | +0.00467 |
| Utah mule deer | 84 | +0.00452 | −0.00320 | −0.00336 |
| Utah bighorn sheep | 85 | −0.00283 | **+0.01686** | **+0.02917** |
| Wyoming mule deer | 54 | +0.00143 | −0.00737 | −0.00398 |
| Wyoming bighorn sheep | 23 | +0.00183 | −0.00652 | −0.01509 |

The 95% animal-bootstrap intervals for **all six q×radiation source means include zero**. Therefore the earlier strong q/daylength association cannot simply be restated as a uniform, similarly signed within-animal shortwave response after the declared adjustment. This does **not** establish absence of a temporal weather effect; the supported source context is a selected, conditionally observed fraction of original source animals.

For **Utah bighorn**, the original source animal summaries were positive for both solar radial cosine slope **+0.016863** (source animal bootstrap95% **[+0.006364,+0.026466]**) and tangential sine slope **+0.029168** (**[+0.017655,+0.041263]**), with 117 original bighorn animals eligible for heading slopes. Daily q×shortwave was **−0.002831** (**[−0.007521,+0.001377]**) in 85 eligible original sheep animals. These component results are exploratory after prior observation of solar data and across multiple taxa/axes. Their original animal bootstrap intervals are not adjusted familywise and do not identify a solar avoidance or orientation strategy.

In **Wyoming bighorn**, source-wide cosine **−0.006520** (bootstrap95% [−0.032300,+0.018927]), sine **−0.015089** ([−0.032962,+0.003215]) and q **+0.001831** ([−0.002104,+0.005453]) all had sign-uncertain source-animal intervals. The difference between source conditional bighorn cosine coefficients is not formally shown to be nonzero. In particular, it would be incorrect to interpret "not significant in Wyoming" as an independently confirmed region×radiation interaction.

### 3.3 Geographic/time-specific conditional movement geometry is heterogeneous within an original source

The fixed two q-support-ranked Wyoming bighorn original 1° gridyear contexts had **opposite source-conditional sine responses**: first **−0.06470** (source animal bootstrap95% [−0.09735,−0.03217]) vs second **+0.03686** ([+0.02385,+0.05022]). The source-wide equal-animal aggregation was near zero and uncertain. Contexts may differ in year, route, habitat, and share tagged animals; they are not independent field ecosystems or fully matched same-animal causal counterfactuals. This is a promising **contextual heterogeneity** pattern rather than a world mechanism.

The new paired source animal q-vs-cos-response Spearman correlations were +0.143 (Utah pronghorn), −0.070 (Utah wapiti), +0.125 (Utah mule deer), −0.046 (Utah bighorn), −0.350 (Wyoming mule deer), and +0.341 (Wyoming bighorn), with no valid world-clade or independent source population p-value. Thus the coefficients are *measurable on the same original tagged animals*, but their covariance is itself source- and taxon-specific and not yet externally confirmed.

## 4. Discussion and limits

This empirical stage resolves the earlier source-aggregation nonoverlap **within two named original GPS studies**: the same individual/date eligibility can indeed support both diel q response and Sun-relative vector response associated with shortwave anomalies. However, neither original source sample provides a uniform signed daypart response to adjusted shortwave. Utah bighorn displays a positive directional conditional association on both compass-relative vector axes that has no comparably signed source-level support in the second regional bighorn sample. Conditional direction effects themselves can reverse across original source contexts. These differences are biologically interesting precisely because they prevent a simplistic "same daylight means same direction" law; they do not isolate solar cues from geography.

The statistical structure is **within animal across weather samples → original source/taxon summaries → two original field studies**, not independent world-population/clade replication. Individual GPS activity was inferred from minimum observed straight-line step velocity, not time spent actively foraging, and weather fields are reanalysis products. The direction response is geometrically affected by reference-frame rotation, original migration corridors and source weather. Solar altitude and fine-clock nuisance controls reduce predictable geometry but cannot remove all residual time-of-day and habitat effects. Missing original individual q-days and source-gridyear selection can bias eligibility.

The distinct historical **48-source** solar-radiation negative-direction replication failure stays **NOT REPLICATED**. The present two-field-study q/Sun/ERA5 results are a different subset, different conditional estimand and exploratory outcome and **cannot retrospectively overturn** that external no-replication conclusion.

## 5. Next stronger biological tests

1. **Fixed-tag same-source context contrasts:** For taxon/source cohorts with both selected original contexts, measure within-tagged-animal *paired differences* of q/cos/sin weather coefficients (rather than only separate source-animal bootstrap means), with identical temporal-support and sampling rules, to determine whether the Wyoming tangential sign flip survives paired identity restriction.
2. **External biome/guild frame:** Select previously unexamined original public GPS study archives using only source taxonomy and sampling coverage, not favorable Sun/radiation effect signs; preserve source/camera detection differences and license restrictions. Fit the same observation estimand across new source populations, not pooled original steps.
3. **Exposure realism:** Compare true hourly weather provider version stability, actual local habitat shade/terrain aspect, rain/cloud state and if possible movement altitude to distinguish modeled radiation opportunity from solar signal choice.
4. **Complementary seasonal axes:** Bring back the original matched-photoperiod rising/falling-season comparisons where sufficient day/night support exists, but do not relabel unavailable original season comparisons.

## Files and provenance

- [Original pre-new-weather source contract](../research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_CONTRACT.json)
- [First original GPS and ERA5 source result in six native rows](../research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json)
- [Post-first ecological complete source synthesis](../research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_ECOLOGICAL_SYNTHESIS_20261026.json)
- [Original provider GPS read/MD5/replay + same-animal FWL](../scripts/run_rsds_same_animal_solar_weather_v0.py)
- [Original study ecological synthesis script](../scripts/analyze_rsds_same_animal_source_results_v0.py)
- [Known truth tests](../tests/test_rsds_same_animal_solar_weather_v0.py) and [results guard tests](../tests/test_rsds_same_animal_source_synthesis_v0.py)
- [Original first source replay and ERA5 field studies GitHub Actions](https://github.com/zuizui0223/RSDS/actions/runs/40274075085)

**Conclusion ceiling:** A conditional, source-specific positive solar-vector/weather association exists in observed Utah bighorn GPS and source-gridyear sign heterogeneity exists within Wyoming bighorn; all six source temporal q×radiation conditional group signs remain unresolved. No causal light-avoidance, mechanism-level temporal-direction decoupling, independent biome generality or confirmatory familywise inference is established.
