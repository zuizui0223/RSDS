# Temporal regularity without invariant solar-relative direction: a cross-axis macroecology of animal movement

**RSDS / original-data ecological manuscript v0 / 2026-10-26**
**Status: exploratory quantitative synthesis of original public tagged-animal GPS source-derived data, not a prospectively confirmed global macroecological effect.**

## Summary

Astronomical daylight and solar position are predictable boundary conditions of animal behavior, whereas realized temporal and directional movement depend on weather, geography and biology. Whether the **time-of-movement** dimension and the **solar-relative direction-of-movement** dimension preserve comparable structure across populations is unresolved. We used an original-data GPS research archive in RSDS to reconstruct six native source-study × taxon outcomes from two independent regional ungulate tracking studies (Utah and Wyoming) using exact public study DOI plus original taxon identity, keeping source-specific GPS eligibility and the distinction between the axes. We evaluated the original individual-year slope of night-versus-day GPS displacement contrast against astronomical daylength alongside the original mean cosine of GPS travel direction relative to instantaneous Sun azimuth, decomposed into geographic heading marginal and residual joint circular covariance.

In *Odocoileus hemionus*, the daylength-associated night–day GPS displacement gradient remained positive in both study-source cohorts (+0.005619, +0.005538 per additional daylight hour, both original source-animal bootstrap intervals excluding zero). However, unconditioned source mean solar-relative bearing reversed from +0.023393 (Utah) to −0.048758 (Wyoming), with both original source animal-bootstrap intervals excluding zero. The unconditioned mean difference (+0.072151) is an exact sum of a geographic-heading marginal difference (+0.045886) and joint circular component difference (+0.026265), demonstrating that geographic bearing composition cannot be ignored when comparing relative-to-Sun travel. The *Ovis canadensis* q–daylength effects were unresolved in both studies, and the Wyoming raw solar-relative mean was also unresolved. Relative nighttime displacement ranked higher in mule deer than bighorn in both original study samples, but Wyoming had no admissible same-grid/month sympatry support.

Across nine native taxa with two original solar GPS archives each, raw solar-relative cosine sign agreed in six; joint circular-covariance sign agreed in only four. A separate 48-study irradiance–solar-direction frame had only 20 identifiable radiation coefficients and failed its preregistered external negative-direction replication criterion. Its original study DOIs did not overlap either Utah/Wyoming qualified daylength archive and did not overlap any seasonal-branch study with an estimable original effect, so source-level marginal associations cannot identify a within-animal irradiance × temporal niche coupling.

These observations motivate a falsifiable ecological distinction: **astronomical-calendar-linked relative daytime/nighttime movement can show regional recurrence while actual solar-relative travel direction changes with local spatial opportunities and source sampling**. The current study does not identify causal solar cues, photoperiod memory, wind/thermal mediation, or independence of temporal and directional behavior within an individual.

## Introduction and hypothesis

Solar azimuth, solar altitude, local apparent solar time, daylength and rising/falling seasonal branch are related but non-equivalent time geometries. Air temperature, shortwave radiation, local wind and plant/water seasonal conditions need not be synchronized with any fixed photoperiod. Treating a source's observed solar-relative direction as a proxy for activity timing assumes a coupling that has not been tested.

We distinguish:
- **H_time (temporal recurrence):** Within a taxon, a defined relation of NIGHT−DAY GPS displacement to daylight can be present in multiple independent original field datasets, conditional on the source sampling protocol.
- **H_direction (spatial contingency):** Source-level cosine of movement bearing minus solar azimuth is not generally invariant across geographically separated original cohorts and can be partly organized by different geographic-direction marginals.
- **H_environment (environment-dependent coupling; untested):** Within matched individuals, activity-daypart movement and solar-relative travel covariance respond differently to shortwave anomalies beyond solar altitude and apparent solar time; the pattern could vary with wind and migration/foraging mode.

The paper's exact ecological target is **what is recurrent in the time-use component, what is geographically contingent in the directional component, and what requires joint observation for generalization**.

## Data and sampling boundaries

The parent RSDS research migration retained 390 source-specific ecological code/data/manuscript/contract files with identical original Git blob SHA. We do not analyze published effect tables. The exact input snapshots for this paper are pinned by the SHA manifest in the accompanying contract and output JSON.

Two original multi-ungulate regional survey DOIs provided identical taxon labels *Odocoileus hemionus* and *Ovis canadensis*, as well as two additional Utah taxa. A native time-response table had six study×taxon rows. The original solar GPS matrix had 69 original source×taxon records, 41 analyzable direction-cohort rows and 28 unsupported/HOLD source×taxon records. Join was exact DOI+taxon, never region-name similarity. Each original study counts as one source unit: two taxa from a common source are not two independent field ecosystems.

### Temporal observable

\(q_{i,d} = \log(1+\widetilde{v}_{night,i,d}) - \log(1+\widetilde{v}_{day,i,d})\), in original GPS step km/h. For animal and calendar year, the original slope \(\beta^L_{i,y}\) describes association between \(q\) and modeled astronomic daylength hours. It is not a fraction of behavior truly allocated to night, and the original 30–120 minute interval/duty-cycle rules selectively admit GPS steps. Within source, individual year and individual aggregation follows the original archived source method; uncertainty is conditional source-animal bootstrap.

### Heading observable and exact decomposition

For actual GPS travel bearing \(B\) and instantaneous calculated Sun azimuth \(A\), raw \(Y=\cos(B-A)\). The circular decomposition at original source is:

\[
E[\cos(B-A)]
=\underbrace{E[\cos B]E[\cos A]+E[\sin B]E[\sin A]}_{\text{geographic heading–Sun marginal}}
+\underbrace{\operatorname{Cov}(\cos B,\cos A)+\operatorname{Cov}(\sin B,\sin A)}_{\text{joint circular component}}.
\]

Neither raw nor joint components alone establish that individuals actively orient by Sun, because availability of geographic headings, migration corridors, wind fields and time windows are not randomized. The archive also records an original within-solar-clock circular covariance estimate as a *different* conditioned estimand.

### Descriptive repeatability rule

Source signs and their original animal-bootstrap intervals are shown as source-conditional observational support, **not** a new multiple-comparison inferential test. A repeated sign in two regional source studies is not a world-population law or effect equivalence. Multiple source taxa with shared study or potentially shared recruited individuals are not independent world ecosystems. Signed source differences and their algebraic marginal+joint decomposition are deterministic, not mediated fractions.

## Results

### Repetition of diel daylength association with geographical reversal of solar-relative travel

| Study / taxon | q–daylength (per hour) | Original descriptive 95% CI | Mean cos(B−A) | Original descriptive 95% CI |
|---|---:|---:|---:|---:|
| Utah / *O. hemionus* | +0.005619 | [+0.004388, +0.006921] | +0.023393 | [+0.016264, +0.030673] |
| Wyoming / *O. hemionus* | +0.005538 | [+0.004102, +0.006973] | −0.048758 | [−0.068825, −0.031284] |
| Utah / *O. canadensis* | −0.003412 | [−0.006503, +0.000612] | +0.020778 | [+0.012582, +0.028674] |
| Wyoming / *O. canadensis* | +0.000217 | [−0.005847, +0.010986] | −0.005281 | [−0.038617, +0.024402] |

Within *O. hemionus*, q–daylength slope Utah−Wyoming was +0.00008117; a small point-estimate difference is **not** a confidence-supported statistical equivalence test. Raw sunward cosine Utah−Wyoming was +0.07215108, exactly +0.04588622 geographic marginal +0.02626486 joint circular component. The 9-taxon two-source audit showed 6/9 raw signs and 4/9 joint signs agreeing, a descriptive signature of context-sensitive geometry and/or behavior.

The q night−day difference between mule deer and bighorn was +0.0439119 in Utah, +0.0155126 in Wyoming. In the latter independent Wyoming study there are zero original qualified shared 1° cell × calendar year-month survey units, so these are *source-group ranks*, not direct competition/sympatry.

### Radiation is a distinct, unreplicated and currently nonoverlapping axis

An independently frozen second source collection did not replicate the initial negative shortwave–direction observational association. Source-equal mean responses per +100 W/m² were −0.019659 (first identifiable N=10) and −0.003361 (external identifiable N=10). All 48 original archives are retained, including 28 undefined shortwave effects. **Zero** of the six temporal–heading survey rows' two study DOIs are in the 48 radiation source DOIs. Two older seasonal-branch *source metadata* DOIs occur in that radiation frame, but both had no identifiable original seasonal-branch response. Therefore a source- or species-aggregated regression of radiation effects on temporal niche response would be a misleading cross-study join. The radiation generality failure remains unchanged.

### What camera time supports, separately

The publicly observed Serengeti/Ireland photo-time analyses show source/photo-time reallocations across matched photoperiod seasonal branches, but unverified camera operating hours and detection effects prevent direct promotion into latent organism time use. We do not pool camera photo hours and GPS night−day step speeds into the same temporal response or claim common species physiological clock without a comparable observation model.

## Discussion and prospective ecology

**Primary biological interpretation:** A light-calendar-associated component of relative day/night GPS displacement can reappear in distinct original field study cohorts without preserving a simple compass-relative distribution of displacement directions. The appropriate generalization target is axis- and hierarchy-specific, not an invariant "solar response" scalar.

**Next original-data analysis (new, source-frozen route):** Directly reconnect to original tagged-individual GPS events under exact provider checksums, observe same individual, source, comparable seasons and movement segments. Let apparent solar-time be τ, Sun altitude h, ERA5 shortwave R, temperature T and 100m wind vector W. Estimate a shortwave anomaly \(R^{*}=R-E[R|h,\tau,\mathrm{calendar},\mathrm{region}]\) with a prespecified flexible but bounded weather/sun model, preserving support. Analyze three responses on the *same GPS steps*: (a) observed step speed/daypart allocation, (b) \( \cos(B-A)\), and (c) \( \sin(B-A)\), with route and directional marginal controls. Aggregate within animal then original study; leave out full original field studies, not randomly held-out GPS fixes. Repeat in distinct movement guild/biome original datasets, keeping all source attrition. Treat R* as observational cloud/radiation contrast, not a causal exposure intervention, and ERA5 at 100m as not individual flight altitude wind.

No multilevel global population/clade coefficient is justified until source identity and geography overlap are audited and enough independent original population/biome samples cover time, direction and weather jointly.

## Source provenance, methods and record

- [Source-scoped 8 original research input blob hashes and freeze contract](../research/RSDS_CROSS_AXIS_GENERALITY_V0_CONTRACT.json)
- [Exact original study-taxon result and cross-scale matrix](../research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json)
- [Reproduction code](../scripts/analyze_rsds_cross_axis_generality_v0.py)
- [Known-truth and fail-closed tests](../tests/test_rsds_cross_axis_generality_v0.py)
- [Broader Japanese ecological program and constraints](../docs/RSDS_JOINT_SOLAR_TEMPORAL_MOVEMENT_MACROECOLOGY_20261026.md)

**Scientific boundaries:** This is exploratory after original results, not a pre-registered inter-axis independence test, does not establish causal photoperiod/solar orientation, and does not estimate what fraction of weather or resource state mediates behavior. Original animal coordinates, track/event timestamps and individual IDs are not redistributed.
