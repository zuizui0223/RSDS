# Time reallocation and solar-relative direction in repeatedly observed ungulates

**RSDS original-GPS manuscript v2 — 26 October 2026. Exploratory, noncausal.**

## Abstract

Whether animal time use and travel direction have comparably repeatable relationships with the Sun is a central question in macroecology. We separated observed day–night GPS displacement from travel bearing relative to instantaneous Sun azimuth, using two published Movebank original GPS studies: Utah (four ungulate taxa; DOI 10.5441/001/1.711) and Wyoming (two taxa; DOI 10.5441/001/1.703). Provider MD5 and previous study-native GPS daylength/velocity effects were replayed before new source directional results were computed. Within each native taxon we selected two 1-degree grid×year observational contexts by maximum overlap of q-qualified original tagged animals, without using response magnitudes or signs.

All six source×taxon groups supported within-tag comparisons of the day–night minimum GPS step displacement difference and solar-relative travel components, totaling **269 paired-tag occurrences** across six source-group comparisons (not 269 independent populations). Descriptive original-animal bootstrap 95% intervals excluded zero for the night–day q contrast in five of six groups (two positive, three negative), but only two of six groups for raw \(\cos(B-A)\) movement, and zero of six for changes in the **joint circular covariance** of geographical heading and solar azimuth. In the two groups with supported changes in raw solar-relative cosine, the geographical bearing/Sun marginal term also changed. Sun-perpendicular sine contrasts varied in sign among three groups with intervals excluding zero. The repeated bighorn sources both showed negative q contrasts, but context ranks refer to sampling support, not any common increase in solar exposure or chronological season.

These observations support recorded within-individual heterogeneity in the **temporal displacement proxy** and show why apparent Sun-relative direction is not interchangeable with an actual solar-orientation response. They neither identify causal photoperiod or solar radiation responses nor provide globally representative ecological generality. Historical ERA5 source-version holds (12/12) and the older externally nonreplicated radiation direction hypothesis remain intact.

## Original-data methods

The registered first source results used a timestamp-only selection of one original fix per UTC hour and 0.5–2-hour accepted GPS intervals, with source-native daypart and displacement criteria and the previous q and daylength means exactly replayed. An observation day contributed q only when at least three accepted DAY and three NIGHT intervals were present in an original animal-year with at least ten such days. We defined \(q=\log(1+\widetilde v_{\rm night})-\log(1+\widetilde v_{\rm day})\), where v is observed straight-line displacement speed (km/h), **not latent activity duration**.

Per source×native taxon, two original source 1-degree grid×calendar-year contexts were selected to maximize actual repeated tag identities with sufficient predeclared q days. Context rank2 minus rank1 reflects recorded GPS support, **not high-versus-low radiation, spring-versus-autumn, a climatic treatment, or a shared between-study gradient**. For the same original animals and q-supported days, daylight directional steps required movement ≥0.1 km, Sun altitude ≥10° and source-valid step intervals. Direction was expressed as \(Y_{\parallel}=\cos(B-A)\) and \(Y_{\perp}=\sin(B-A)\), with GPS bearing B and computed instantaneous Sun azimuth A.

The identity
\[
E[\cos(B-A)]=E[\cos B]E[\cos A]+E[\sin B]E[\sin A]
+\operatorname{Cov}(\cos B,\cos A)+\operatorname{Cov}(\sin B,\sin A)
\]
decomposes a raw Sun-relative mean into geographical heading–Sun **marginal composition** and a **joint circular covariance**. Original tagged animals were the paired resampling units (2,000 draws; nominal descriptive 95% intervals). No GPS fixes were treated as independent population replications; two source studies remain only two ecological field studies.

## Results

| Original study/taxon | Paired tags | Δq | Δraw cos(B−A) | Δgeographic marginal | Δjoint circular |
|---|---:|---:|---:|---:|---:|
| Utah / *Antilocapra americana* | 33 | +0.01748 | +0.02531 | +0.02113 | +0.00418 |
| Utah / *Cervus elaphus* | 50 | +0.01318 | −0.00177 | −0.01117 | +0.00940 |
| Utah / *Odocoileus hemionus* | 101 | −0.01111 | +0.02374 | +0.00788 | +0.01587 |
| Utah / *Ovis canadensis* | 27 | −0.02979 | +0.05850 | +0.03133 | +0.02717 |
| Wyoming / *Odocoileus hemionus* | 37 | +0.00329 | +0.00337 | +0.00211 | +0.00125 |
| Wyoming / *Ovis canadensis* | 21 | −0.01120 | −0.01522 | +0.00706 | −0.02228 |

For q, source-tag bootstrap intervals excluded zero in five groups (all except Wyoming mule deer). Utah bighorn q delta was −0.02979, CI [−0.04483, −0.01631], while Wyoming bighorn was −0.01120, CI [−0.01664, −0.00549]. These contrasts represent different source coverage-ranked contexts and are **not replicated contrasts along a controlled environment axis**.

Raw solar-relative cosine intervals excluded zero in Utah pronghorn and bighorn (both positive). In those same two groups the marginal term intervals also excluded zero. All six joint-circular delta intervals included zero. Perpendicular sine differences were positive for Utah pronghorn and mule deer, negative for Utah elk (three nominal interval signs), again contradicting a single common signed heading rule. All source-axis intervals are unadjusted for multiple comparisons and conditional on the sampled original tags.

## Interpretation and limitations

These original data advance from comparing different animals' source averages to **matching the same original tag between two observed contexts**. Differences in q are apparent within individuals but may reflect different sampling calendars, geographical routes, collar duty cycles, daylight availability, or other unmeasured factors. The time and bearing endpoints share individual and qualifying date, but not necessarily identical GPS intervals: q pools day/night displacement and can include valid segments outside the specific grid used to select solar-bearing daylight steps. Thus the result is not a controlled within-habitat environmental experiment.

The raw cosine identity shows precisely where directional geographical marginals enter. A zero-spanning joint interval **does not establish a zero covariance change** or prove that the marginal term causally explains the raw shift. Biological solar compass, solar radiation avoidance, thermoregulation, clock entrainment, and a difference in globally transferable temporal versus spatial niche organization remain unproven.

The prior ERA5/Open-Meteo source replication audit produced 12/12 serialized weather JSON hash mismatches and no valid paired radiation effect. This GPS-only route does **not** silently substitute re-fetched weather. A separate source-frozen comparison of numerical UTC-hourly values, units and astronomical geometry (rather than serialized metadata), followed by prospectively specified field-study validation, is necessary to test whether actual shortwave anomalies contribute to time and direction jointly.

## Evidence and reproducibility

The [original first-result matrix](../research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT_20261026.json), [frozen contract](../research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_CONTRACT.json), [first-result interpretation](../research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_ECOLOGY_SYNTHESIS_20261026.json), [original-source analysis](../scripts/run_rsds_same_tag_gps_geometry_v2.py), [source matching tests](../tests/test_rsds_same_tag_gps_geometry_v2.py) and [Japanese ecological synthesis](../docs/RSDS_SAME_TAG_GPS_GEOMETRY_ECOLOGY_20261026.md) record complete source attribution and inferential limits.

No original tagged animal IDs, source GPS trajectories, exact event coordinates or timestamps, individual response values or restricted original source files have been distributed.
