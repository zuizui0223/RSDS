# RSDS — 同一個体の時間利用と太陽相対移動に対する日射・風の条件付き関係（原GPS v0）

**2026-10-26 / Original Utah & Wyoming publicly archived Movebank GPS / New post-prior-result exploratory empirical route.**

## 現在の問い

RSDSで前に確認した「ミュールジカの日長に対する夜昼相対変位は2調査で同方向だが、太陽相対移動方向の生平均は逆符号」という対照は、個体内の日射反応を直接測定していない。

そこで原公表GPSを再度読み出し、**同じ元タグ個体、同じ適格観測日、同じGPS原移動区間**から、時間ニッチ代理量と方向ニッチを再計量する。日射量の気象による残余変動がこれら2つの移動量とどう共変するのか、元個体内で条件付きに記述する。

元の38〜48研究の太陽方位／日射結果を上書きしない。新しい独立検証標本ではなく、以前の観測結果から課題を定式化した**別ルートの探索的原データ再解析**である。

## Source frames — no effect-based source choice

| Original study | Original Movebank public source DOI | Native tracked taxa | Prior qualified result to replay |
|---|---|---|---|
| Utah UDWR 2019–2020 | 10.5441/001/1.711 | *Antilocapra americana*, *Cervus elaphus*, *Odocoileus hemionus*, *Ovis canadensis* | Each species' native q NIGHT−DAY and individual/year q~astronomical daylength, original source animal counts and means |
| Wyoming 2019–2020 | 10.5441/001/1.703 | *Odocoileus hemionus*, *Ovis canadensis* | The original separate two-source q NIGHT−DAY and q~daylength, original counts and means |

We require exact original published source bundle UUID, size and provider-MD5; replay original full unselected native taxon results before new hourly ERA5 is requested. The original source q/night/day is observed *minimum straight-line GPS step displacement speed*, not inferred time spent truly active.

## Three coupled yet distinct original-individual outcomes

**(1) Within-source q on the same accepted animal observation dates.** After original one-fix-per-UTC-hour selection, an individual date qualifies only when ≥3 native source DAY and ≥3 NIGHT original 0.5–2-hour GPS intervals have valid speed (<60 km/h), and original animal/calendar year has ≥10 qualifying dates. Daily q = log1p(median nighttime GPS step speed) − log1p(median daytime GPS step speed).

**(2) Solar vector and shortwave on the same individual, date, originally qualified daylight GPS segments.** Use instantaneous Sun azimuth A, solar altitude h, apparent solar clock τ, travel bearing B (GPS displacement, not body axis); responses are cos(B−A) and sin(B−A). Require movement ≥0.1 km, sunlight altitude ≥10° at endpoints/midpoint and original interval 0.5–2 h. Link the same qualifying source gridyear to UTC-hourly explicitly ERA5 100m wind, surface shortwave and 2m air temperature, with time and units validated and only provider SHA256 released.

**(3) Conditional response with 1h solar clock / geometry / wind.** Select at most the top two source-original 1° grid×year frames per native taxon by *independent animals and original qualified q-day counts* before seeing new weather/solar effects; no replacement. Within animal/year/month×one-hour apparent-solar-clock blocks (≥6 original daylight GPS intervals on ≥3 solar dates), double-residualize R and the solar-heading responses on solar altitude and its square, sin/cos solar azimuth, daylength, 2m temperature, wind flow-to/Sun angle and wind speed. Calculate signed individual observational partial slopes per +100 W/m² modeled shortwave only when residual weather SD≥20 W/m² and enough source blocks/dates.

For the **same animal/dates** where ≥3 such geometrically qualified direction steps and original q are available, take mean event shortwave residual and relate daily q to the residual within animal/calendar-month, also accounting for daylength, mean 2m temperature and mean 100m wind. Require ≥20 original paired days, ≥2 calendar months and residual variation support. This is the tested "joint" feature; it does *not* mean night heading angles or radiative skin exposure were measured.

We report q weather coefficient and cos/sin heading weather coefficients for the **same qualified individual**; additionally report the signed Spearman association of the paired individual coefficients only where enough original tags support them.

## Ecological hierarchy and inference limit

The original source×taxon includes many tagged animals, but **only two regional field-study source IDs**. Multiple selected gridyears from one source share individuals and years. Therefore: aggregate within an original animal across selected contexts, then equal original animals within a taxon/source; bootstrap original animal IDs (never GPS steps). Original source weather analysis and causal radiant heating are different hypotheses; original brightness, wind at ground/100m, daylength, migration, vegetation phenology, micro-shade and species' movement strategies remain interdependent.

Source q days, solar-bearing daylight moving steps, the 1° source-gridyear, and fine solar-clock blocks impose eligibility constraints and can yield "UNSUPPORTED" despite strong original overall source q. We retain every original native taxon with insufficient-source-weather/geometry/temporal supports and never assume an unestimated effect is ecological zero.

### Original code and first outcomes

- [Frozen scientific contract](../research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_CONTRACT.json)
- [Original source GPS replay and same-tagg​ed-animal analysis code](../scripts/run_rsds_same_animal_solar_weather_v0.py)
- [Known-truth test cases and source-gates](../tests/test_rsds_same_animal_solar_weather_v0.py)
- [Native source first-results JSON](../research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json) — generated only after the first two source jobs finish; do not invent source weather effects beforehand.
- [Primary original time×direction study](../manuscript/RSDS_SOLAR_TIME_DIRECTION_MACROECOLOGY_v0_20261026.md)

**Do not interpret a positive or negative newly estimated weather slope as a solar compass, thermoregulation, competition, an experimentally manipulated weather effect, or a world-population parameter.** A new independent original field-source validation frame is required before that type of generality claim.
