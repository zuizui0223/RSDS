# RSDS Wave3 original GPS ecology: source-specific night/day movement, not a global direction

**2026-11-15 | First frozen 51-source original-data result, followed by derivative-only source-native bootstrap CI audit.** The initial original GPS result is immutable at Git blob `4644b8ed38ae74dae345df8d90ea182f5dce4ee3` and was produced by original 16-source first-run [Action #43702456789](https://github.com/zuizui0223/RSDS/actions/runs/43702456789). The original first 51-source aggregate was copied from that run and frozen by [Action #43702720110](https://github.com/zuizui0223/RSDS/actions/runs/43702720110). The source CI audit was performed **only from first-run aggregate artifacts**, with no new GPS reads, via [Action #43703095918](https://github.com/zuizui0223/RSDS/actions/runs/43703095918).

## Original source-selection denominators

The fixed 2026-10-19 Movebank catalog had 437 published packages, 259 GPS/GNSS source candidates and 120 source-only metadata-priority candidates. Earlier Wave1+Wave2 selected 35 distinct original repository study IDs, yielding complete original registered-native q for 10. This Wave3 continued the **original fixed metadata SHA ranking** at single-native ranks 25–40 (16 new sources, quota chosen after earlier outcomes) and preserved all original source file holds.

Among 16 new selected source studies: **5 original-file HOLD** (three >200MB/unknown bytes, two ambiguous primary files), **1 additional original GPS/schema HOLD**, **4 parsed GPS sources with no q native support**, and **6 sources with original q and minimum independent tagged animals**. All 16 original source IDs stayed in the denominator; no sources/taxa were changed after signs were seen.

Combining 35 old and 16 new source units gives 51 source-study IDs, 66 registered source-native rows, 16 source-complete q estimates and 12 source-complete daylength estimates. **35/51 lack full q**, not zero q. Original human-location privacy restriction, original 48-study external solar-azimuth/signed irradiance **NOT REPLICATED**, and old 12 ERA5 source identity **HOLD** remain unchanged.

## Original source-native q point estimates AND tagged-animal bootstrap 95% intervals

`q` is the original source-native, animal-equal difference in **log1p median observed GPS straight-line step velocity** between NIGHT and DAY under the frozen observation protocol. This is *not time spent active*, fundamental chronotype, true solar irradiance selection or causal daylength response.

| Newly source-complete original species | q point | Original individual bootstrap 95% | Nominal source interval classification |
| --- | ---: | --- | --- |
| *Phoebastria irrorata* | −0.55209 | [−0.69092, −0.39621] | negative |
| *Odocoileus virginianus* | −0.00638 | [−0.01749, +0.00486] | includes zero |
| *Tapirus terrestris* | +0.04340 | [+0.03273, +0.05516] | positive |
| *Numenius arquata* | −0.04980 | [−0.08882, −0.01091] | negative |
| *Branta leucopsis* | −0.02611 | [−0.06660, +0.01469] | includes zero |
| *Tockus deckeni* | −0.15084 | [−0.18628, −0.11599] | negative |

Among these 6 originally supported sources: **five negative and one positive point estimate**, but **three nominal negative CI, one nominal positive CI, two CI including zero**. The original older 10 complete-source respondents had 7 positive and 3 negative point estimates. All **16** original full-source respondents now have **8 positive and 8 negative point signs**. No world-level 50:50 conclusion is licensed: source completeness and different species/source designs are informative, missing source effects are not zero, and the source-quota choice was adaptive.

The unweighted mean across six new q respondents was −0.12363; leave-one-original-source-out descriptive means ranged **−0.15704 to −0.03794**, still negative. Across all 16 q respondents the unweighted source mean was approximately **+0.00339**. These are respondent-only descriptions; not a 51-source, 120-frame, or global ecological estimand. They are not formally unbiased estimates and use **no** missingness correction. Old Wave1+Wave2's response-only poststratified q estimate +0.07568 is *a different weighting scheme* and must not be treated as the same as this new source-unweighted mean.

## An independent biological point: position and seasonal variation need not coincide

Within the six new q-supported source taxa, four also had separately supported original q–astronomical-daylength **slope** estimates:

| Species | Original q slope per h of daylength | Individual bootstrap 95% | q mean interval |
| --- | ---: | --- | --- |
| *Odocoileus virginianus* | −0.04300 | [−0.07409, −0.01251] | q mean interval includes zero |
| *Tapirus terrestris* | +0.01254 | [−0.00333, +0.02939] | q mean positive |
| *Numenius arquata* | +0.00595 | [−0.00944, +0.02185] | q mean negative |
| *Branta leucopsis* | +0.08000 | [+0.01391, +0.16333] | q mean interval includes zero |

**Two of four point slope signs differ from their overall q mean point sign**, but only two slopes exclude zero, and none of these comparisons estimates a cross-source, within-calendar source covariance. In particular, both *O. virginianus* and *B. leucopsis* have nominal daylength slopes excluding zero while their **mean q interval includes zero**. This supports treating **day-night displacement niche position** and **seasonal modulation** as distinct descriptive biological axes; it does not demonstrate causal responses to solar photoperiod, photoreceptors or competition.

## What is a genuine positive result

- **Source-level ecological heterogeneity:** Original GPS time-niche speed contrasts differ in their estimated signs, with several independent source bootstrap intervals on opposite sides of zero. A single across-taxon positive-night direction is not confirmed.
- **Different temporal axes:** source-wide NIGHT–DAY mean contrast and its original within-animal/year astronomical-daylength association are not interchangeable.
- **Ecologically consequential measurement selection:** original GPS data volume and admissible segment count do not guarantee source q, which requires repeated independent tagged animals across both dayparts (previous original 35-source observation audit). All missing sources are explicit.

**Not established:** a global species proportion with nocturnal movement, true activity time, solar-dependent navigation, taxonomic cause, within-species geographic portability, latitude/irradiance causation, or standard meta-analytic independent ecological populations.

## Original and audit result receipts

- [Source-only pre-event metadata rank manifest](../research/RSDS_437_WAVE3_16_SOURCE_METADATA_ONLY_FIRST_RESULT_20261115.json)
- [Pre-GPS original provider file MD5 receipt](../research/RSDS_437_WAVE3_16_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1.json)
- [Frozen first 35+16 source original GPS aggregate](../research/RSDS_437_WAVE3_PLUS_PREVIOUS_35_ORIGINAL_GPS_51_SOURCE_META_FIRST_RESULT_20261115.json)
- [First artifact source-level original tagged-animal bootstrap CI post-outcome audit](../research/RSDS_WAVE3_FIRST_ORIGINAL_SOURCE_Q_CI_POSTOUTCOME_RESULT_20261115.json)
- [Post-outcome source artifact digests and scientific contract](../research/RSDS_WAVE3_FIRST_ORIGINAL_SOURCE_Q_CI_AUDIT_V1_CONTRACT.json)

### Next ecological priority

To establish a transferable biological pattern rather than an across-study source mixture, fix **within-species across-study matched original-source pairs** and compare q and daylength slope on a shared biological calendar/observation design, separately from the 51-source ecology descriptive inventory. Original source-only metadata must adjudicate wild versus managed/domestic domains *without seeing q signs*. Any truly held-out new population effect requires a genuinely new pre-event source selection, not retroactive promotion of these post-outcome summary contrasts.
