"""Full frozen four new published multi-native Movebank study temporal ecology.

One DOI/study ID is ONE field study. All ten original registered native
study×taxon records and every pre-event source MD5/schematic/temporal HOLD
are preserved. A common original recorded CALENDAR date across species is
not experimental sympatry, the same experienced weather or an independent
individual-level interaction.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE
from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import (
    fixed_sources, NEG, SOURCE_KEY
)

OUTPUT=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_GPS_SHARED_CALENDAR_V1_FIRST_RESULT_20261029.json"


def aggregate(sources):
    spec,orig=fixed_sources()
    if len(sources)!=4 or set(sources)!=set(range(4)):
        raise ValueError("Original source study sampling frame not complete")
    rows=[];taxa=[];states=Counter();source_effects=Counter()
    for i,(s,e) in enumerate(zip(spec["sources"],orig)):
        x=sources[i]
        expected=spec["sources"][i]
        if (x.get("source_study_index")!=i or
                x.get("version")!=SOURCE_KEY or
                x.get("original_study_DOI")!=expected["doi"] or
                x.get("registered_native_taxa")!=expected["registered_taxa"] or
                x.get("previous_48_source_negative_directional_radiation_replication_UNCHANGED")!=NEG
                or not x.get("previous_12_ERA5_source_version_HOLDS_UNCHANGED")):
            raise ValueError("Unregistered source study/DOI/native taxa or negative historical result")
        if e["source_status"]=="ORIGINAL_SOURCE_METADATA_HOLD" and \
           x["status"]!="PRE_EVENT_SOURCE_METADATA_HOLD":
            raise ValueError("Pre-event source metadata HOLD was reclassified as biology")
        natives=x["original_native_taxon_results"]
        if len(natives)!=len(expected["registered_taxa"]) or \
           {z["taxon"] for z in natives}!=set(expected["registered_taxa"]):
            raise ValueError("Original source native taxa incomplete")
        status=x["status"]; states[status]+=1
        if status=="ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED":
            if not x["original_MD5_verified"] or not x["source_file_metadata_checked_again"]:
                raise ValueError("GPS effect opened before provider source original file MD5 verified")
            support=x["shared_source_calendar_support"]
            if not isinstance(support,dict) or support["registered_source_taxa"]!=len(natives):
                raise ValueError("Original shared-calendar source support missing")
            days=support["source_calendar_dates_shared_by_all_registered_taxa"]
            q_status={}
            for z in natives:
                if z["status"]!="ORIGINAL_NATIVE_Q_EVALUATED":
                    raise ValueError("Source original native q status not evaluated")
                for key,number in (("original_q_full_source_tag_effect","original_source_animals_with_diel_q"),
                                   ("original_q_daylength_full_source_tag_effect","original_source_animals_with_q_daylength")):
                    val=z[key]
                    if val is not None and val["n_individuals"]!=z[number]:
                        raise ValueError("Original source published q native support count impossible")
                if z["same_study_all_native_common_calendar_q"]["qualified_original_tagged_animals"]!=support["source_registered_native_q_tag_counts"][z["taxon"]]:
                    raise ValueError("Original native species q support differs from common calendar")
                for key in ("same_study_all_native_common_calendar_q",
                            "same_study_all_native_common_calendar_daylength_slope"):
                    k=z[key]
                    n=k["qualified_original_tagged_animals"]
                    if n<5:
                        if k["nominal_interval_sign"]!="UNAVAILABLE" or k["mean"] is not None or k["source_original_tag_bootstrap95"] is not None:
                            raise ValueError("Unsupported native source q/daylength effect became numeric")
                    else:
                        ci=k["source_original_tag_bootstrap95"]
                        if ci is None or len(ci)!=2 or not ci[0]<=k["mean"]<=ci[1]:
                            raise ValueError("Invalid source original-tagged animal bootstrap interval")
                q_status[z["taxon"]]=z["same_study_all_native_common_calendar_q"]["nominal_interval_sign"]
                source_effects["native_q_"+q_status[z["taxon"]]]+=1
                source_effects["native_daylength_"+z["same_study_all_native_common_calendar_daylength_slope"]["nominal_interval_sign"]]+=1
            pairs=x["all_native_study_pair_shared_calendar_results"]
            expected_pairs={tuple(z) for z in itertools.combinations(expected["registered_taxa"],2)}
            if len(pairs)!=len(expected_pairs) or {(p["native_first"],p["native_second"]) for p in pairs}!=expected_pairs:
                raise ValueError("Original fixed native source pair roster incomplete")
            for p in pairs:
                if p["same_source_calendar_dates_across_ALL_taxa"]!=days:
                    raise ValueError("Native study pair compared different source date populations")
                q=p["source_original_animal_q_pair"]
                a=next(z for z in natives if z["taxon"]==p["native_first"])
                b=next(z for z in natives if z["taxon"]==p["native_second"])
                if q["native_tag_support"]!=[
                        a["same_study_all_native_common_calendar_q"]["qualified_original_tagged_animals"],
                        b["same_study_all_native_common_calendar_q"]["qualified_original_tagged_animals"]]:
                    raise ValueError("Different tagged individual support in native q pair")
                if min(q["native_tag_support"])<5 and q["first_minus_second_q"] is not None:
                    raise ValueError("Unsupported pair was converted to numeric ecological difference")
                source_effects["pair_q_"+q["nominal_interval_sign"]]+=1
                season=p["source_common_day_q_difference_daylength"]
                if days and season is None:
                    raise ValueError("Source has common dates but missing pair photoperiod")
                if season:
                    if season["source_shared_calendar_days"]!=days:
                        raise ValueError("Original source pair calendar daylength support disagrees")
                    source_effects["pair_within_month_"+season["within_month_reference_daylength"]["nominal_interval_sign"]]+=1
        else:
            if status not in ("PRE_EVENT_SOURCE_METADATA_HOLD","ORIGINAL_SOURCE_OR_NATIVE_SHARED_CALENDAR_HOLD","SOURCE_NO_VALID_GPS_FIXES"):
                raise ValueError("Source outcome state not in frozen set")
        for z in natives:
            taxa.append({"study":expected["doi"],"native_taxon":z["taxon"],"source_status":status,
                         "native_original_q_tag_count":z.get("original_source_animals_with_diel_q"),
                         "shared_source_q":z.get("same_study_all_native_common_calendar_q"),
                         "shared_source_q_daylength":z.get("same_study_all_native_common_calendar_daylength_slope")})
        rows.append({
            "source_study_rank":i,
            "published_DOI":expected["doi"],
            "original_study_id":expected["study_id"],
            "guild_metadata_hint":expected["guild_hint"],
            "source_registered_native_taxa":expected["registered_taxa"],
            "source_status":status,
            "source_hold_error_type":x.get("source_hold_error_type"),
            "source_provider_original_MD5_replayed":x["original_MD5_verified"],
            "original_native_source_q_and_daylength":natives,
            "original_all_native_common_calendar_support":x["shared_source_calendar_support"],
            "all_source_native_pairs_original_q_and_season":x["all_native_study_pair_shared_calendar_results"]
        })
    if len(rows)!=4 or len(taxa)!=10:
        raise ValueError("Four original multi-native studies/10 native cohort rows lost")
    return {
        "version":"RSDS_FOUR_NEW_MULTI_NATIVE_STUDY_Q_SEASON_SHARED_CALENDAR_V1_FIRST_RESULT",
        "date":"2026-10-29",
        "stage":"SOURCE_ORIGINAL_EVENTS_OPENED_ONLY_AFTER_FOUR_MD5_SNAPSHOTS_FIXED",
        "fixed_original_published_studies":4,
        "registered_original_study_native_taxon_rows":10,
        "study_status_counts":dict(states),
        "supported_native_and_source_pair_nominal_sign_counts":dict(source_effects),
        "all_frozen_source_study_original_results":rows,
        "all_source_study_native_q_support":taxa,
        "prior_12_ERA5_version_source_HOLDS_UNCHANGED":True,
        "prior_48_source_shortwave_direction_NEGATIVE_EXTERNAL_REPLICATION_UNCHANGED":NEG,
        "scientific_claim_boundary":[
            "Original sources and native species were fixed on SOURCE metadata, before opening original animal temporal q data. The four published sources are purposive, not a representative world random sample.",
            "The geese source was earlier inspected for solar-heading observations, so this is new only for its same-calendar temporal ecology route, not untouched independent external all-phenotype validation.",
            "Source common astronomical calendar dates across registered taxa do NOT prove their spatial co-occurrence or actual individual encounters, common microclimate or biological competition.",
            "Study-level reference astronomical daylength derives from median native source GPS latitudes, not direct per-animal experienced sunlight. Equatorial island sources may have inadequate daylight range and valid q without estimable slope.",
            "GPS q is log1p(NIGHT minimum straight-line recorded step speed) minus log1p(DAY counterpart), not actual activity-hour budget. Instrument interval and individual tags/day support can distort rankings.",
            "Source native taxon q and photoperiod slopes may use distinct qualified tag sets. No unsupported effect or climate slope converted to zero.",
            "Month-block nominal 95% uncertainty is post-selection exploratory, cannot identify causal photoperiod from date, weather, breeding and seasonal behavior.",
            "Each original study DOI is one replication unit; multiple native taxa and pairwise comparisons within it are dependent, no P-value from independent GPS fixes.",
            "48-study shortwave radiation directional external negative replication and 12 old ERA5 serialized JSON weather version HOLDs remain unchanged.",
            "No original animal tags, private GPS coordinates/time, source original individual-year, calendar dates or individual coefficients publicly redistributed."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input",type=Path,required=True)
    cli.add_argument("--out",type=Path,default=OUTPUT)
    args=cli.parse_args()
    input={}
    for i in range(4):
        f=args.input/("study_"+str(i)+"_original_multinative_q_v1.json")
        if not f.is_file():
            raise ValueError("Missing one fixed original field study; no source replacement: "+str(i))
        input[i]=json.loads(f.read_text(encoding="utf-8"))
    out=aggregate(input)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_FOUR_MULTI_NATIVE_ORIGINAL_TEMPORAL_ATLAS",json.dumps({
        "studies":out["fixed_original_published_studies"],
        "study_native_taxa":out["registered_original_study_native_taxon_rows"],
        "source_states":out["study_status_counts"],
        "sign_status":out["supported_native_and_source_pair_nominal_sign_counts"],
        "prior_ext_negative":out["prior_48_source_shortwave_direction_NEGATIVE_EXTERNAL_REPLICATION_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
