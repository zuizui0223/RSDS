"""Collect four separate source GPS joint covariance artifacts WITHOUT taxa pooling.

The only direct source-to-source matched species pair is mule deer-minus-
bighorn in Utah and Wyoming. Both are previously outcome-aware selected.
This code exposes the descriptive comparison but never assigns global
biological meta-analysis significance to four purposively sampled studies.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.run_rsds_original_joint_covariance_v1 import EXPECTED
from scripts.rsds_original_joint_calendar_covariance_v1 import (
    B,NEGATIVE,SOURCE_SHA,PARENT_MONTH,quadratic
)

STUDIES=("Illinois","Wyoming","Utah","Geese")
METHODS=("paired_month_only","paired_month_plus_animal_weights")
PREFIX="q_within_month_slope_per_hour:"

def aggregate(objects):
    if set(objects)!=set(STUDIES):
        raise ValueError("Four original DOI slots required; no source substitution")
    seen=set()
    digest=[]
    for study in STUDIES:
        obj=objects[study]
        expected=EXPECTED[study]
        if (obj.get("status")!="ORIGINAL_SOURCE_JOINT_COVARIANCE_ESTIMATED"
            or obj.get("original_DOI")!=expected["doi"]
            or obj.get("original_source_GPS_MD5_verified") is not True
            or obj.get("original_native_q_daylength_effect_replayed") is not True
            or obj.get("original_source_result_SHA")!=SOURCE_SHA[study][1]
            or obj.get("original_four_month_result_SHA")!=PARENT_MONTH[1]
            or obj.get("historical_external_shortwave_nonreplication_UNCHANGED")!=NEGATIVE
            or obj.get("historic_12_ERA5_raw_JSON_HOLD_UNCHANGED") is not True
            or obj.get("no_individual_date_location_or_weight_draws_exported") is not True
            or obj.get("original_source_shared_calendar_days")!=expected["days"]
            or len(obj.get("source_species_in_fixed_order",[]))!=expected["taxa"]):
            raise ValueError("Original source metadata/bio/replay/HOLD changed: "+study)
        if obj["original_DOI"] in seen:
            raise ValueError("Repeated source DOI as independent study")
        seen.add(obj["original_DOI"])
        species=obj["source_species_in_fixed_order"]
        p=2*len(species)
        for mode in METHODS:
            z=obj["covariance_sensitivity"][mode]
            if z["replicates"]!=B or len(z["drawn_parameter_order"])!=p:
                raise ValueError("Missing joint source covariance replication")
            if z["drawn_parameter_order"]!=(
                ["q_date_equal:"+s for s in species]
                +[PREFIX+s for s in species]):
                raise ValueError("Joint covariance dimension/order changed")
            V=z["joint_sampling_covariance_matrix"]
            if len(V)!=p or any(len(row)!=p for row in V):
                raise ValueError("Mismatched source covariance matrix")
            for i in range(p):
                if V[i][i]<0 or not math.isfinite(V[i][i]):
                    raise ValueError("Invalid source covariance diagonal")
                for j in range(p):
                    if not math.isclose(V[i][j],V[j][i],abs_tol=1e-12,rel_tol=0):
                        raise ValueError("Source native covariance not symmetric")
            for k in range(p):
                c=[0.0]*p
                c[k]=1;c[(k+1)%p]=-1
                if quadratic(V,c)<-1e-10:
                    raise ValueError("Not a PSD source covariance")
            if len(z["species_contrasts"])!=(3 if len(species)==3 else 1):
                raise ValueError("Source pair was dropped/reselected")
            if any(not v.get("source_month_slope_prior_replayed")
                   for v in z["species_contrasts"]):
                raise ValueError("Historical source month slope was reclassified")
        digest.append({
            "study":study,"doi":obj["original_DOI"],
            "native_taxa":species,
            "native_row_count":len(species),
            "calendar_days":obj["original_source_shared_calendar_days"],
            "eligible_month_days":obj["eligible_common_month_days"],
            "eligible_month_blocks":obj["eligible_common_year_month_blocks"],
            "source_effects":obj["date_equal_estimates"],
            "source_joint_covariance":obj["covariance_sensitivity"],
        })
    a=objects["Utah"]["covariance_sensitivity"]
    b=objects["Wyoming"]["covariance_sensitivity"]
    comparisons={}
    for mode in METHODS:
        ut=a[mode]["species_contrasts"]
        wy=b[mode]["species_contrasts"]
        if len(ut)!=1 or len(wy)!=1 or [
            ut[0]["first_taxon"],ut[0]["second_taxon"]] != [
            "Odocoileus hemionus","Ovis canadensis"] or [
            wy[0]["first_taxon"],wy[0]["second_taxon"]] != [
            "Odocoileus hemionus","Ovis canadensis"]:
            raise ValueError("Cannot compare different native ungulate pairs")
        q_ut=ut[0]["source_date_equal_q_difference"]
        q_wy=wy[0]["source_date_equal_q_difference"]
        sl_ut=ut[0]["source_month_slope_difference_per_hour"]
        sl_wy=wy[0]["source_month_slope_difference_per_hour"]
        # Two distinct original source DOIs are treated as conditionally
        # independent *for this descriptive variance proxy only*.
        comparisons[mode]={
            "source_pair":"Odocoileus hemionus MINUS Ovis canadensis",
            "Utah_q":q_ut,"Wyoming_q":q_wy,
            "Utah_minus_Wyoming_q":q_ut-q_wy,
            "conditional_source_independent_q_difference_se_proxy":
                math.hypot(ut[0]["joint_bootstrap_q_difference_se"],
                           wy[0]["joint_bootstrap_q_difference_se"]),
            "Utah_month_slope":sl_ut,"Wyoming_month_slope":sl_wy,
            "Utah_minus_Wyoming_month_slope":sl_ut-sl_wy,
            "conditional_source_independent_slope_difference_se_proxy":
                math.hypot(ut[0]["joint_bootstrap_month_slope_difference_se"],
                           wy[0]["joint_bootstrap_month_slope_difference_se"]),
            "not_a_two_study_global_meta_analysis":True,
            "post_outcome_selected_species_and_studies":True,
        }
    return {
        "version":"RSDS_FOUR_ORIGINAL_JOINT_MONTH_COVARIANCE_ATLAS_V1",
        "stage":"POST_OUTCOME_FOUR_SOURCE_GPS_DESCRIPTIVE_DEPENDENCE",
        "study_DOIs":4,"study_native_rows":10,
        "source_pair_rows":8,
        "sources":digest,
        "same_native_pair_two_DOI_descriptive_covariance_sensitivity":comparisons,
        "original_48_source_negative_shortwave_external_nonreplication_UNCHANGED":NEGATIVE,
        "original_12_ERA5_raw_JSON_HOLD_UNCHANGED":True,
        "no_raw_GPS_individuals_locations_dates_or_replicates_in_artifact":True,
        "no_world_population_effect_or_prospective_external_replication":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    data={}
    for study in STUDIES:
        path=args.directory/("joint_covariance_"+study+".json")
        data[study]=json.loads(path.read_text(encoding="utf-8"))
    result=aggregate(data)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    display={
        "studies":result["study_DOIs"],"taxa":result["study_native_rows"],
        "pair_rows":result["source_pair_rows"],
        "sources":[
            {"study":s["study"],"days":s["eligible_month_days"],
             "blocks":s["eligible_month_blocks"],
             "effects":s["source_effects"],
             "within_native_month_corr_q":(
                s["source_joint_covariance"]["paired_month_only"][
                    "joint_sampling_covariance_matrix"][0][1]
            )}
            for s in result["sources"]
        ],
        "mule_bighorn":result[
            "same_native_pair_two_DOI_descriptive_covariance_sensitivity"],
        "no_global_significance":True
    }
    print("RSDS_JOINT_COVARIANCE_ATLAS",json.dumps(display,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
