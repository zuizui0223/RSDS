"""Freeze ALL untouched original public GPS source studies for six prior q-repeat taxa.

This selection is outcome-guided at the TAXON level (six already discovered taxa),
not at SOURCE level. No NEW provider file listing or animal/GPS row is queried.
All six previously repeated taxa are kept, including one with zero new candidates.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import (
    HERE,pinned,study_id_versions,standard_date
)

CONTRACT=HERE/"research/RSDS_SIX_TAXON_UNTOUCHED_SOURCE_SELECTION_CONTRACT_20261115.json"
VERSION="RSDS_SIX_PREVIOUS_NATIVE_REPEATS_ALL_UNTOUCHED_SOURCE_METADATA_V1"


def normalize_doi(v):
    if not isinstance(v,str):
        raise ValueError("Missing DOI")
    v=v.strip().lower()
    for prefix in ("https://doi.org/","http://doi.org/","doi:"):
        if v.startswith(prefix):
            return v[len(prefix):]
    return v


def fixed_input():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    expect={
        "version":VERSION,
        "expected_new_source_identity_slots":6,
        "expected_unseen_catalog_match_packages_regardless_eligibility":16,
        "previously_examined_source_studies":63,
        "previously_source_native_registered_rows":88,
        "previously_source_native_q_supported":38,
        "max_original_provider_file_bytes":200_000_000,
        "new_provider_file_metadata_opened":False,
        "new_original_animal_GPS_rows_opened":0,
        "control_legacy_48_signed_irradiance_external_nonreplication":"NOT_REPLICATED",
        "control_12_legacy_ERA5_raw_SHA":"HOLD",
    }
    if any(c.get(k)!=v for k,v in expect.items()):
        raise ValueError("Frozen source-only science/source-count contract changed")
    old=pinned(c["earlier_source_native_atlas_path"],
               c["earlier_source_native_atlas_git_blob"])
    if old["distinct_original_archived_source_DOIs_and_study_IDs"]!=63 or old["all_registered_original_source_native_rows"]!=88 or old["all_original_source_native_q_supported_rows"]!=38:
        raise ValueError("Older 63-source original-native result changed")
    six=sorted(r["scientific_name"] for r in old["same_species_source_study_q_comparison_candidates"])
    if six!=sorted(c["fixed_six_species_from_previous_result"]):
        raise ValueError("Previously observed six species labels changed")
    catalog=pinned(c["frozen_catalog_path"],c["frozen_catalog_git_blob"])
    if catalog.get("package_count")!=437 or len(catalog["packages"])!=437:
        raise ValueError("Original Movebank 437 catalog changed")
    return c,catalog,old


def old_source_ids(c):
    w1=pinned(c["wave1_original_source_selection_path"],c["wave1_original_source_selection_git_blob"])
    w2=pinned(c["wave2_original_source_selection_path"],c["wave2_original_source_selection_git_blob"])
    w3=pinned(c["wave3_original_source_selection_path"],c["wave3_original_source_selection_git_blob"])
    l=pinned(c["legacy_original_12_native_path"],c["legacy_original_12_native_git_blob"])
    old=[]
    old.extend((s["source_DOI_label"],s["original_repository_study_id"])
               for s in w1["pre_event_selected_sources_in_rank_order"])
    old.extend((s["original_DOI_labels"][0],s["original_repository_study_id"])
               for s in w2["all_original_wave2_sources_in_rank_order"])
    old.extend((s["original_DOI_labels"][0],s["original_repository_study_id"])
               for s in w3["all_wave3_selected_source_metadata_in_SHA_rank_order"])
    by={}
    for s in l["all_original_study_native_taxon_rows_including_HOLD"]:
        by.setdefault(normalize_doi(s["original_study_doi"]),s["original_study_doi"])
    if len(old)!=51 or len(by)!=12:
        raise ValueError("Previous 51 + legacy 12 source histories incomplete")
    old.extend((v,None) for v in by.values())
    if len({normalize_doi(d) for d,_ in old})!=63:
        raise ValueError("Repeated original DOI in prior 63 source results")
    return old


def select():
    c,catalog,old_native=fixed_input()
    fixed=old_source_ids(c)
    old_doi={normalize_doi(d) for d,s in fixed}
    older_study_ids={str(i) for _,i in fixed if i is not None}
    by_doi={}
    for p in catalog["packages"]:
        for doi in p.get("dois") or []:
            norm=normalize_doi(doi)
            if norm in by_doi:
                raise ValueError("Ambiguous catalog DOI")
            by_doi[norm]=p
    for doi in old_doi:
        if doi not in by_doi:
            raise ValueError("Previously original study DOI missing from fixed catalog")
        for x in by_doi[doi].get("repository_study_ids") or []:
            older_study_ids.add(str(x))
    if len(older_study_ids)!=63:
        raise ValueError("63 original study IDs are not independent keys")
    taxa=sorted(c["fixed_six_species_from_previous_result"])
    wanted=set(taxa)
    versions=study_id_versions(catalog["packages"])
    match=[p for p in catalog["packages"] if wanted.intersection(p.get("scientific_names") or [])
           and not any(normalize_doi(d) in old_doi for d in p.get("dois") or [])
           and not any(str(i) in older_study_ids for i in p.get("repository_study_ids") or [])]
    if len(match)!=16:
        raise ValueError("All unexamined six-taxon catalog identity opportunity changed")
    rows=[]
    for p in sorted(match,key=lambda p:(normalize_doi(p["dois"][0]) if p.get("dois") else "",p["package_uuid"])):
        names=p.get("scientific_names") or []
        idlist=p.get("repository_study_ids") or []
        dois=p.get("dois") or []
        selected=sorted(wanted.intersection(names))
        issues=[]
        if not set(p.get("sensor_names") or []).intersection({"GPS","GNSS"}):
            issues.append("NOT_ORIGINAL_GPS_GNSS_SENSOR")
        if p.get("licenses")!=["CC0 1.0 Universal"]:
            issues.append("NON_CC0_SINGLE_LICENSE")
        if type(p.get("reported_animal_count")) is not int or p["reported_animal_count"]<10:
            issues.append("INSUFFICIENT_METADATA_ANIMAL_COUNT")
        if type(p.get("reported_location_count")) is not int or p["reported_location_count"]<10000:
            issues.append("INSUFFICIENT_METADATA_FIX_COUNT")
        if len(idlist)!=1 or len(dois)!=1:
            issues.append("AMBIGUOUS_OR_MISSING_STUDY_ID_DOI")
        if not names or len(set(names))!=len(names) or any(not n.strip() for n in names):
            issues.append("AMBIGUOUS_NATIVE_LABELS")
        if len(idlist)==1 and idlist[0] in versions and versions[idlist[0]][0]["package_uuid"]!=p["package_uuid"]:
            issues.append("NON_CANONICAL_STUDY_ID_PACKAGE_VERSION")
        if len(idlist)==1 and str(idlist[0]) in older_study_ids:
            issues.append("ORIGINAL_STUDY_ID_PREVIOUSLY_EVALUATED")
        row={
            "original_source_DOI":dois[0] if len(dois)==1 else None,
            "original_source_study_id":idlist[0] if len(idlist)==1 else None,
            "original_package_UUID":p["package_uuid"],
            "original_public_title":p["title"],
            "original_native_taxa":names,
            "matched_prior_six_taxa":selected,
            "original_sensor_labels":p.get("sensor_names") or [],
            "original_license_labels":p.get("licenses") or [],
            "reported_animal_count_metadata_only":p.get("reported_animal_count"),
            "reported_location_count_metadata_only":p.get("reported_location_count"),
            "source_status":"FROZEN_NEW_SOURCE_ID_FOR_FILE_PREFLIGHT" if not issues else "METADATA_EXCLUDED_NO_SOURCE_ACCESS",
            "metadata_exclusion_reasons":issues,
            "provider_original_file_metadata_queried":False,
            "original_GPS_animal_events_opened":0,
        }
        rows.append(row)
    chosen=[x for x in rows if x["source_status"]=="FROZEN_NEW_SOURCE_ID_FOR_FILE_PREFLIGHT"]
    if len(chosen)!=6 or len({x["original_source_study_id"] for x in chosen})!=6 or len({normalize_doi(x["original_source_DOI"]) for x in chosen})!=6:
        raise ValueError("Six new source IDs or unique study DOI count changed")
    tally=Counter(s for row in chosen for s in row["matched_prior_six_taxa"])
    expected=c["expected_new_source_studies_by_taxon"]
    if any(tally.get(s,0)!=expected[s] for s in taxa):
        raise ValueError("Previously unseen eligible sources per six taxa changed")
    chosen.sort(key=lambda x:normalize_doi(x["original_source_DOI"]))
    for ix,x in enumerate(chosen):
        x["source_index"]=ix
    return {
        "version":VERSION,
        "stage":"FIRST_SOURCE_ONLY_SIX_TAXON_ORIGINAL_GPS_STUDY_IDENTITIES_BEFORE_NEW_FILE_ACCESS",
        "catalog_snapshot_date":"2026-10-19",
        "catalog_original_packages":437,
        "previous_study_DOI_count":63,
        "fixed_six_previous_native_q_repeat_species":taxa,
        "new_source_candidates_with_same_taxa_even_if_metadata_ineligible":len(rows),
        "new_original_metadata_qualified_source_studies":len(chosen),
        "remaining_metadata_unqualified_or_inapplicable_candidates":len(rows)-len(chosen),
        "new_qualified_source_counts_by_six_species":{s:tally.get(s,0) for s in taxa},
        "new_original_source_metadata_selection_rows":chosen,
        "all_six_taxa_unexamined_metadata_candidate_decisions_including_exclusions":rows,
        "all_new_source_selection_based_only_on_public_catalog_metadata":True,
        "taxon_selection_was_after_prior_biological_outcomes":True,
        "no_source_was_selected_based_on_new_animal_GPS_or_q_sign":True,
        "new_GPS_rows_opened":0,
        "new_provider_file_metadata_opened":False,
        "old_48_source_sunlight_direction_heldout_not_replicated":True,
        "old_12_ERA5_raw_SHA_hold":True,
        "source_population_independence_not_verified":True,
        "claim_ceiling":c["interpretation"],
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    d=select()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_SOURCE_ONLY_REPEAT_TAXA",json.dumps({
        "prior_sources":63,
        "target_taxa":6,
        "source_slots":len(d["new_original_source_metadata_selection_rows"]),
        "hold_or_excluded":d["remaining_metadata_unqualified_or_inapplicable_candidates"],
        "source_studies":[{"species":x["matched_prior_six_taxa"],"doi":x["original_source_DOI"]} for x in d["new_original_source_metadata_selection_rows"]],
        "no_GPS_or_original_file_metadata":True
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
