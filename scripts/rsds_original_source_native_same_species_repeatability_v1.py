"""Cross-study SAME-SPECIES original source-native GPS q audit, without new GPS.

Read only old 12 original source-native effects, 51 selected new source IDs
and previously SHA-frozen original source aggregate artifacts. Critically,
a multi-species SOURCE mean q is NEVER assigned to one named species.
All repeated taxa are identified from the ENTIRE archived source roster,
including unavailable Q, before looking at point signs.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from statistics import mean

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned

CONTRACT = HERE / "research/RSDS_ORIGINAL_SOURCE_NATIVE_CROSS_STUDY_REPEATABILITY_AUDIT_V1_CONTRACT.json"
VERSION = "RSDS_ORIGINAL_SOURCE_NATIVE_CROSS_STUDY_REPEATABILITY_AUDIT_V1"
NEG = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def norm_doi(s):
    if not isinstance(s, str) or not s.strip():
        raise ValueError("Missing source DOI")
    x = s.strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if x.startswith(prefix):
            return x[len(prefix):]
    return x


def q_effect(effect, old=False):
    """A source-native (NOT source-equal-all-species) original tag bootstrap."""
    if effect is None or effect.get("mean") is None:
        return {
            "available": False,
            "q": None,
            "original_source_tagged_animals": None,
            "source_tag_bootstrap95": None,
            "point_sign": "UNAVAILABLE",
            "nominal_source_CI_sign": "UNAVAILABLE",
        }
    v = effect["mean"]
    ci = effect.get("tag_bootstrap95") if old else effect.get("animal_bootstrap_95_percentile")
    n = effect.get("source_tags") if old else effect.get("n_individuals")
    if (not isinstance(n, int) or n < 5
        or not isinstance(v, (int, float)) or not math.isfinite(v)
        or not isinstance(ci, list) or len(ci) != 2
        or not all(isinstance(z, (int, float)) and math.isfinite(z) for z in ci)
        or ci[0] > ci[1]):
        raise ValueError("Invalid original source-native bootstrap support")
    sign = "POSITIVE" if v > 0 else "NEGATIVE" if v < 0 else "ZERO"
    support = "POSITIVE" if ci[0] > 0 else "NEGATIVE" if ci[1] < 0 else "INCLUDES_ZERO"
    return {
        "available": True,
        "q": float(v),
        "original_source_tagged_animals": n,
        "source_tag_bootstrap95": list(map(float, ci)),
        "point_sign": sign,
        "nominal_source_CI_sign": support,
    }


def from_wave3_ci(obj):
    if obj.get("available") is False:
        return q_effect(None)
    if obj.get("available") is not True:
        raise ValueError("Source Wave3 original first CI availability invalid")
    return q_effect({
        "mean": obj["mean"],
        "n_individuals": obj["n_individuals"],
        "animal_bootstrap_95_percentile": obj["bootstrap_95_percentile"],
    })


def make_row(source, taxon, q, daylight, study_complete=None):
    doi = norm_doi(source["doi"])
    if not taxon or not isinstance(taxon, str):
        raise ValueError("Blank source-native scientific name")
    return {
        "original_DOI": doi,
        "original_repository_study_id": str(source["study_id"]),
        "scientific_name_original_label": taxon,
        "source_wave": source["wave"],
        "registered_native_count": source["n_native"],
        "source_complete_all_native_q": study_complete,
        "original_native_q": q,
        "original_native_daylength": daylight,
    }


def build(source_directory):
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c.get("version") != VERSION
        or c.get("original_wave2_first_GPS_action_run_id") != 43692435152
        or c.get("original_wave1_all_native_effects_git_blob") != "57029f30b65991118325c7d0cec1d8b3a36f4427"
        or c.get("original_wave3_first_CI_audit_git_blob") != "57d8acf1f573e7573527b7cac9d21a76958a666c"
        or c.get("original_first_51_study_meta_git_blob") != "4644b8ed38ae74dae345df8d90ea182f5dce4ee3"
        or c.get("old_48_source_solar_irradiance_direction_replicated") is not False
        or c.get("old_12_ERA5_provider_raw_SHA_HOLD_unchanged") is not True
        or c.get("no_new_animal_GPS_or_provider_file_reads") is not True):
        raise ValueError("Original source scientific contract/source route changed")
    catalog = pinned(c["original_2026_10_19_movebank_catalog_path"],
                     c["original_2026_10_19_movebank_catalog_git_blob"])
    legacy = pinned(c["prior_12_original_source_native_studies_path"],
                    c["prior_12_original_source_native_studies_git_blob"])
    wave1 = pinned(c["original_wave1_all_native_effects_path"],
                   c["original_wave1_all_native_effects_git_blob"])
    w1 = pinned(c["original_wave1_selection_path"],
                c["original_wave1_selection_git_blob"])
    w2 = pinned(c["original_wave2_selection_path"],
                c["original_wave2_selection_git_blob"])
    w3 = pinned(c["original_wave3_selection_path"],
                c["original_wave3_selection_git_blob"])
    ci3 = pinned(c["original_wave3_first_CI_audit_path"],
                 c["original_wave3_first_CI_audit_git_blob"])
    meta = pinned(c["original_first_51_study_meta_path"],
                  c["original_first_51_study_meta_git_blob"])
    if (catalog.get("package_count") != 437 or
        len(catalog["packages"]) != 437 or
        len(legacy["all_original_study_native_taxon_rows_including_HOLD"]) != 22 or
        legacy["all_original_studies"] != 12 or
        len(wave1["all_first_wave_registered_native_group_effects"]) != 27 or
        len(w1["pre_event_selected_sources_in_rank_order"]) != 16 or
        len(w2["all_original_wave2_sources_in_rank_order"]) != 19 or
        len(w3["all_wave3_selected_source_metadata_in_SHA_rank_order"]) != 16 or
        len(ci3["new_16_source_original_native_ci_and_support_rows"]) != 16 or
        meta["original_source_study_units"] != 51 or
        len(meta["all_51_source_study_rows"]) != 51 or
        meta["original_source_study_unit_all_registered_native_rows"] != 66):
        raise ValueError("Original source study and registered native denominators changed")

    by_doi={}
    for p in catalog["packages"]:
        for item in p.get("dois") or []:
            doi=norm_doi(item)
            if doi in by_doi:
                raise ValueError("Original public DOI repeated across package records")
            by_doi[doi]=p
    canonical={}
    def register(wave,doi,study_id,taxa):
        k=norm_doi(doi)
        if k in canonical:
            raise ValueError("Source DOI repeated between historical and SHA-ranked waves")
        p=by_doi.get(k)
        if p is None:
            raise ValueError("Source DOI not found in fixed original 437 public catalog")
        if str(study_id) not in [str(x) for x in p["repository_study_ids"]]:
            raise ValueError("Selected original study-ID not equal catalog DOI")
        if set(taxa) != set(p["scientific_names"]):
            raise ValueError("Original public source native roster changed")
        canonical[k]={
            "doi": k, "study_id": str(study_id), "n_native":len(taxa),
            "wave":wave, "taxa":list(taxa),
            "original_catalog_title":p["title"],
            "original_catalog_package_uuid":p["package_uuid"],
        }
    for s in w1["pre_event_selected_sources_in_rank_order"]:
        register("wave1",s["source_DOI_label"],s["original_repository_study_id"],s["scientific_names"])
    for s in w2["all_original_wave2_sources_in_rank_order"]:
        register("wave2",s["original_DOI_labels"][0],s["original_repository_study_id"],s["scientific_name_labels"])
    for s in w3["all_wave3_selected_source_metadata_in_SHA_rank_order"]:
        register("wave3",s["original_DOI_labels"][0],s["original_repository_study_id"],s["scientific_name_labels"])
    prior = defaultdict(list)
    for z in legacy["all_original_study_native_taxon_rows_including_HOLD"]:
        prior[norm_doi(z["original_study_doi"])].append(z)
    if len(prior) != 12:
        raise ValueError("Original prior source DOI count changed")
    for doi, arr in prior.items():
        p=by_doi.get(doi)
        if p is None or len(p["repository_study_ids"]) != 1:
            raise ValueError("Original prior DOI/study ID metadata missing")
        register("legacy12",doi,p["repository_study_ids"][0],
                 [x["native_taxon"] for x in arr])
    if len(canonical) != 63 or len({x["study_id"] for x in canonical.values()}) != 63:
        raise ValueError("Original 12+51 study DOI / study-ID distinct accounting failed")

    rows=[]
    for doi,source in canonical.items():
        if source["wave"] == "legacy12":
            for original in prior[doi]:
                rows.append(make_row(
                    source, original["native_taxon"], q_effect(original["q"], old=True),
                    q_effect(original["daylength"], old=True), None
                ))
    wave1_by_doi=defaultdict(list)
    for original in wave1["all_first_wave_registered_native_group_effects"]:
        wave1_by_doi[norm_doi(original["source_study_DOI"])].append(original)
    if len(wave1_by_doi)!=16:
        raise ValueError("Wave1 original native source ID count changed")
    for s in canonical.values():
        if s["wave"]!="wave1":
            continue
        arr=wave1_by_doi[s["doi"]]
        if len(arr)!=s["n_native"] or set(x["native_taxon"] for x in arr)!=set(s["taxa"]):
            raise ValueError("Wave1 original per-native data missing")
        mi=next((x for x in meta["all_51_source_study_rows"][:16]
                 if norm_doi(x["source_DOI_labels"][0])==s["doi"]),None)
        if mi is None:
            raise ValueError("Original 35-source study ID not mapped to native effects")
        for original in arr:
            rows.append(make_row(s,original["native_taxon"],
                                 q_effect(original["q"]),q_effect(original["daylength"]),
                                 mi["complete_source_native_q"]))

    artifact_rows=c["original_wave2_19_source_aggregates_first_artifact_digests"]
    if len(artifact_rows)!=19 or [x["index"] for x in artifact_rows]!=list(range(19)):
        raise ValueError("Wave2 original first artifact SHA digests absent/out-of-order")
    files=sorted(source_directory.glob("wave2_source_*.json"))
    if len(files)!=19:
        raise ValueError("All nineteen original Wave2 first source-aggregate JSONs are required")
    w2_original={}
    for file in files:
        raw=json.loads(file.read_text(encoding="utf-8"))
        ix=raw.get("wave2_index")
        if not isinstance(ix,int) or ix not in range(19) or ix in w2_original:
            raise ValueError("Duplicated original Wave2 first source artifact")
        w2_original[ix]=raw
    if len(w2_original)!=19:
        raise ValueError("Incomplete Wave2 first source original source ID results")
    for i in range(19):
        source=w2["all_original_wave2_sources_in_rank_order"][i]
        rec=w2_original[i]
        doi=norm_doi(source["original_DOI_labels"][0])
        s=canonical[doi]
        prev=meta["all_51_source_study_rows"][i+16]
        if (rec.get("version")!="RSDS_WAVE2_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1"
            or rec.get("source_original_repository_study_id")!=s["study_id"]
            or norm_doi(rec.get("source_original_DOI_labels")[0])!=doi
            or rec.get("original_registered_native_taxa")!=s["taxa"]
            or rec.get("source_metadata_stratum")!=prev["source_stratum"]
            or rec.get("original_source_GPS_event_rows_opened")!=prev["original_source_event_row_count"]
            or rec.get("old_48_signed_solar_external_nonreplication_UNCHANGED")!=NEG
            or rec.get("old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED") is not True
            or rec.get("original_raw_GPS_dates_tag_IDs_coordinates_exported") is not False
            or rec.get("human_mobility_GPS_data_opened") is not False):
            raise ValueError("Wave2 SHA-first source artifact mismatched original 35-source meta")
        native=rec["native_group_results"]
        if len(native)!=s["n_native"] or [x["taxon"] for x in native]!=s["taxa"]:
            raise ValueError("Wave2 original registered native source roster changed")
        qvals=[q_effect(x["q"])["q"] for x in native]
        dayvals=[q_effect(x["daylength"])["q"] for x in native]
        validq=[v for v in qvals if v is not None]
        validday=[v for v in dayvals if v is not None]
        if (bool(prev["complete_source_native_q"])!=(len(validq)==s["n_native"])
            or bool(prev["complete_source_daylength"])!=(len(validday)==s["n_native"])
            or prev["native_q_supported"]!=len(validq)
            or prev["native_daylength_supported"]!=len(validday)
            or (prev["source_q_mean_equal_registered_native"] is not None and
                abs(mean(validq)-prev["source_q_mean_equal_registered_native"])>1e-12)):
            raise ValueError("Wave2 source mean q wrongly assigned to per-native or changed")
        for item in native:
            rows.append(make_row(s,item["taxon"],q_effect(item["q"]),
                                 q_effect(item["daylength"]),prev["complete_source_native_q"]))

    ci3_by={norm_doi(x["source_DOI_labels"][0]):x for x in
            ci3["new_16_source_original_native_ci_and_support_rows"]}
    for s in canonical.values():
        if s["wave"]!="wave3":
            continue
        original=ci3_by[s["doi"]]
        if original["original_registered_native_taxa"]!=s["taxa"]:
            raise ValueError("Wave3 original source-native CI label mismatch")
        mi=next((x for x in meta["all_51_source_study_rows"][35:]
                 if norm_doi(x["source_DOI_labels"][0])==s["doi"]),None)
        if mi is None:
            raise ValueError("Wave3 original source absent from frozen 51 study results")
        oq=from_wave3_ci(original["native_q"])
        od=from_wave3_ci(original["native_daylength"])
        if (bool(mi["complete_source_native_q"])!=oq["available"]
            or bool(mi["complete_source_daylength"])!=od["available"]
            or (oq["available"] and abs(oq["q"]-mi["source_q_mean_equal_registered_native"])>1e-12)
            or (od["available"] and abs(od["q"]-mi["source_daylength_mean_equal_registered_native"])>1e-12)):
            raise ValueError("Wave3 source-native original q/daylength changed")
        rows.append(make_row(s,s["taxa"][0],oq,od,mi["complete_source_native_q"]))

    if len(rows)!=88 or len({(x["original_DOI"],x["scientific_name_original_label"]) for x in rows})!=88:
        raise ValueError("Original 12(22 natives) + 51(66 natives) source-native frame changed")
    registered=defaultdict(list)
    for row in rows:
        row["original_catalog_title"]=canonical[row["original_DOI"]]["original_catalog_title"]
        registered[row["scientific_name_original_label"]].append(row)

    repeated=[]
    for taxon,native in sorted(registered.items()):
        if len(native)<2:
            continue
        have=[x for x in native if x["original_native_q"]["available"]]
        signs=[x["original_native_q"]["point_sign"] for x in have]
        bs=[x["original_native_q"]["nominal_source_CI_sign"] for x in have]
        distinct_doIs=len(set(x["original_DOI"] for x in native))
        same= len(have)>=2 and len(set(signs))==1 and signs[0] in ("POSITIVE","NEGATIVE")
        all_supported=same and all(z==signs[0] for z in bs)
        repeated.append({
            "scientific_name":taxon,
            "registered_distinct_source_DOI":distinct_doIs,
            "source_native_q_estimable_DOI":len(have),
            "source_native_q_nominal_source_CI_excludes_zero":sum(x in ("NEGATIVE","POSITIVE") for x in bs),
            "q_point_signs":dict(sorted(Counter(signs).items())),
            "q_original_source_CI_signs":dict(sorted(Counter(bs).items())),
            "same_q_point_sign_within_supported_archives":same if len(have)>=2 else None,
            "all_supported_source_CI_same_sign":all_supported if len(have)>=2 else None,
            "source_native_rows":native,
            "study_id_unique_but_original_tag_or_population_independence_verified":False,
        })
    repeated.sort(key=lambda x:(-x["source_native_q_estimable_DOI"],
                               -x["registered_distinct_source_DOI"],x["scientific_name"]))
    possible=[x for x in repeated if x["source_native_q_estimable_DOI"]>=2]
    count=Counter(x["source_wave"] for x in rows)
    if count!=Counter({"legacy12":22,"wave1":27,"wave2":23,"wave3":16}):
        raise ValueError("Source/native support and wave affiliation changed")
    return {
        "version":VERSION,
        "stage":"POST_OUTCOME_SOURCE_NATIVE_SAME_SPECIES_ORIGINAL_ARCHIVE_SYNTHESIS",
        "original_public_archive_2026_10_19_packages":437,
        "distinct_original_archived_source_DOIs_and_study_IDs":63,
        "source_studies_by_wave":{"legacy12":12,"wave1":16,"wave2":19,"wave3":16},
        "all_registered_original_source_native_rows":88,
        "all_original_source_native_q_supported_rows":sum(x["original_native_q"]["available"] for x in rows),
        "all_original_source_native_daylength_supported_rows":sum(x["original_native_daylength"]["available"] for x in rows),
        "repeated_species_labels_with_two_or_more_distinct_DOIs":len(repeated),
        "repeated_species_labels_with_two_or_more_estimable_source_native_q":len(possible),
        "same_species_source_study_q_comparison_candidates":possible,
        "all_repeated_source_native_species_including_original_HOLD":repeated,
        "source_native_first_wave2_original_artifact_digest_roster":artifact_rows,
        "source_native_original_51_meta_git_blob":c["original_first_51_study_meta_git_blob"],
        "previous_48_source_solar_irradiance_direction_replicated":False,
        "original_12_ERA5_raw_SHA_HOLD":True,
        "original_GPS_events_reopened":0,
        "individual_tag_ID_geographic_overlap_verified":False,
        "same_source_population_independence_verified":False,
        "source_pair_tag_reuse_possible":True,
        "original_source_nominal_CI_not_multistudy_inference":True,
        "ecological_claim_ceiling":c["population_independence_ceiling"],
        "study_population_sampling":"12 previously selected original source DOIs plus 51 metadata-priority Wave1–3 selected; outcome-adaptive quotas and nonrandom source completeness; do not count duplicate taxa within source as independent studies.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    out=build(a.directory)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_SAME_SPECIES_CROSS_ARCHIVE",json.dumps({
        "distinct_original_source_DOI":out["distinct_original_archived_source_DOIs_and_study_IDs"],
        "native_rows":out["all_registered_original_source_native_rows"],
        "repeated_name_labels":out["repeated_species_labels_with_two_or_more_distinct_DOIs"],
        "at_least_two_q_archives":out["repeated_species_labels_with_two_or_more_estimable_source_native_q"],
        "candidate_species":[{
            "taxon":x["scientific_name"],
            "sources":x["source_native_q_estimable_DOI"],
            "point_signs":x["q_point_signs"],
            "source_bootstrap_support":x["q_original_source_CI_signs"]}
            for x in out["same_species_source_study_q_comparison_candidates"]
        ],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
