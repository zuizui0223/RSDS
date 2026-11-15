"""Source-pinned original GPS runner for four existing co-calendar field studies.

Original GPS accessed only in ephemeral Action tempfile. Publisher UUID,
byte count and MD5 are checked BEFORE opening tracking responses. Previously
archived native effects and shared calendar support must replay. An error
returns explicit HOLD and never replaces the source or native species.
No tracking events/individual IDs/dates/coordinates are exported.
"""
from __future__ import annotations

import argparse
import json
import math
import tempfile
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import (
    original_files,checksum_and_download
)
from scripts import run_rsds_illinois_three_native_shared_calendar_v1 as ill
from scripts import run_rsds_utah_coobserved_photoperiod_v1 as ut
from scripts import run_rsds_wyoming_coobserved_photoperiod_v1 as wy
from scripts import run_rsds_four_new_multinative_original_gps_temporal_v1 as multi
from scripts.rsds_original_joint_calendar_covariance_v1 import (
    B,SEED,NEGATIVE,frozen_parent,summarize_joint
)

EXPECTED = {
    "Illinois":{"doi":"10.5441/001/1.649","days":169,"months":9,"taxa":3},
    "Wyoming":{"doi":"10.5441/001/1.703","days":193,"months":8,"taxa":2},
    "Utah":{"doi":"10.5441/001/1.711","days":199,"months":10,"taxa":2},
    "Geese":{"doi":"10.5441/001/1.fk899541","days":120,"months":5,"taxa":3},
}

def source_file(package_uuid,bitstream_uuid,expected_bytes,md5,path):
    sources=[x for x in original_files(package_uuid)
             if x["bitstream_uuid"]==bitstream_uuid]
    if len(sources)!=1 or sources[0]["bytes"]!=expected_bytes:
        raise ValueError("Original publisher UUID/bytes mismatch")
    original=sources[0]
    expected=md5.lower()
    if original["checksum"]["value"].lower()!=expected:
        raise ValueError("Original provider MD5 metadata mismatch")
    receipt=checksum_and_download(original,path)
    if receipt.get("provider_checksum_verified") is not True:
        raise ValueError("Original source bytes did not match MD5")
    return original

def source_illinois(path,prior):
    spec,previous,old_solar=ill.freeze()
    source=spec["source_original_study"]
    source_file(source["package_uuid"],source["bitstream_uuid"],
                source["bytes"],source["md5"],path)
    tracks,coverage,total=ill.load_source_taxa(path,source["file_name"],list(ill.TAXA))
    ill.check_original_source_coverage(tracks,source,old_solar,total)
    for sp in ill.TAXA:
        result,_=ill.one_taxon(sp,tracks[sp])
        ref=previous["source_native_taxa"][sp]
        for label,key in (("day_night_speed_contrast","q"),
                          ("daylength_association","q_daylength")):
            v=result[label]
            old=ref[key]
            if (v is None or old["mean"] is None or
                v["n_individuals"]!=old["source_tags"] or
                not math.isclose(v["mean"],old["mean"],abs_tol=1e-12,rel_tol=0)):
                raise ValueError("Historical Illinois original native q replay changed")
    if prior["status"]!="THREE_NATIVE_SHARED_SOURCE_DATE_EVALUATED":
        raise ValueError("Previously locked three native Illinois source not evaluable")
    pool={sp:ill.private_tag_calendar(tracks[sp],"original") for sp in ill.TAXA}
    co,dates,_=ill.triple_common_calendar(pool)
    if {s:len(co[s]) for s in ill.TAXA}!=prior["source_all_three_native_support"][
        "native_original_tag_count_on_all_three_shared_source_dates"]:
        raise ValueError("Illinois source native tag support changed")
    return source["doi"],co,dates,list(ill.TAXA),ill.REFERENCE

def source_ungulate(region,path,prior):
    library = ut if region=="Utah" else wy
    spec,prev,prev_delta=library.frozen()
    old,new,_,_=library.verify_frame()
    source,taxa,reference,_=library.source_reference(region,old,new)
    if source["doi"]!=EXPECTED[region]["doi"]:
        raise ValueError("Original two-source DOI changed")
    source_file(source["package_uuid"],source["original_bitstream_uuid"],
                source["original_bytes"],source["original_md5"],path)
    if region=="Utah":
        tracks,coverage,total=library.extract_time_only_original_source(path,taxa)
    else:
        tracks,coverage,total=library.streaming_time_selected_csv(path,taxa)
    library.source_baseline_check(tracks,reference)
    library.replay_original_triple(tracks,prev,prev_delta)
    if (prior["status"]!="STUDY_COMMON_DATE_EVALUATED"
        or prior["both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed"] is not True):
        raise ValueError("Archived original ungulate study results not evaluated")
    selected={sp:library.as_three_original_calendar(tracks[sp]) for sp in library.TAXA}
    common,dates,_=library.matched_two_taxon_calendar(selected)
    if {sp:len(common[sp]) for sp in library.TAXA}!=prior["shared_support"][
        "original_tags_with_10_days_on_both_taxa_shared_calendar"]:
        raise ValueError("Original shared-date ungulate native tag counts changed")
    return source["doi"],common,dates,list(library.TAXA),library.REFERENCE_LATITUDE

def source_geese(path,prior):
    spec,manifest=multi.fixed_sources()
    study=spec["sources"][3]
    file=manifest[3]
    receipt=file["source_original_file"]
    if (file["source_status"]!="SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE"
        or study["doi"]!=EXPECTED["Geese"]["doi"]
        or sorted(study["registered_taxa"]) != sorted(prior_row["taxon"] for
              x in prior["all_frozen_source_study_original_results"]
              if x["published_DOI"]==study["doi"]
              for prior_row in x["original_native_source_q_and_daylength"])):
        raise ValueError("Original source-only registered geese taxonomy changed")
    source_file(study["package_uuid"],receipt["bitstream_uuid"],
                receipt["bytes"],receipt["checksum"]["value"],path)
    taxa=study["registered_taxa"]
    tracks,coverage,total=multi.load_source_taxa(path,receipt["filename"],taxa)
    parent=next(x for x in prior["all_frozen_source_study_original_results"]
                if x["published_DOI"]==study["doi"])
    if (parent["source_status"]!="ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED"
        or not parent["source_provider_original_MD5_replayed"]):
        raise ValueError("Historical geese original GPS was not evaluated")
    old_native={x["taxon"]:x for x in parent["original_native_source_q_and_daylength"]}
    original={}
    for sp in taxa:
        current,_=multi.one_taxon(sp,tracks[sp])
        prior_native=old_native[sp]
        for k,r in (("day_night_speed_contrast","original_q_full_source_tag_effect"),
                    ("daylength_association","original_q_daylength_full_source_tag_effect")):
            one=current[k]
            frozen=prior_native[r]
            if (one is None or frozen is None or
                one["n_individuals"]!=frozen["n_individuals"] or
                not math.isclose(one["mean"],frozen["mean"],abs_tol=1e-12,rel_tol=0)):
                raise ValueError("Original 3-geese native q/daylength first effect changed")
        original[sp],_=multi.per_tag_year_days(tracks[sp],"ORIGINAL_1H")
    latitude=multi.source_reference_latitude(tracks)
    if latitude is None: raise ValueError("Original geese study missing astronomical latitude")
    common,dates,_=multi.common_original_calendar(original)
    support=parent["original_all_native_common_calendar_support"]
    if (len(dates)!=support["source_calendar_dates_shared_by_all_registered_taxa"]
        or {sp:len(common[sp]) for sp in taxa}!=support["source_registered_native_q_tag_counts"]):
        raise ValueError("Historical geese date/native tag support changed")
    converted={
        sp:{tag:{yr:{d:{"original":float(q[0])}
                    for d,q in dates.items()} for yr,dates in years.items()}
            for tag,years in common[sp].items()}
        for sp in taxa
    }
    return study["doi"],converted,dates,taxa,latitude

def compute(study,nboot=B,seed=SEED):
    info=EXPECTED.get(study)
    if info is None: raise ValueError("Only four prior co-calendar studies")
    meta={
        "version":"RSDS_ORIGINAL_GPS_JOINT_MONTH_ANIMAL_COVARIANCE_V1",
        "study":study,"original_DOI":info["doi"],
        "status":"SOURCE_COVARIANCE_HOLD",
        "original_source_GPS_MD5_verified":False,
        "original_native_q_daylength_effect_replayed":False,
        "historic_shortwave_negative_external_nonreplication_UNCHANGED":NEGATIVE,
        "historic_ERA5_raw_JSON_12_HOLD_UNCHANGED":True,
        "private_events_exported":False,
    }
    stage="FROZEN_PREVIOUS_BLOBS"
    try:
        previous,expected=frozen_parent(study)
        with tempfile.TemporaryDirectory(prefix="rsds_private_gps_") as directory:
            local=Path(directory)/("original.zip" if study=="Utah" else
                                   "original.csv.gz" if study=="Geese" else "original.csv")
            stage="ORIGINAL_PROVIDER_UUID_BYTES_MD5_AND_REPLAY"
            if study=="Illinois":
                doi,pool,dates,taxa,lat=source_illinois(local,previous)
            elif study in ("Utah","Wyoming"):
                doi,pool,dates,taxa,lat=source_ungulate(study,local,previous)
            else:
                # Preserve original declared ZIP/CSV format; provider's source
                # filename determines the actual original container suffix.
                _,entries=multi.fixed_sources()
                filename=entries[3]["source_original_file"]["filename"]
                local=Path(directory)/("original"+(
                    ".zip" if filename.lower().endswith(".zip") else
                    ".csv.gz" if filename.lower().endswith(".csv.gz") else ".csv"))
                doi,pool,dates,taxa,lat=source_geese(local,previous)
            meta["original_source_GPS_MD5_verified"]=True
            meta["original_native_q_daylength_effect_replayed"]=True
            stage="FIXED_SOURCE_DATE_MONTH_AND_COVARIANCE"
            if (doi!=info["doi"] or len(taxa)!=info["taxa"] or
                len(dates)!=info["days"] or
                len({(yr,d.month) for yr,d in dates})!=info["months"]):
                raise ValueError("Fixed source native/calendar support diverged")
            result=summarize_joint(study,doi,pool,dates,taxa,lat,
                                   previous,expected,nboot=nboot,seed=seed)
        result.update({
            "original_source_GPS_MD5_verified":True,
            "original_native_q_daylength_effect_replayed":True,
            "private_events_exported":False,
        })
        return result
    except Exception as exc:
        meta["hold_stage"]=stage
        meta["hold_error_type"]=type(exc).__name__
        # Never publish the exception string: provider messages may contain
        # private paths, original locations, native animal IDs or timestamps.
        return meta

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",choices=tuple(EXPECTED),required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=compute(a.study)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_JOINT_SOURCE_COVARIANCE",json.dumps({
        "study":a.study,"status":result["status"],
        "MD5":result["original_source_GPS_MD5_verified"],
        "replayed":result["original_native_q_daylength_effect_replayed"],
        "month_blocks":result.get("eligible_common_year_month_blocks"),
        "hold_stage":result.get("hold_stage"),
        "hold_type":result.get("hold_error_type"),
    },sort_keys=True),flush=True)
    # These four studies were historically feasible; a new HOLD must not
    # silently pass source CI or be substituted with a convenient source.
    if result["status"]!="ORIGINAL_SOURCE_JOINT_COVARIANCE_ESTIMATED":
        raise SystemExit(2)

if __name__ == "__main__":
    main()
