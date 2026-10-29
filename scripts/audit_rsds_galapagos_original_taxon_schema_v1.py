"""Galapagos original Movebank taxon/schema label diagnostic only.

IMPORTANT: The prior source's exact registered species parser matched ZERO native
taxon rows, but the original publisher 137MB source MD5 was verified. This
reads only the ORIGINAL HEADER and candidate taxonomy COLUMN from the already
provider-checksummed original CSV. The row parser does not inspect or retain
the other values (individual tags, GPS position, event dates, q, orientation).
No source substitution, synonym rescue or biological q reclassification.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)

CONTRACT=HERE/"research/RSDS_GALAPAGOS_SOURCE_ONLY_TAXON_SCHEMA_V1_CONTRACT.json"
FIRST=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_GPS_SHARED_CALENDAR_V1_FIRST_RESULT_20261029.json"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
EXPECTED_TAXON=(
    "Chelonoidis donfaustoi","Chelonoidis hoodensis","Chelonoidis porteri"
)


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    b=FIRST.read_bytes()
    sha=hashlib.sha1(("blob "+str(len(b))+"\0").encode("ascii")+b).hexdigest()
    previous=json.loads(b.decode("utf-8"))
    if (c["version"]!="RSDS_GALAPAGOS_PUBLISHER_ORIGINAL_TAXON_LABEL_AUDIT_V1"
            or sha!=c["frozen_original_source_first_result_git_blob_SHA"]
            or previous["fixed_original_published_studies"]!=4
            or previous["registered_original_study_native_taxon_rows"]!=10
            or not previous["prior_12_ERA5_version_source_HOLDS_UNCHANGED"]
            or previous["prior_48_source_shortwave_direction_NEGATIVE_EXTERNAL_REPLICATION_UNCHANGED"]!=NEG):
        raise ValueError("Previously fixed first source studies, taxon counts or external solar status changed")
    source=previous["all_frozen_source_study_original_results"][0]
    if (source["published_DOI"]!=c["source_study"]["doi"]
            or source["source_status"]!="SOURCE_NO_VALID_GPS_FIXES"
            or {z["taxon"] for z in source["original_native_source_q_and_daylength"]}!=set(EXPECTED_TAXON)
            or any(z["original_source_metadata_taxonomy_coverage"]["original_taxon_rows"]!=0
                   for z in source["original_native_source_q_and_daylength"])):
        raise ValueError("Galapagos previous exact-taxonomy source HOLD reclassified")
    if tuple(c["source_study"]["registered_taxa"])!=EXPECTED_TAXON:
        raise ValueError("Original publisher metadata native taxon frame was replaced")
    return c


def normal(s):
    return "-".join(s.lower().replace("_","-").replace(" ","-").split("-"))


def taxonomy_only_file(path,contract):
    """Only the selected taxon label field is ever accessed from raw CSV row."""
    csv.field_size_limit(1_000_000)
    with path.open("r",newline="",encoding="utf-8-sig") as stream:
        parser=csv.reader(stream)
        header=next(parser,None)
        if not header:
            return {"taxonomy_state":"ORIGINAL_CSV_HEADER_MISSING","rows":0}
        names=[normal(x.strip()) for x in header]
        permitted=[normal(x) for x in contract["allowed_taxon_column_names"]]
        labels=[x for x in permitted if x in names]
        require=["timestamp","individual-local-identifier","location-lat","location-long"]
        gp_schema={k:k in names for k in require}
        if len(labels)!=1:
            return {
                "taxonomy_state":"PUBLISHER_TAXON_COLUMN_MISSING_OR_AMBIGUOUS",
                "source_original_total_columns":len(header),
                "requested_taxonomy_field_candidates_present":labels,
                "original_tracking_GPS_expected_schema_column_present":gp_schema,
                "taxonomy_rows_scanned":0,"no_original_animal_positions_read":True
            }
        taxonfield=labels[0]
        idx=names.index(taxonfield)
        counts=Counter()
        missing=0
        n=0
        maxrows=contract["limits"]["max_rows"]
        maxlabels=contract["limits"]["max_taxon_labels"]
        for row in parser:
            n+=1
            if n>maxrows:
                raise ValueError("Original publisher source >2m rows; no partial taxon diagnostic")
            label=row[idx].strip() if idx<len(row) else ""
            if not label:
                missing+=1
                continue
            # Only plausible public TAXON labels. Do not accidentally emit
            # source animal IDs stored under an incorrectly named field.
            if len(label)>100 or not all(ch.isalpha() or ch in " ._-()/0123456789" for ch in label):
                raise ValueError("Unexpected non-taxonomic source field; keep taxonomy label HOLD")
            counts[label]+=1
            if len(counts)>maxlabels:
                return {
                    "taxonomy_state":"HIGH_CARDINALITY_SOURCE_TAXON_LABEL_HOLD",
                    "source_original_taxon_column":taxonfield,
                    "observed_distinct_labels_lower_bound":len(counts),
                    "source_native_labels_exported":False,
                    "no_original_animal_positions_read":True
                }
    declared=set(contract["source_study"]["registered_taxa"])
    return {
        "taxonomy_state":"PUBLIC_SOURCE_TAXON_LABEL_COUNTS_AUDITED",
        "source_original_total_columns":len(header),
        "source_original_taxon_column":taxonfield,
        "original_tracking_GPS_expected_schema_column_present":gp_schema,
        "taxonomy_rows_scanned":n,
        "rows_with_empty_source_taxon_label":missing,
        "distinct_original_taxon_labels":len(counts),
        "original_taxon_row_counts":[
            {"raw_source_taxon_label":x,
             "source_original_rows":counts[x],
             "exact_registered_taxon_name_match":x in declared}
            for x in sorted(counts)
        ],
        "total_rows_exact_registered_taxon_name_match":
            sum(counts[x] for x in counts if x in declared),
        "original_taxon_labels_not_automatically_reinterpreted":True,
        "no_original_animal_positions_read":True
    }


def run(out):
    c=frozen()
    src=c["source_study"]
    result={
        "version":"RSDS_GALAPAGOS_PUBLISHER_ORIGINAL_SOURCE_TAXON_SCHEMA_V1_FIRST_RESULT",
        "date":"2026-10-29",
        "original_study_DOI":src["doi"],
        "original_source_UUID":src["package_uuid"],
        "original_MD5_verified":False,
        "source_file_choice_fixed_before_original_q":True,
        "original_q_and_photoperiod_first_result_UNCHANGED":True,
        "taxonomy_state":"SOURCE_PROVIDER_HOLD",
        "taxonomy_diagnosis":None,
        "old_12_ERA5_raw_JSON_SHA_HOLDS_UNCHANGED":True,
        "old_48_signed_radiation_direction_external_nonreplication_UNCHANGED":NEG,
        "no_individual_GPS_tag_time_coordinate_or_effect_exported":True
    }
    out.mkdir(parents=True,exist_ok=True)
    file=out/"private_galapagos_original_taxon_only.csv"
    try:
        matches=[x for x in original_files(src["package_uuid"])
                 if x["bitstream_uuid"]==src["bitstream_uuid"]]
        if (len(matches)!=1 or matches[0]["name"]!=src["original_file_name"]
                or matches[0]["bytes"]!=src["file_bytes"]
                or matches[0]["checksum"]!=src["checksum"]):
            raise ValueError("Pre-event original tortoise file provider MD5/size/UUID no longer matches")
        status=checksum_and_download(matches[0],file)
        if not status["provider_checksum_verified"]:
            raise ValueError("Original publisher file full MD5 not verified")
        result["original_MD5_verified"]=True
        output=taxonomy_only_file(file,c)
        result["taxonomy_diagnosis"]=output
        result["taxonomy_state"]=output["taxonomy_state"]
    except Exception as err:
        result["taxonomy_state"]="TAXON_SCHEMA_PROVIDER_OR_SOURCE_HOLD"
        result["taxonomy_diagnostic_error_type"]=type(err).__name__
    finally:
        file.unlink(missing_ok=True)
    target=out/"galapagos_publisher_taxon_only_schema.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",
                      encoding="utf-8")
    print("RSDS_GALAPAGOS_ORIGINAL_TAXON_SCHEMA",json.dumps({
        "doi":src["doi"],
        "MD5":result["original_MD5_verified"],
        "state":result["taxonomy_state"],
        "taxon_labels":(result["taxonomy_diagnosis"] or {}).get("distinct_original_taxon_labels"),
        "original_rows":(result["taxonomy_diagnosis"] or {}).get("taxonomy_rows_scanned"),
        "previous_q_first_result":"UNCHANGED"
    },sort_keys=True),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    run(args.out)


if __name__=="__main__":
    main()
