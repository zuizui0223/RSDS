"""Source-only preflight of six frozen public Movebank GPS studies.

NO original animal time, location, speed, movement outcome or tag is read here.
No candidate is replaced based on provider file status or later outcome.
"""
from __future__ import annotations
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE, original_files

CONTRACT = HERE / "research/RSDS_NEW_SIX_ORIGINAL_GPS_V0_CONTRACT.json"
CATALOG = HERE / "research/seasonal_database_catalog/2026-10-19/movebank/catalog.json"
MATRIX = HERE / "research/GLOBAL_SOLAR_AZIMUTH_COMBINED_V2_ALL_ORIGINAL_SOURCE_TAXON_MATRIX_20261025.csv"
RADIATION_A = HERE / "research/SOLAR_AZIMUTH_24_SOURCE_MULTI_POPULATION_ERA5_V0_FIRST_ALL_RESULTS_20261025.json"
RADIATION_B = HERE / "research/SOLAR_RADIATION_EXTERNAL_24_SOURCE_FIRST_ALL_RESULTS_20261025.json"
PREVIOUS_24 = HERE / "research/GLOBAL_SOLAR_AZIMUTH_NEW24_PROVIDER_FROZEN_SOURCE_FILE_MANIFEST_20261025.json"
FORMATS = (".csv", ".zip", ".csv.gz")


def doi(s):
    return str(s).lower().removeprefix("doi:").removeprefix("https://doi.org/")


def validate_roster():
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    a = json.loads(CATALOG.read_text(encoding="utf-8"))
    if c["version"] != "RSDS_NEW_SIX_ORIGINAL_PUBLIC_GPS_JOINT_TEMPORAL_DIRECTION_V0":
        raise ValueError("Frozen pre-event scientific sampling route changed")
    if len(a["packages"]) != 437 or len(c["original_fixed_sources"]) != 6:
        raise ValueError("Frozen source study frame no longer complete")
    by_uuid = {x["package_uuid"]: x for x in a["packages"]}
    sources = c["original_fixed_sources"]
    if len({s["id"] for s in sources}) != 6 or len({s["study_id"] for s in sources}) != 6:
        raise ValueError("Six source study IDs are not unique")
    used_dois = set()
    used_ids = set()
    for x in json.loads(RADIATION_A.read_text(encoding="utf-8"))["source_groups"]:
        used_dois.update(map(doi, x["source_dois"]))
        used_ids.update(map(str, x["source_study_ids"]))
    for x in json.loads(RADIATION_B.read_text(encoding="utf-8"))["all_external_first_source_results"]:
        used_dois.update(map(doi, x["original_source_doIs"]))
        used_ids.update(map(str, x["source_study_ids"]))
    for x in json.loads(PREVIOUS_24.read_text(encoding="utf-8"))["source_file_entries"]:
        used_dois.update(map(doi, x["original_dois"]))
        used_ids.update(map(str, x["source_study_ids"]))
    with MATRIX.open(newline="", encoding="utf-8") as f:
        used_dois.update(doi(x["doi"]) for x in csv.DictReader(f))
    for s in sources:
        p = by_uuid.get(s["id"])
        if p is None:
            raise ValueError("Selected source not in pre-existing original 437 registry")
        if (p["scientific_names"] != [s["taxon"]]
                or p["repository_study_ids"] != [s["study_id"]]
                or s["doi"] not in set(map(doi, p["dois"]))
                or s["reported_animals"] != p["reported_animal_count"]
                or s["reported_fixes"] != p["reported_location_count"]
                or not any(x in ("GPS","GNSS") for x in p["sensor_names"])):
            raise ValueError("Original pre-event source-only species, DOI, study or sensor no longer matches")
        if s["doi"] in used_dois or s["study_id"] in used_ids:
            raise ValueError("Prospective six-study source accidentally reused a known GPS source study")
        if not p["licenses"]:
            raise ValueError("Original public source missing license metadata")
    return c, by_uuid


def one_source(s, provider):
    entry = {
        "package_uuid": s["id"], "original_dois": [s["doi"]],
        "source_study_ids": [s["study_id"]], "original_species_labels": [s["taxon"]],
        "original_licenses": provider["licenses"],
        "registry_reported_animals": s["reported_animals"],
        "registry_reported_locations": s["reported_fixes"],
        "source_recruitment_region_hint": s["region"],
        "locomotor_guild_hint": s["guild"],
        "animal_event_rows_opened_by_this_preflight": 0,
        "source_metadata_status": "SOURCE_METADATA_HOLD",
        "original_file": None,
    }
    try:
        native = original_files(s["id"])
        usable = [x for x in native
                  if str(x.get("name") or "").lower().endswith(FORMATS)
                  and "reference-data" not in x["name"].lower()
                  and "__macosx" not in x["name"].lower()
                  and not x["name"].rsplit("/",1)[-1].startswith("._")]
        if len(usable) != 1:
            entry["metadata_hold_reason"] = "NO_UNAMBIGUOUS_SINGLE_PRIMARY_ORIGINAL_FILE"
            entry["original_data_file_candidates"] = len(usable)
            return entry
        x = usable[0]
        cap = 200_000_000
        check = x.get("checksum") or {}
        if not isinstance(x.get("bytes"), int) or x["bytes"] <= 0 or x["bytes"] > cap:
            entry["metadata_hold_reason"] = "ORIGINAL_FILE_ABSENT_OR_OVER_FROZEN_200MB_BUDGET"
            return entry
        if (check.get("checkSumAlgorithm", "").upper() not in ("MD5","SHA256","SHA-256")
                or not check.get("value")):
            entry["metadata_hold_reason"] = "ORIGINAL_FILE_PROVIDER_CHECKSUM_UNSUPPORTED"
            return entry
        entry["original_file"] = {
            "bitstream_uuid": x["bitstream_uuid"],
            "filename": x["name"],
            "bytes": x["bytes"],
            "checksum": check,
        }
        entry["source_metadata_status"] = "ORIGINAL_FILE_METADATA_ADMISSIBLE"
    except Exception as ex:
        entry["metadata_hold_reason"] = type(ex).__name__
    return entry


def run():
    contract, registry = validate_roster()
    sources = contract["original_fixed_sources"]
    with ThreadPoolExecutor(max_workers=3) as pool:
        fs = {pool.submit(one_source, x, registry[x["id"]]): x["id"] for x in sources}
        by_id = {fs[f]: f.result() for f in as_completed(fs)}
    rows = [by_id[x["id"]] for x in sources]
    return {
        "version": "RSDS_NEW_SIX_ORIGINAL_GPS_FILE_METADATA_V0",
        "freeze_date": "2026-10-26",
        "fixed_source_studies": 6,
        "source_metadata_admissible": sum(x["source_metadata_status"] == "ORIGINAL_FILE_METADATA_ADMISSIBLE" for x in rows),
        "source_metadata_HOLD": sum(x["source_metadata_status"] != "ORIGINAL_FILE_METADATA_ADMISSIBLE" for x in rows),
        "original_animal_GPS_rows_read_by_preflight": 0,
        "source_entries": rows,
        "no_source_replacement_or_effect_selection": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("RSDS_NEW_SIX_METADATA", json.dumps({
        "studies": result["fixed_source_studies"],
        "usable": result["source_metadata_admissible"],
        "holds": result["source_metadata_HOLD"],
        "opened_gps": result["original_animal_GPS_rows_read_by_preflight"]
    }),flush=True)


if __name__ == "__main__":
    main()
