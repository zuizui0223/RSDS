"""Aggregate two fixed original tagged-animal context sensitivity runs.

This preserves every original study×taxon and unavailable result. Contexts
are ORIGINAL GPS-coverage ranked, not environmental high/low.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_CONTRACT.json"
BASELINE = ROOT / "research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def aggregate(utah, wyoming):
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    if (contract["version"] != "RSDS_SAME_TAG_CONTEXT_PAIR_V1"
            or baseline["frozen_external_radiation_nonreplication_UNCHANGED"] != NEGATIVE):
        raise ValueError("Original scientific source and historical validation not fixed")
    sample = {"Utah": utah, "Wyoming": wyoming}
    if any(x.get("region") != region or x.get("version") != contract["version"]
           for region, x in sample.items()):
        raise ValueError("Wrong source tagged-animal identity")
    expected = {s["study"]: s for s in contract["source_frame"]}
    if set(expected) != set(sample):
        raise ValueError("Missing original published study")
    rows = []
    states = {}
    for region in ("Utah", "Wyoming"):
        source = sample[region]
        if source["original_source_doi"] != expected[region]["doi"]:
            raise ValueError("Wrong original Movebank DOI for this region")
        taxa = source["source_taxa"]
        if {x["taxon"] for x in taxa} != set(expected[region]["taxa"]):
            raise ValueError("Incomplete 6 source×taxon frame")
        states[region] = source["state"]
        for item in taxa:
            effect = item.get("paired_tag_effect")
            if effect is not None:
                for axis, x in effect["per_axis"].items():
                    if x["status"] not in ("PAIRED_ESTIMABLE", "INSUFFICIENT_PAIRED_ANIMALS"):
                        raise ValueError("Undefined paired original source status")
                    if x["status"] == "PAIRED_ESTIMABLE" and x["paired_source_animals"] < 5:
                        raise ValueError("Effect reported without 5 original distinct paired tags")
            rows.append({
                "original_region": region,
                "source_doi": source["original_source_doi"],
                "taxon": item["taxon"],
                "original_native_q_replay_before_ERA5":
                    item.get("original_source_replay_success", False),
                "status": item.get("new_paired_sensitivity_status", "SOURCE_HOLD"),
                "selected_contexts_exact_original_replay":
                    item.get("frozen_first_contexts_verified", False),
                "original_contexts": item.get("context_coverage", []),
                "paired_same_tag_effect": effect,
            })
    if len(rows) != 6:
        raise ValueError("Not the two original studies / six original native source taxa")
    paired_n = sum(row["paired_same_tag_effect"] is not None
                   for row in rows)
    available_by_axis = {
        axis: sum(
            row["paired_same_tag_effect"] is not None and
            row["paired_same_tag_effect"]["per_axis"][axis]["status"] == "PAIRED_ESTIMABLE"
            for row in rows
        )
        for axis in ("sun_cos", "sun_sin", "q")
    }
    return {
        "version": "RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT",
        "date": "2026-10-26",
        "source_first_result": "RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026",
        "analysis_stage": "POST_SOURCE_RESULT_PAIRED_TAG_CONTEXT_SENSITIVITY",
        "source_original_studies": 2,
        "native_study_taxon_rows": len(rows),
        "source_original_states": states,
        "source_taxa_with_two_comparable_original_contexts": paired_n,
        "number_of_source_taxa_with_estimate_by_axis": available_by_axis,
        "all_native_source_taxon_results": rows,
        "heldout_negative_radiation_nonreplication_UNCHANGED": NEGATIVE,
        "claim_ceiling": [
            "Difference is rank2-minus-rank1 ORIGINAL q-support-ranked gridyear context, not brighter-versus-darker, early-versus-late, or randomly assigned environment.",
            "Same tagged animal in two contexts is not a controlled field experiment and cannot establish sunlight, wind or daylength causality.",
            "Across-source taxon comparisons may share landscape, monitoring schemes and calibration.",
            "No per-individual original identities, exact source-grid coordinate, original GPS fix or event timestamps are distributed.",
            "Source weather SHA holds are not ecological zeros. Any archived first v0 source weather change is retained and not silently substituted.",
            "Original v0 original results and previously preregistered external negative radiation nonreplication are unchanged."
        ]
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input", type=Path, required=True)
    cli.add_argument("--out", type=Path, default=ROOT / "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT_20261026.json")
    args = cli.parse_args()
    files = {}
    for name in ("utah", "wyoming"):
        path = args.input / ("source_" + name + "_paired_context.json")
        if not path.is_file():
            raise ValueError("Original fixed source result absent: " + name)
        files[name] = json.loads(path.read_text(encoding="utf-8"))
    result = aggregate(files["utah"], files["wyoming"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
                        encoding="utf-8")
    print("RSDS_SAME_TAG_PAIR_ATLAS",
          json.dumps({
              "studies": result["source_original_studies"],
              "study_taxa": result["native_study_taxon_rows"],
              "qualified_axis_counts": result["number_of_source_taxa_with_estimate_by_axis"],
              "source_states": result["source_original_states"],
              "heldout_verdict": result["heldout_negative_radiation_nonreplication_UNCHANGED"]
          }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
