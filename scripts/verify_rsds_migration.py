"""Verify migrated ODSP solar/season ecology assets without any network or telemetry.

This scientific integrity guard checks every original Git blob (SHA-1), source
manifest record, per-file receipt and explicitly frozen external negative-
direction nonreplication. RSDS science remains distinct from ODSP prediction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance/ODSP_SOLAR_SEASON_MIGRATION_MANIFEST_20261026.json"
RECEIPT = ROOT / "provenance/RSDS_TRANSFER_VERIFICATION_20261026.json"

EXPECTED_SOURCE_COMMIT = "d58fb759416fc458f69ef61ea79b7b1d4164f2c5"
EXPECTED_NEGATIVE_EXTERNAL = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def blob_sha_bytes(data: bytes) -> str:
    """Original Git SHA-1 object hash of the raw unmodified file bytes."""
    return hashlib.sha1(("blob " + str(len(data)) + "\0").encode("ascii") + data).hexdigest()


def validate(root: Path = ROOT) -> dict:
    root = root.resolve()
    manifest = json.loads((root / MANIFEST.relative_to(ROOT)).read_text(encoding="utf-8"))
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_text(encoding="utf-8"))
    if manifest["source_commit"] != EXPECTED_SOURCE_COMMIT:
        raise ValueError("The original source commit identity was changed")
    if manifest["source_repository"] != "zuizui0223/odsp" or manifest["target_repository"] != "zuizui0223/RSDS":
        raise ValueError("Source or target repository identity changed")
    expected = manifest["files"]
    if len(expected) != 390 or manifest["text_items"] != 390:
        raise ValueError("Original frozen 390-item transfer frame differs")
    if receipt["origin_commit"] != EXPECTED_SOURCE_COMMIT or not receipt["all_manifest_items_verified"]:
        raise ValueError("Migration receipt is not a passed frozen-source receipt")
    if receipt["exact_original_blob_sha_matches"] != len(expected):
        raise ValueError("Receipt no longer covers all source files")
    if len({p["target"] for p in expected}) != len(expected):
        raise ValueError("Duplicate original destination path")
    receipts = {r["path"]: r for r in receipt["files"]}
    if len(receipts) != len(expected):
        raise ValueError("Receipt has an incomplete, duplicate or extra file path")
    verified_bytes = 0
    for record in expected:
        relative = Path(record["target"])
        full_path = (root / relative).resolve()
        if not full_path.is_relative_to(root) or not full_path.is_file():
            raise ValueError("Original source file missing from RSDS or unsafe path: " + str(relative))
        raw = full_path.read_bytes()
        if len(raw) != record["bytes"] or blob_sha_bytes(raw) != record["blob_sha"]:
            raise ValueError("ODSP original blob SHA mismatch: " + str(relative))
        rr = receipts.get(record["target"])
        if rr != {"path": record["target"], "git_blob_sha": record["blob_sha"],
                  "verified_bytes": record["bytes"], "status": "exact_match"}:
            raise ValueError("Source-identity provenance receipt mismatch: " + str(relative))
        verified_bytes += len(raw)
    if verified_bytes != manifest["source_text_bytes"] or verified_bytes != receipt["exact_original_source_bytes"]:
        raise ValueError("Migrated source-byte totals have changed")
    if not receipt["historic_ODSP_workflows_inert_under_archive"] or receipt["original_animal_GPS_fixes_copied"]:
        raise ValueError("Unsafe original workflow/raw GPS migration mode")

    original_nonreplication = json.loads(
        (root / "research/SOLAR_RADIATION_TWO_24_SOURCE_NONREPLICATION_DIAGNOSTIC_FIRST_RESULT_20261025.json")
        .read_text(encoding="utf-8")
    )
    if original_nonreplication["prior_frozen_external_replication_decision_UNCHANGED"] != EXPECTED_NEGATIVE_EXTERNAL:
        raise ValueError("Original failed heldout negative-direction replication was reclassified")
    if original_nonreplication["all_48_previously_frozen_source_studies_included"] != 48:
        raise ValueError("Original external study frame changed")

    for file in ("README.md", "docs/RSDS_ECOLOGICAL_SYNTHESIS_20261026.md",
                 "LICENSE", "odsp/__init__.py"):
        if not (root / file).is_file():
            raise ValueError("RSDS handoff/compatibility component missing: " + file)
    return {
        "frozen_ODSP_commit": EXPECTED_SOURCE_COMMIT,
        "verified_files": len(expected),
        "verified_bytes": verified_bytes,
        "external_solar_radiation_directional_verdict": EXPECTED_NEGATIVE_EXTERNAL,
        "rsds_is_scientific_primary_and_old_results_unchanged": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    print(json.dumps(validate(args.root), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
