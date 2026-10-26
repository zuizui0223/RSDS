#!/usr/bin/env python3
"""Import *all relevant tracked* ODSP solar/sun/time/season ecology into RSDS.

Run inside a checked-out RSDS repository, with a SEPARATE checked-out ODSP tree.
Preserves raw bytes, source provenance and ALL positive, null, negative, HOLD
and unavailable results. Does NOT download or republish animal GPS records.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

SOURCE_REPO = "https://github.com/zuizui0223/odsp"
SOURCE_BRANCH = "ecology/solar-irradiance-direction-macro-gradient-v0"
SOURCE_COMMIT = "d58fb759416fc458f69ef61ea79b7b1d4164f2c5"
PATTERN = re.compile(
    r"solar|sun|season|photoperiod|circadian|era5|movebank|wyoming|utah|"
    r"serengeti|kruger|ireland|calendar|duck|equinoct|daylength|uljin|"
    r"irradiance|diel|astronom|temporal|etosha|masolele", re.I
)
EXTRA = {
    "docs/odsp_ungulates_two_independent_source_20261022.md",
    "docs/same_daylength_day_night_displacement_signatures_20261019.md",
    "scripts/collect_movebank_seasonal_catalog.py",
    "tests/test_movebank_seasonal_catalog.py",
}
ODSP_TEMPORAL_DEPENDENCIES = {
    "odsp/serengeti_terminal.py",
    "odsp/temporal_information.py",
    "odsp/temporal_partition.py",
    "odsp/temporal_crossfit.py",
    "odsp/crossfitted_transferability.py",
    "odsp/grouped_transferability.py",
    "odsp/niche_geometry.py",
    "odsp/transferability.py",
}
MANDATORY = {
    "manuscript/GLOBAL_SOLAR_AZIMUTH_ORIGINAL_GPS_MACROECOLOGY_v2_20261025.md",
    "manuscript/SOLAR_SEASONAL_HYSTERESIS_ECOLOGY_PAPER_v0.md",
    "research/SOLAR_RADIATION_TWO_24_SOURCE_NONREPLICATION_DIAGNOSTIC_FIRST_RESULT_20261025.json",
    "research/SOLAR_RADIATION_EXTERNAL_24_SOURCE_FIRST_ALL_RESULTS_20261025.json",
    "research/SOLAR_AZIMUTH_24_SOURCE_MULTI_POPULATION_ERA5_V0_FIRST_ALL_RESULTS_20261025.json",
    "research/SOLAR_IRRADIANCE_DIRECTION_MACRO_GRADIENT_V0_CONTRACT.json",
    "scripts/analyze_solar_irradiance_direction_macro_gradient_v0.py",
    "tests/test_solar_irradiance_direction_macro_gradient_v0.py",
    "research/seasonal_database_catalog/2026-10-19/movebank/catalog.json",
}
ROOT_RELATED = (
    "N2_SERENGETI_POST_TERMINAL_RUNNER_REPAIR.json",
    "N2_SERENGETI_RECOVERY_ATTEMPT2_FAILURE.json",
    "N2_SERENGETI_TECHNICAL_RECOVERY_CONTRACT.json",
    "N2_SERENGETI_TEMPORAL_TERMINAL_DECISION.json",
    "N2_SERENGETI_TEMPORAL_TERMINAL_RECEIPT.json",
    "N2_SERENGETI_TERMINAL_CLOSEOUT_CONTRACT.json",
)


def read_tracked_files(source: Path) -> list[str]:
    raw = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z"])
    return sorted(x.decode("utf-8") for x in raw.split(b"\x00") if x)


def source_is_eligible(path: str) -> bool:
    if path == "LICENSE":
        return True  # Preserve the originating code's MIT copyright notice.
    if path in ODSP_TEMPORAL_DEPENDENCIES or path in EXTRA or path in ROOT_RELATED:
        return True
    if path.startswith("research/seasonal_database_catalog/"):
        return True  # Includes the source registry's original metadata-page tar.gz.
    if path.startswith(".github/workflows/"):
        return PATTERN.search(path) is not None
    prefix = path.partition("/")[0]
    if prefix not in {"docs", "manuscript", "research", "scripts", "tests"}:
        return False
    return PATTERN.search(path) is not None


def target_path(source_path: str) -> str:
    if source_path.startswith(".github/workflows/"):
        return "archive/odsp_workflows/" + source_path.rsplit("/", 1)[-1]
    if source_path == "LICENSE":
        return "LICENSE"
    return source_path


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def migrate(source: Path, dest: Path) -> dict:
    source = source.resolve()
    dest = dest.resolve()
    if source == dest or source in dest.parents or dest in source.parents:
        raise ValueError("ODSP and RSDS must be separate working trees")
    actual = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if actual != SOURCE_COMMIT:
        raise ValueError(f"Refuse mutable source: expected {SOURCE_COMMIT}, got {actual}")
    names = read_tracked_files(source)
    selected = [name for name in names if source_is_eligible(name)]
    if not MANDATORY.issubset(set(selected)):
        raise ValueError("Required original ecological evidence/code/registry is absent")
    if len(selected) < 370:
        raise ValueError(f"Unexpectedly incomplete scientific transfer: {len(selected)}")
    seen = set()
    entries = []
    for name in selected:
        target = target_path(name)
        if target in seen:
            raise ValueError("Two original files mapped to one RSDS path")
        seen.add(target)
        old = source / name
        new = dest / target
        if not old.is_file() or old.is_symlink():
            raise ValueError(f"Unrecognized original file type: {name}")
        payload = old.read_bytes()
        new.parent.mkdir(parents=True, exist_ok=True)
        new.write_bytes(payload)
        if new.read_bytes() != payload:
            raise ValueError(f"Byte-for-byte transfer failed: {name}")
        entries.append({
            "odsp_path": name,
            "rsds_path": target,
            "bytes": len(payload),
            "sha256": sha256(payload),
        })
    manifest = {
        "schema": "rsds-odsp-solar-season-research-transfer-v1",
        "source_repository": SOURCE_REPO,
        "source_branch": SOURCE_BRANCH,
        "source_commit": SOURCE_COMMIT,
        "transfer_date_utc": "2026-10-26",
        "source_selection": "All tracked solar, sun, time, season, radiation, ERA5, "
            "Movebank, calendar, photoperiod and Serengeti research files; "
            "ODSP temporal algorithm dependencies; full published discovery catalogs.",
        "original_animal_telemetry_redistributed": False,
        "archived_workflows_location": "archive/odsp_workflows/",
        "files": entries,
        "file_count": len(entries),
        "total_bytes": sum(x["bytes"] for x in entries),
        "note": "Source snapshots include positive, negative, unavailable, source-HOLD "
            "and exploratory findings; archived ODSP workflows are NOT RSDS active CI.",
    }
    out = dest / "migration/ODSP_TRANSFER_MANIFEST.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    # Only a scoped compatibility namespace. We do not import full unrelated ODSP.
    pkg = dest / "odsp/__init__.py"
    pkg.write_text(
        '"""RSDS-scoped legacy ODSP temporal-ecology dependencies.\n'
        'These modules preserve previously frozen original temporal estimands.\n'
        'General ODSP prediction/forecast modules have not been migrated.\n'
        '"""\n', encoding="utf-8"
    )
    print("RSDS_TRANSFER", json.dumps({
        "files": len(entries), "bytes": manifest["total_bytes"],
        "source_commit": SOURCE_COMMIT,
        "workflow_archive_count": sum(
            x["rsds_path"].startswith("archive/odsp_workflows/") for x in entries),
        "negative_and_hold_originals_retained": True,
        "original_gps_events_copied": False,
    }, sort_keys=True))
    return manifest


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--dest", default=".", type=Path)
    args = p.parse_args()
    migrate(args.source, args.dest)
