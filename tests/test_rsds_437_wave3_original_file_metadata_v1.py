"""Wave3 original-provider metadata-only screen and privacy/hold unit tests."""
from __future__ import annotations

import unittest
from unittest.mock import patch

from scripts import rsds_437_wave3_original_file_metadata_v1 as m


def original_file(name="original_gps.csv", size=12_000_000, md5="a" * 32):
    return {
        "name": name,
        "bytes": size,
        "bitstream_uuid": "known-uuid",
        "checksum": {"checkSumAlgorithm": "MD5", "value": md5},
    }


class Wave3SourceProviderMetadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.selected = m.freeze()

    def test_source_rank_and_frozen_file_preflight(self):
        self.assertEqual(len(self.selected), 16)
        self.assertEqual(
            [r["source_rank_1_based"] for r in self.selected], list(range(25, 41))
        )
        self.assertEqual(
            len({r["original_repository_study_id"] for r in self.selected}), 16
        )

    def test_admissible_MD5_file_only(self):
        s = self.selected[0]
        with patch.object(m, "original_files", return_value=[original_file()]) as fetch:
            report = m.inspect(s)
            fetch.assert_called_once_with(s["package_uuid"])
        self.assertEqual(report["source_status"], "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        self.assertEqual(report["original_file"]["checksum"]["value"], "a" * 32)
        self.assertEqual(report["original_file"]["bytes"], 12_000_000)
        self.assertTrue(report["provider_original_file_list_obtained"])
        self.assertEqual(report["animal_GPS_events_opened"], 0)

    def test_human_location_never_queries_file_listing(self):
        human = dict(self.selected[0], scientific_name_labels=["Homo sapiens"])
        with patch.object(m, "original_files", side_effect=RuntimeError("BAD_LOCATION_ACCESS")) as fetch:
            report = m.inspect(human)
            fetch.assert_not_called()
        self.assertEqual(report["source_status"], "HUMAN_LOCATION_PRIVACY_HOLD")
        self.assertFalse(report["provider_original_file_list_obtained"])
        self.assertEqual(report["animal_GPS_events_opened"], 0)

    def test_over_budget_missing_checksum_and_provider_failure_are_HOLD(self):
        s = self.selected[0]
        with patch.object(m, "original_files", return_value=[original_file(size=200_000_001)]):
            report = m.inspect(s)
        self.assertEqual(report["source_status"], "SOURCE_METADATA_HOLD")
        self.assertEqual(report["source_hold_reason"],
                         "FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES")
        with patch.object(m, "original_files", return_value=[original_file(md5="wrong")]):
            report = m.inspect(s)
        self.assertEqual(report["source_status"], "SOURCE_METADATA_HOLD")
        self.assertEqual(report["source_hold_reason"],
                         "SOURCE_PROVIDER_MD5_OR_BITSTREAM_UUID_MISSING")
        with patch.object(m, "original_files", side_effect=OSError("sensitive internal URL")):
            report = m.inspect(s)
        self.assertEqual(report["source_hold_reason"], "PROVIDER_FILE_METADATA_ERROR_OSError")
        self.assertNotIn("sensitive internal URL", str(report))

    def test_fixed_roster_complete_and_no_GPS_download(self):
        with patch.object(m, "original_files", return_value=[original_file()]) as get:
            d = m.evaluate()
        self.assertEqual(get.call_count, 16)
        self.assertEqual(d["source_study_slots"], 16)
        self.assertEqual(d["new_wave3_single_native"], 16)
        self.assertEqual(d["total_selected_original_source_studies_after_wave3"], 51)
        self.assertEqual(len(d["source_rows"]), 16)
        self.assertEqual(d["source_status_counts"], {
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE": 16
        })
        self.assertEqual([r["original_repository_study_id"] for r in d["source_rows"]],
                         [r["original_repository_study_id"] for r in self.selected])
        self.assertTrue(d["provider_file_metadata_only"])
        self.assertEqual(d["all_new_wave3_animal_GPS_events_opened"], 0)
        self.assertFalse(d["human_GPS_or_position_data_retrieved"])
        self.assertEqual(d["original_48_source_radiation_negative_external_replication_UNCHANGED"],
                         m.NEG)
        self.assertTrue(d["old_12_ERA5_original_raw_SHA_HOLD_UNCHANGED"])


if __name__ == "__main__":
    unittest.main()
