"""Strict second-wave provider metadata screen and human privacy HOLD tests."""
from __future__ import annotations
import unittest
from unittest.mock import patch

from scripts import rsds_437_wave2_original_file_metadata_v1 as m

def source(name="sample-gps.csv",size=12000000):
    return {"name":name,"bytes":size,"bitstream_uuid":"known-binary-uuid",
            "checksum":{"checkSumAlgorithm":"MD5","value":"a"*32}}

class Wave2OriginalFileAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.first=m.freeze()

    def test_previous_wave_and_new_identity_gate(self):
        self.assertEqual(len(self.first),19)
        self.assertEqual([x["wave2_index"] for x in self.first],list(range(19)))
        self.assertEqual(sum("Homo sapiens" in x["scientific_name_labels"]
                             for x in self.first),1)
        self.assertEqual(len({x["original_repository_study_id"] for x in self.first}),19)

    def test_human_source_never_queries_provider_content(self):
        human=next(r for r in self.first if "Homo sapiens" in r["scientific_name_labels"])
        with patch.object(m,"original_files",side_effect=RuntimeError("HUMAN_DATA_ACCESSED")) as f:
            result=m.inspect(human)
            f.assert_not_called()
        self.assertEqual(result["source_status"],"HUMAN_LOCATION_PRIVACY_HOLD")
        self.assertFalse(result["provider_original_file_list_obtained"])
        self.assertEqual(result["animal_GPS_events_opened"],0)
        self.assertEqual(result["source_hold_reason"],"NO_HUMAN_TRAJECTORY_SOURCE_ACCESS")

    def test_unambiguous_MD5_source_file_only(self):
        nonhuman=next(r for r in self.first if "Homo sapiens" not in r["scientific_name_labels"])
        with patch.object(m,"original_files",return_value=[source()]) as f:
            r=m.inspect(nonhuman)
            f.assert_called_once_with(nonhuman["package_uuid"])
        self.assertTrue(r["provider_original_file_list_obtained"])
        self.assertEqual(r["source_status"],"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        self.assertEqual(r["original_file"]["checksum"]["value"],"a"*32)
        self.assertEqual(r["animal_GPS_events_opened"],0)

    def test_provider_HOLD_does_not_reselect(self):
        nonhuman=next(r for r in self.first if "Homo sapiens" not in r["scientific_name_labels"])
        with patch.object(m,"original_files",return_value=[source(size=200000001)]):
            r=m.inspect(nonhuman)
        self.assertEqual(r["source_status"],"SOURCE_METADATA_HOLD")
        self.assertEqual(r["source_hold_reason"],
                         "FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES")
        with patch.object(m,"original_files",side_effect=OSError("provider unavailable")):
            r=m.inspect(nonhuman)
        self.assertEqual(r["source_status"],"SOURCE_METADATA_HOLD")
        self.assertEqual(r["source_hold_reason"],"PROVIDER_FILE_METADATA_ERROR_OSError")
        self.assertNotIn("provider unavailable",str(r))
        self.assertEqual(r["animal_GPS_events_opened"],0)

    def test_all_nineteen_statuses_include_one_human_HOLD(self):
        with patch.object(m,"original_files",return_value=[source()]):
            d=m.evaluate()
        self.assertEqual(d["source_study_slots"],19)
        self.assertEqual(len(d["source_rows"]),19)
        self.assertEqual(d["source_status_counts"].get("HUMAN_LOCATION_PRIVACY_HOLD"),1)
        self.assertEqual(d["source_status_counts"].get(
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"),18)
        self.assertFalse(d["human_GPS_or_position_data_retrieved"])
        self.assertEqual(d["all_new_wave2_animal_GPS_events_opened"],0)
        self.assertTrue(d["old_12_ERA5_original_raw_SHA_HOLD_UNCHANGED"])

if __name__=="__main__":
    unittest.main()
