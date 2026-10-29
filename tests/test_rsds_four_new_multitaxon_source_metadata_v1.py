"""Original Movebank source metadata ONLY, before new shared-day species q outcomes."""
import unittest
from unittest.mock import patch

from scripts.preflight_rsds_four_new_multitaxon_gps_v1 import (
    source_registry, source_only_one, MAX_BYTES
)


class FourOriginalMultitaxonStudyFileFreezeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract, cls.catalog=source_registry()

    def test_four_distinct_original_field_studies_ten_native_taxa(self):
        v=self.contract["sources"]
        self.assertEqual(len(v),4)
        self.assertEqual(sum(len(x["registered_taxa"]) for x in v),10)
        self.assertEqual(len({x["study_id"] for x in v}),4)
        self.assertEqual(len({x["doi"] for x in v}),4)
        self.assertEqual([x["doi"] for x in v],[
            "10.5441/001/1.2cp86266",
            "10.5441/001/1.3fv21n7m",
            "10.5441/001/1.mf903197",
            "10.5441/001/1.fk899541"
        ])
        self.assertTrue(all(x["reported_fixes"]>=150000 for x in v))
        self.assertTrue(all(len(set(x["registered_taxa"]))>=2 for x in v))

    def test_original_provider_single_checksum_tracking_csv_without_reading_animal_events(self):
        x=self.contract["sources"][0]
        file={"name":"source-original-tracking.csv","bitstream_uuid":"test-uuid",
              "bytes":1050,
              "checksum":{"checkSumAlgorithm":"MD5","value":"f"*32}}
        with patch("scripts.preflight_rsds_four_new_multitaxon_gps_v1.original_files",
                   return_value=[file]):
            result=source_only_one(x,self.catalog[x["package_uuid"]])
        self.assertEqual(result["source_status"],"SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE")
        self.assertEqual(result["source_original_file"]["bitstream_uuid"],"test-uuid")
        self.assertEqual(result["original_animal_events_opened_before_file_freeze"],0)
        self.assertNotIn("animal_id",result)

    def test_companion_accelerometer_excluded_before_tracking_GPS_selection(self):
        x=self.contract["sources"][2]
        files=[{"name":"research-gps.csv","bitstream_uuid":"gps",
                "bytes":5000,"checksum":{"checkSumAlgorithm":"MD5","value":"a"*32}},
               {"name":"research-acc1of5.csv.zip","bitstream_uuid":"acc",
                "bytes":100,"checksum":{"checkSumAlgorithm":"MD5","value":"b"*32}}]
        with patch("scripts.preflight_rsds_four_new_multitaxon_gps_v1.original_files",
                   return_value=files):
            out=source_only_one(x,self.catalog[x["package_uuid"]])
        self.assertEqual(out["source_status"],"SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE")
        self.assertEqual(out["source_original_file"]["bitstream_uuid"],"gps")

    def test_ambiguous_source_does_not_pick_favorable_tracking_file(self):
        x=self.contract["sources"][1]
        files=[{"name":"a.csv","bitstream_uuid":"a","bytes":100,
                "checksum":{"checkSumAlgorithm":"MD5","value":"a"*32}},
               {"name":"b.csv","bitstream_uuid":"b","bytes":100,
                "checksum":{"checkSumAlgorithm":"MD5","value":"b"*32}}]
        with patch("scripts.preflight_rsds_four_new_multitaxon_gps_v1.original_files",
                   return_value=files):
            r=source_only_one(x,self.catalog[x["package_uuid"]])
        self.assertEqual(r["source_status"],"ORIGINAL_SOURCE_METADATA_HOLD")
        self.assertIsNone(r["source_original_file"])

    def test_source_file_cap_fixed_before_behavior_outcome(self):
        x=self.contract["sources"][3]
        file={"name":"source.csv","bitstream_uuid":"primary","bytes":MAX_BYTES+1,
              "checksum":{"checkSumAlgorithm":"MD5","value":"d"*32}}
        with patch("scripts.preflight_rsds_four_new_multitaxon_gps_v1.original_files",
                   return_value=[file]):
            result=source_only_one(x,self.catalog[x["package_uuid"]])
        self.assertEqual(result["source_status"],"ORIGINAL_SOURCE_METADATA_HOLD")
        self.assertEqual(result["hold_reason"],"FROZEN_200MB_ORIGINAL_FILE_BUDGET_OR_MISSING_SIZE")


if __name__=="__main__":
    unittest.main()
