"""Known-truth tests before any of six new original GPS event bodies are opened."""
import unittest
from unittest.mock import patch

from scripts.preflight_rsds_new_six_original_gps_v0 import (
    doi, one_source, validate_roster
)


class SixFrozenPublicSourceMetadataTests(unittest.TestCase):
    def test_six_new_study_uuid_doi_taxa_and_absence_from_older_studies(self):
        spec, registry = validate_roster()
        selected = spec["original_fixed_sources"]
        self.assertEqual(len(selected), 6)
        self.assertEqual(len({x["id"] for x in selected}), 6)
        self.assertEqual(len({x["study_id"] for x in selected}), 6)
        self.assertTrue(all(x["reported_animals"] >= 10 for x in selected))
        self.assertTrue(all(x["reported_fixes"] >= 50_000 for x in selected))
        self.assertEqual({x["taxon"] for x in selected},
                         {"Meleagris gallopavo", "Lynx rufus",
                          "Vulpes vulpes","Prionailurus viverrinus",
                          "Ardea alba","Phascolarctos cinereus"})
        self.assertEqual(len(registry), 437)

    def test_full_original_checksummed_single_file_is_eligible(self):
        c, registry = validate_roster()
        s = c["original_fixed_sources"][0]
        f = {"bitstream_uuid":"original-record","name":"original_gps.csv",
             "bytes":42_000,"checksum":{"checkSumAlgorithm":"MD5",
                                         "value":"a"*32}}
        with patch("scripts.preflight_rsds_new_six_original_gps_v0.original_files",
                   return_value=[f]):
            r=one_source(s,registry[s["id"]])
        self.assertEqual(r["source_metadata_status"], "ORIGINAL_FILE_METADATA_ADMISSIBLE")
        self.assertEqual(r["animal_event_rows_opened_by_this_preflight"],0)
        self.assertEqual(r["original_file"]["bitstream_uuid"],"original-record")

    def test_ambiguous_original_source_is_not_rescued(self):
        c, registry = validate_roster()
        s = c["original_fixed_sources"][0]
        f = {"bitstream_uuid":"x","name":"original.csv","bytes":10,
             "checksum":{"checkSumAlgorithm":"MD5","value":"b"*32}}
        with patch("scripts.preflight_rsds_new_six_original_gps_v0.original_files",
                   return_value=[f,f]):
            r=one_source(s,registry[s["id"]])
        self.assertEqual(r["source_metadata_status"],"SOURCE_METADATA_HOLD")
        self.assertEqual(r["original_data_file_candidates"],2)
        self.assertIsNone(r["original_file"])

    def test_gps_file_priority_excludes_original_accelerometer_companions(self):
        c, registry = validate_roster()
        s = c["original_fixed_sources"][0]
        gps={"bitstream_uuid":"gps-only","name":"Bobcat tracking-gps.csv",
             "bytes":50_000,"checksum":{"checkSumAlgorithm":"MD5","value":"c"*32}}
        acc={"bitstream_uuid":"acc-data","name":"Bobcat tracking-acc1of9.csv.zip",
             "bytes":75_000,"checksum":{"checkSumAlgorithm":"MD5","value":"d"*32}}
        with patch("scripts.preflight_rsds_new_six_original_gps_v0.original_files",
                   return_value=[acc,gps]):
            got=one_source(s,registry[s["id"]])
        self.assertEqual(got["source_metadata_status"],"ORIGINAL_FILE_METADATA_ADMISSIBLE")
        self.assertEqual(got["original_file"]["bitstream_uuid"],"gps-only")

    def test_too_large_or_missing_checksum_does_not_trigger_download(self):
        c, registry=validate_roster()
        s=c["original_fixed_sources"][0]
        f={"bitstream_uuid":"x","name":"source.csv",
           "bytes":200_000_001,"checksum":{"checkSumAlgorithm":"MD5","value":"b"*32}}
        with patch("scripts.preflight_rsds_new_six_original_gps_v0.original_files",
                   return_value=[f]):
            r=one_source(s,registry[s["id"]])
        self.assertEqual(r["source_metadata_status"],"SOURCE_METADATA_HOLD")
        f["bytes"]=1000
        f["checksum"]=None
        with patch("scripts.preflight_rsds_new_six_original_gps_v0.original_files",
                   return_value=[f]):
            r=one_source(s,registry[s["id"]])
        self.assertEqual(r["source_metadata_status"],"SOURCE_METADATA_HOLD")

    def test_doi_url_identifier_normalization(self):
        self.assertEqual(doi("doi:10.5441/001/1.703"),"10.5441/001/1.703")
        self.assertEqual(doi("https://doi.org/10.5441/001/1.703"),
                         "10.5441/001/1.703")


if __name__=="__main__":
    unittest.main()
