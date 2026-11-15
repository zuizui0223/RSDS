"""Pure metadata original-file screening tests for SHA-locked source wave."""
from __future__ import annotations
import copy
import unittest
from unittest.mock import patch

from scripts import rsds_437_wave1_original_file_metadata_v1 as m

def record(name,size=40000000,checksum="f"*32):
    return {"name":name,"bytes":size,"bitstream_uuid":"source-provider-uuid",
            "checksum":{"checkSumAlgorithm":"MD5","value":checksum}}

class PreEventOriginalSourceMetadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest,cls.wave=m.freeze()

    def test_manifest_closed_before_file_metadata(self):
        self.assertEqual(len(self.wave),16)
        self.assertEqual(self.manifest["animal_event_records_accessed"],0)
        self.assertEqual(len({p["original_repository_study_id"] for p in self.wave}),16)
        self.assertEqual(len({p["package_uuid"] for p in self.wave}),16)
        self.assertEqual(sum(p["selection_stratum"]=="multi_native" for p in self.wave),8)

    def test_selection_manifest_failure_if_altered(self):
        with patch.object(m,"SELECTION_BLOB_SHA","0000000000000000000000000000000000000000"):
            with self.assertRaisesRegex(ValueError,"modified"):
                m.freeze()
        with patch.object(m,"FRAME_SHA256","bad"):
            with self.assertRaisesRegex(ValueError,"changed"):
                m.freeze()

    def test_unambiguous_original_tracking_filename(self):
        f,hold,n,rule=m.original_primary([
            record("readme.csv",120),record("telemetry.gps.original.csv",10000),
            record("other metadata.csv",123)])
        self.assertEqual(hold,None)
        self.assertEqual(n,1)
        self.assertEqual(rule,"UNIQUE_GPS_OR_TRACKING_FILENAME_ONLY")
        self.assertEqual(f["bytes"],10000)
        self.assertEqual(f["checksum"]["algorithm"],"MD5")

    def test_exact_10x_relative_file_size_rule(self):
        file,hold,n,rule=m.original_primary([
            record("raw.csv",20000000),record("small.csv",2000000)])
        self.assertIsNone(hold)
        self.assertEqual(file["bytes"],20000000)
        self.assertEqual(rule,"LARGEST_CSV_10X_SMALLER_COMPANION_ONLY")
        file,hold,n,rule=m.original_primary([
            record("raw.csv",18000000),record("small.csv",2000000)])
        self.assertIsNone(file)
        self.assertEqual(hold,"AMBIGUOUS_OR_MISSING_PRIMARY_SOURCE_FILE")

    def test_hold_never_resamples_better_file(self):
        for files,reason in [
            ([record("a.csv"),record("b.zip")],"AMBIGUOUS_OR_MISSING_PRIMARY_SOURCE_FILE"),
            ([record("a.zip",200000001)],"FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES"),
            ([record("a.csv",40000000,checksum="")],"SOURCE_PROVIDER_MD5_OR_BITSTREAM_UUID_MISSING"),
            ([record("a.csv") ,record("b.csv")],"AMBIGUOUS_OR_MISSING_PRIMARY_SOURCE_FILE"),
            ([], "AMBIGUOUS_OR_MISSING_PRIMARY_SOURCE_FILE"),
        ]:
            z,status,n,rule=m.original_primary(files)
            self.assertIsNone(z)
            self.assertEqual(status,reason)

    def test_stage2_provider_api_uses_only_preselected_package_and_hold(self):
        wave=self.wave[0]
        with patch.object(m,"original_files",return_value=[record("GPS-original.csv")]) as ref:
            result=m.inspect(wave)
            self.assertEqual(result["source_status"],"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
            self.assertEqual(result["original_file"]["bytes"],40000000)
            self.assertEqual(result["animal_event_records_opened"],0)
            ref.assert_called_once_with(wave["package_uuid"])
        with patch.object(m,"original_files",side_effect=RuntimeError("provider unavailable")):
            result=m.inspect(wave)
            self.assertEqual(result["source_status"],"SOURCE_METADATA_HOLD")
            self.assertEqual(result["hold_reason"],"SOURCE_PROVIDER_FILE_METADATA_ERROR_RuntimeError")
            self.assertNotIn("provider unavailable",str(result))

    def test_metadata_output_remains_anonymized_and_complete(self):
        with patch.object(m,"original_files",return_value=[record("study-GPS.csv")]):
            r=m.evaluate()
        self.assertEqual(r["original_package_slots"],16)
        self.assertEqual(r["selected_single_native_slots"],8)
        self.assertEqual(r["selected_multi_native_slots"],8)
        self.assertEqual(r["original_GPS_event_records_opened"],0)
        self.assertEqual(r["source_status_counts"]["ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"],16)
        self.assertNotIn("animal_id",str(r))
        self.assertNotIn("latitude",str(r))
        self.assertTrue(r["prior_48_shortwave_direction_external_not_replicated"])

if __name__=="__main__":
    unittest.main()
