"""Before-outcome guards for wave2 original GPS and strict human no-access."""
from __future__ import annotations
import unittest
from unittest.mock import patch

from scripts import run_rsds_437_wave2_original_gps_first_effects_v1 as m

class FirstWave2GPS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.sources=m.freeze_original()

    def test_fixed_19_slots_and_4_protected_preflight_holds(self):
        self.assertEqual(len(self.sources),19)
        from collections import Counter
        counts=Counter(s["source_status"] for s in self.sources)
        self.assertEqual(counts,{"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":15,
                                 "SOURCE_METADATA_HOLD":3,
                                 "HUMAN_LOCATION_PRIVACY_HOLD":1})
        self.assertEqual([s["wave2_index"] for s in self.sources],list(range(19)))
        self.assertEqual(sum(len(s["registered_native_taxa"]) for s in self.sources),23)

    def test_human_source_never_reads_original_provider_API_or_events(self):
        with patch.object(m,"original_files",side_effect=RuntimeError("SHOULD_NEVER_RUN")) as p:
            x=m.process(10)
            p.assert_not_called()
        self.assertEqual(x["source_science_status"],"HUMAN_LOCATION_PRIVACY_HOLD")
        self.assertEqual(x["source_hold_stage"],"HUMAN_SOURCE_PRIVACY_NO_ACCESS")
        self.assertEqual(x["original_source_GPS_event_rows_opened"],0)
        self.assertFalse(x["human_mobility_GPS_data_opened"])
        self.assertEqual(len(x["native_group_results"]),1)
        self.assertIsNone(x["native_group_results"][0]["q"])

    def test_file_HOLDs_never_trigger_original_download_or_dropped_taxa(self):
        indices=[i for i,s in enumerate(self.sources) if s["source_status"]=="SOURCE_METADATA_HOLD"]
        self.assertEqual(indices,[1,17,18])
        with patch.object(m,"original_files",side_effect=RuntimeError("NO_PROVIDER")) as p:
            for i in indices:
                x=m.process(i)
                self.assertEqual(x["source_science_status"],"SOURCE_FILE_METADATA_HOLD")
                self.assertEqual(len(x["native_group_results"]),
                                 len(self.sources[i]["registered_native_taxa"]))
                self.assertTrue(all(y["q"] is None for y in x["native_group_results"]))
            p.assert_not_called()

    def test_provider_MD5_error_stops_without_response_based_source_selection(self):
        with patch.object(m,"original_files",side_effect=OSError("no provider")):
            r=m.process(0)
        self.assertEqual(r["source_science_status"],"SOURCE_GPS_HOLD")
        self.assertEqual(r["source_hold_stage"],"ORIGINAL_UUID_BYTES_MD5")
        self.assertEqual(r["source_hold_reason"],"OSError")
        self.assertNotIn("no provider",str(r))
        self.assertFalse(r["original_file_MD5_bytes_verified"])
        self.assertEqual(r["native_group_results"][0]["status"],"SOURCE_WIDE_FAIL_CLOSED")

    def test_out_of_manifest_source_no_entry(self):
        with self.assertRaisesRegex(ValueError,"nineteen"):
            m.process(-1)
        with self.assertRaisesRegex(ValueError,"nineteen"):
            m.process(19)
        with self.assertRaisesRegex(ValueError,"nineteen"):
            m.process(20)

    def test_reclassified_old_external_shortwave_failure_stops(self):
        with patch.object(m,"NEG","REPLICATED"):
            with self.assertRaisesRegex(ValueError,"source/GPS contract altered"):
                m.freeze_original()

if __name__=="__main__":
    unittest.main()
