"""Strict pre-outcome source privacy and fail-closed controls for Wave3."""
from __future__ import annotations

import unittest
from unittest.mock import patch

from scripts import run_rsds_437_wave3_original_gps_first_effects_v1 as m


class OriginalWave3GPSGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract, cls.sources = m.freeze_original()

    def test_frozen_sha_source_ROSTER_11_MD5_candidate_5_HOLD(self):
        self.assertEqual(len(self.sources), 16)
        self.assertEqual(self.contract["original_file_max_bytes"], 200_000_000)
        self.assertEqual(
            [r["original_single_native_source_rank_1_based"] for r in self.sources],
            list(range(25,41)),
        )
        self.assertEqual(sum(x["source_status"] == "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
                             for x in self.sources), 11)
        self.assertEqual(sum(x["source_status"] == "SOURCE_METADATA_HOLD"
                             for x in self.sources), 5)
        self.assertEqual(len({r["original_repository_study_id"] for r in self.sources}),16)
        self.assertTrue(all("Homo sapiens" not in r["registered_native_taxa"]
                            for r in self.sources))

    def test_every_prior_source_HOLD_is_never_downloaded(self):
        indices = [i for i,r in enumerate(self.sources) if r["source_status"]=="SOURCE_METADATA_HOLD"]
        self.assertEqual(indices,[2,4,6,8,11])
        with patch.object(m, "original_files", side_effect=RuntimeError("NO_PROVIDER")) as p:
            for i in indices:
                result = m.process(i)
                self.assertEqual(result["source_science_status"],"SOURCE_FILE_METADATA_HOLD")
                self.assertIsNone(result["native_group_results"][0]["q"])
                self.assertFalse(result["original_file_MD5_bytes_verified"])
                self.assertEqual(result["original_source_GPS_event_rows_opened"],0)
                self.assertEqual(result["original_registered_native_taxa"],
                                 self.sources[i]["registered_native_taxa"])
            p.assert_not_called()

    def test_changed_provider_file_holds_original_source_fail_closed(self):
        i = next(i for i,r in enumerate(self.sources)
                 if r["source_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        with patch.object(m, "original_files", side_effect=OSError("remote private URL")):
            result = m.process(i)
        self.assertEqual(result["source_science_status"], "SOURCE_GPS_HOLD")
        self.assertEqual(result["source_hold_stage"], "ORIGINAL_UUID_BYTES_MD5")
        self.assertEqual(result["source_hold_reason"], "OSError")
        self.assertEqual(result["native_group_results"][0]["status"],"SOURCE_WIDE_FAIL_CLOSED")
        self.assertNotIn("remote private URL",str(result))

    def test_only_16_original_source_slots(self):
        for ix in (-1, 16, 20):
            with self.assertRaisesRegex(ValueError, "sixteen"):
                m.process(ix)

    def test_historical_negative_external_solar_is_not_turned_positive(self):
        with patch.object(m, "NEG", "REPLICATED"):
            with self.assertRaisesRegex(ValueError,"contract altered"):
                m.freeze_original()


if __name__ == "__main__":
    unittest.main()
