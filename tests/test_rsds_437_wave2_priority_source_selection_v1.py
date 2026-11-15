"""Second-wave fixed-source selection and absence of before-access event peeking."""
from __future__ import annotations
import copy
import unittest
from unittest.mock import patch
from scripts import rsds_437_wave2_priority_source_selection_v1 as m

class Wave2SourcePriority(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.select()

    def test_original_437_frame_and_unopened_wave2_event(self):
        d=self.result
        self.assertEqual(d["wave2_new_original_source_count"],19)
        self.assertEqual(d["original_candidate_source_frame"],{
            "single_native":109,"multi_native":11})
        self.assertEqual(d["two_wave_combined_counts"],{
            "single_native":24,"multi_native":11})
        self.assertEqual(d["all_previously_unselected_metadata_candidates_remaining"],85)
        self.assertEqual(d["original_wave2_GPS_events_opened"],0)
        self.assertEqual(d["original_wave2_provider_files_opened"],0)
        self.assertTrue(d["old_ERA5_original_12_raw_SHA_HOLD"])
        self.assertFalse(d["old_48_signed_sunlight_external_replicated"])

    def test_all_19_studies_do_not_overlap_prior_16(self):
        _,prior=m.freeze()
        prior_ids={x["original_repository_study_id"] for x in
                   prior["pre_event_selected_sources_in_rank_order"]}
        rows=self.result["all_original_wave2_sources_in_rank_order"]
        self.assertEqual(len(rows),19)
        self.assertEqual(len({x["original_repository_study_id"] for x in rows}),19)
        self.assertEqual(len({x["package_uuid"] for x in rows}),19)
        self.assertFalse(prior_ids & {r["original_repository_study_id"] for r in rows})

    def test_wave2_rank_positions_are_existing_wave1_sha(self):
        rows=self.result["all_original_wave2_sources_in_rank_order"]
        self.assertEqual([z["stratum_rank_1_based"] for z in rows[:16]],
                         list(range(9,25)))
        self.assertEqual([z["stratum_rank_1_based"] for z in rows[16:]],
                         [9,10,11])
        self.assertEqual([r["stratum"] for r in rows[:16]],["single_native"]*16)
        self.assertEqual([r["stratum"] for r in rows[16:]],["multi_native"]*3)
        self.assertTrue(all(len(x["source_lottery_SHA256"])==64 for x in rows))
        self.assertTrue(all(x["original_GPS_animal_events_opened"]==0
                            and x["original_source_file_metadata_opened"] is False
                            for x in rows))

    def test_deterministic_and_contract_changed_blocks(self):
        self.assertEqual(m.select(),self.result)
        with patch.object(m,"BASE_SHA","0000000000000000000000000000000000000000"):
            with self.assertRaisesRegex(ValueError,"quotas/history changed"):
                m.freeze()
        with patch.object(m,"pinned",side_effect=ValueError("source SHA modified")):
            with self.assertRaisesRegex(ValueError,"source SHA modified"):
                m.freeze()

if __name__=="__main__":
    unittest.main()
