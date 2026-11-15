"""Pinned original DOI-grain q-measurability discriminability test."""
from __future__ import annotations
import unittest
from unittest.mock import patch

from scripts import rsds_35_source_q_measurability_rank_diagnostic_v1 as m


class StudyGrainMeasurability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.run()

    def test_discrimination_uses_original_study_not_fix_count(self):
        r=self.result
        self.assertEqual(r["original_35_source_frame"],35)
        self.assertEqual(r["GPS_parsed_original_sources"],25)
        self.assertEqual(r["GPS_native_source_q_complete"],10)
        self.assertEqual(r["original_replayed_source_native_rows"],34)
        for item in r["original_measurement_scores"].values():
            self.assertEqual(item["source_study_N_with_defined_score"],
                             item["source_unique_study_IDs"])
            self.assertTrue(item[
                "not_validated_out_of_source_prediction_or_selection_correction"])
            self.assertLessEqual(item["source_study_N_with_defined_score"],25)

    def test_first_original_ratios_near_nonseparating_with_effect_complete_q(self):
        r=self.result["original_measurement_scores"]
        self.assertEqual(r["full_original_segment_acceptance_fraction"][
            "source_study_N_with_defined_score"],24)
        self.assertAlmostEqual(r["full_original_segment_acceptance_fraction"][
            "all_parsed_sources"]["AUC_HIGHSCORE_PREDICTS_MEASURABLE"],
            0.5428571428571428,places=12)
        self.assertAlmostEqual(r["time_0p5_to_2h_fraction"][
            "all_parsed_sources"]["AUC_HIGHSCORE_PREDICTS_MEASURABLE"],
            0.7071428571428572,places=12)
        self.assertAlmostEqual(r["time_0p5_to_2h_fraction"][
            "within_original_wave_descriptive"]["wave1"][
            "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],.75,places=12)
        self.assertAlmostEqual(r["time_0p5_to_2h_fraction"][
            "within_original_wave_descriptive"]["wave2"][
            "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],2/3,places=12)

    def test_qualified_tag_effect_is_definition_not_independent_predictor(self):
        r=self.result
        self.assertGreater(r["original_measurement_scores"][
            "qualified_q_tags_definition_aware"]["all_parsed_sources"][
                "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],.9)
        self.assertIn("NOT an independent predictor",r["important_negative_control"])
        self.assertTrue(r["original_48_sunlight_external_negative_nonreplication_UNCHANGED"])
        self.assertTrue(r["original_12_ERA5_raw_JSON_HOLD_UNCHANGED"])

    def test_rank_probability_known_truth_and_unidentifiable(self):
        examples=[
            {"q_complete":True,"score":2.0},
            {"q_complete":True,"score":3.0},
            {"q_complete":False,"score":1.0},
            {"q_complete":False,"score":4.0}
        ]
        obj=m.rank_probability(examples)
        self.assertEqual(obj["pairwise_source_comparisons"],4)
        self.assertEqual(obj["AUC_HIGHSCORE_PREDICTS_MEASURABLE"],.5)
        examples[0]["score"]=4
        examples[3]["score"]=4
        self.assertAlmostEqual(m.rank_probability(examples)[
            "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],.625)
        self.assertIsNone(m.rank_probability(examples[:2])[
            "AUC_HIGHSCORE_PREDICTS_MEASURABLE"])

    def test_frozen_first_Git_blob_modified_fails_closed(self):
        with patch.object(m,"FIRST",(
             "research/RSDS_35_SOURCE_GPS_OBSERVATION_FUNNEL_FIRST_RESULT_20261115.json",
             "0"*40)):
            with self.assertRaisesRegex(ValueError,"Previously fixed source/blob changed"):
                m.run()

if __name__=="__main__":
    unittest.main()
