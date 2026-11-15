"""Source-level original-snapshot evidence tests for q-sign fragility."""
from __future__ import annotations

import unittest
from unittest.mock import patch
from scripts import rsds_437_source_measurability_q_sign_fragility_v1 as m

class SourceSamplingAndSignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=m.run()

    def test_all_original_source_gates_and_old_nulls(self):
        s=self.d["source_study_attrition"]
        self.assertEqual(s["selected_original_source_studies"],35)
        self.assertEqual(s["provider_file_or_privacy_HOLD"],7)
        self.assertEqual(s["original_GPS_source_wide_HOLD"],3)
        self.assertEqual(s["original_GPS_parsed_but_no_native_q"],12)
        self.assertEqual(s["source_has_at_least_one_native_q"],13)
        self.assertEqual(s["source_complete_all_registered_native_q"],10)
        self.assertEqual(s["source_complete_all_native_daylength"],8)
        self.assertEqual(s["partial_native_q_sources"],3)
        self.assertEqual(s["multi_native_same_calendar_month_covariance_estimated"],0)
        self.assertEqual(s["remaining_untouched_single_native_metadata_source_candidates"],85)
        self.assertTrue(self.d["old_48_source_external_solar_direction_negative_nonreplicated"])
        self.assertTrue(self.d["old_12_ERA5_raw_json_source_hold"])

    def test_respondent_q_sign_reversal_is_not_global(self):
        x=self.d["response_only_poststratified_q"]
        self.assertAlmostEqual(x["wave1_original_16_source_response_only"],-.0477935014187846,places=12)
        self.assertAlmostEqual(x["wave1_plus_wave2_35_source_response_only"],.0756846129419663,places=12)
        self.assertTrue(x["sign_reversal_after_wave2_sources"])
        self.assertEqual(x["source_n_with_full_q"],10)
        self.assertTrue(x["not_population_estimand_not_confirmatory"])
        self.assertEqual(self.d["source_effect_signs"],
                         {"respondent_full_source_q_positive":7,
                          "respondent_full_source_q_negative":3,
                          "respondent_full_source_q_zero":0})

    def test_leave_one_out_and_two_out_are_source_level_not_native_level(self):
        x=self.d["response_only_poststratified_q"]
        self.assertEqual(len(x["leave_one_source_out_records"]),10)
        self.assertAlmostEqual(x["leave_one_source_out_range"][0],
                               .021550332593264866,places=12)
        self.assertAlmostEqual(x["leave_one_source_out_range"][1],
                               .1070515422845755,places=12)
        self.assertEqual(x["leave_two_sensitivity_count"],44)
        self.assertAlmostEqual(x["leave_two_source_out_range"][0],
                               -.004395193192772982,places=12)
        worst=x["most_negative_two_source_removal"]
        self.assertEqual(set(worst["held_out_source_DOIs"]),{
                         "doi:10.5441/001/1.291","doi:10.5441/001/1.212g53s7"})
        self.assertLess(worst["remaining_respondent_only_poststratified_q"],0)

    def test_missing_effect_range_condition_is_not_ci(self):
        x=self.d["response_only_poststratified_q"]
        self.assertAlmostEqual(x["hypothetical_missing_effect_sensitivity_for_35_selected"][0],
                               -.09501132019904394,places=12)
        self.assertAlmostEqual(x["hypothetical_missing_effect_sensitivity_for_35_selected"][1],
                               .35688378182026836,places=12)
        y=self.d["finite_120_metadata_frame_source_measurability"]["q"]
        self.assertEqual(y["by_source_stratum"],{"single_native":8,"multi_native":2})
        self.assertEqual(y["strict_finite_frame_identified_measurability_bounds"],
                         [10/120,95/120])
        self.assertTrue(y["not_a_world_prevalence_or_estimate_of_causal_measurement"])

    def test_parent_Git_blob_tamper_causes_failure(self):
        with patch.object(m,"FIRST",("research/RSDS_437_WAVE1_PLUS_WAVE2_35_ORIGINAL_GPS_FIRST_META_RESULT_20261115.json","0000000000000000000000000000000000000000")):
            with self.assertRaisesRegex(ValueError,"Previously fixed source/blob changed"):
                m.run()

if __name__=="__main__":
    unittest.main()
