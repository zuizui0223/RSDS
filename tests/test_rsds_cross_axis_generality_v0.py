"""Pinned original field-study temporal × solar-heading generality checks."""
import json
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.analyze_rsds_cross_axis_generality_v0 import (
    DECISION_FROZEN, HEADING, NATIVE_TAXA, RAD_A, RAD_B, RAD_DECISION, ROOT,
    SEASON_A, SEASON_B, TIME, csv_read, json_read, normalize_doi,
    choose_six_originals, pair_time_and_direction, radiation_time_overlap,
    nine_taxon_repeatability, support_sign, run,
)


class CrossAxisOriginalEcologyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()

    def test_six_original_source_taxa_link_by_same_public_DOI_not_label_guess(self):
        rows = self.result["joint_time_direction_source_taxon_atlas"]
        self.assertEqual(len(rows), 6)
        self.assertEqual(len({(r["source_doi"], r["taxon"]) for r in rows}), 6)
        self.assertEqual(len({r["source_doi"] for r in rows}), 2)
        self.assertEqual(set(r["taxon"] for r in rows), set(NATIVE_TAXA))
        self.assertEqual({r["source_doi"] for r in rows},
                         {"10.5441/001/1.711", "10.5441/001/1.703"})

    def test_mule_deer_same_signed_temporal_assoc_but_opposing_solar_raw(self):
        rows = {r["source_region"]: r for r in self.result["joint_time_direction_source_taxon_atlas"]
                if r["taxon"] == "Odocoileus hemionus"}
        self.assertAlmostEqual(rows["Utah"]["daylength_slope_per_hour"], 0.005619415769068772, places=14)
        self.assertAlmostEqual(rows["Wyoming"]["daylength_slope_per_hour"], 0.005538249061142167, places=14)
        self.assertEqual(rows["Utah"]["daylength_sign"], "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP")
        self.assertEqual(rows["Wyoming"]["daylength_sign"], "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP")
        self.assertAlmostEqual(rows["Utah"]["mean_sunward_cosine"], 0.023393157844424365, places=14)
        self.assertAlmostEqual(rows["Wyoming"]["mean_sunward_cosine"], -0.04875792375635666, places=14)
        self.assertEqual(rows["Utah"]["solar_raw_sign"], "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP")
        self.assertEqual(rows["Wyoming"]["solar_raw_sign"], "NEGATIVE_SOURCE_ANIMAL_BOOTSTRAP")
        q = next(x for x in self.result["independent_study_cross_regional_taxon_pairs"]
                 if x["taxon"] == "Odocoileus hemionus")
        self.assertTrue(q["both_positive_daylength_CIs"])
        self.assertTrue(q["both_raw_direction_CIs_have_opposite_signs"])
        self.assertAlmostEqual(q["raw_sunward_mean_difference_Utah_minus_Wyoming"],
                               q["solar_marginal_component_difference"] +
                               q["solar_joint_component_difference"], places=12)

    def test_bighorn_has_unresolved_wyoming_direction_and_daylength(self):
        p = next(x for x in self.result["independent_study_cross_regional_taxon_pairs"]
                 if x["taxon"] == "Ovis canadensis")
        self.assertFalse(p["both_positive_daylength_CIs"])
        self.assertFalse(p["both_raw_direction_CIs_have_opposite_signs"])
        rows = {r["source_region"]: r for r in self.result["joint_time_direction_source_taxon_atlas"]
                if r["taxon"] == "Ovis canadensis"}
        self.assertEqual(rows["Utah"]["daylength_sign"], "UNRESOLVED_SIGN")
        self.assertEqual(rows["Wyoming"]["daylength_sign"], "UNRESOLVED_SIGN")
        self.assertEqual(rows["Utah"]["solar_raw_sign"], "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP")
        self.assertEqual(rows["Wyoming"]["solar_raw_sign"], "UNRESOLVED_SIGN")

    def test_nightday_species_ranking_repeats_but_not_sympatric(self):
        p = self.result["independent_study_cross_regional_taxon_pairs"][0]
        self.assertGreater(p["night_day_q_mule_minus_bighorn_Utah"], 0)
        self.assertGreater(p["night_day_q_mule_minus_bighorn_Wyoming"], 0)
        self.assertAlmostEqual(p["night_day_q_mule_minus_bighorn_Utah"], 0.043911916, places=7)
        self.assertAlmostEqual(p["night_day_q_mule_minus_bighorn_Wyoming"], 0.015512570, places=7)

    def test_repeated_taxon_directed_circular_covariance_is_not_universal(self):
        s = self.result["nine_taxon_solar_heading_repeatability"]
        self.assertEqual(s["n_taxa_with_exactly_two_archived_sources"], 9)
        self.assertEqual(s["same_raw_cosine_sign"], 6)
        self.assertEqual(s["opposite_raw_cosine_sign"], 3)
        self.assertEqual(s["same_joint_circular_covariance_sign"], 4)
        self.assertEqual(s["opposite_joint_circular_covariance_sign"], 5)

    def test_radiation_sample_does_not_invent_cross_axis_individual_study(self):
        a = self.result["sunlight_weather_temporal_observation_overlap"]
        self.assertEqual(a["radiation_original_source_archives"], 48)
        self.assertEqual(a["radiation_estimable_shortwave_slopes"], 20)
        self.assertEqual(a["radiation_undefined_source_slopes_NOT_ZERO"], 28)
        self.assertEqual(a["shared_DOIs_with_diel_daylength_studies"], 0)
        self.assertEqual(a["shared_DOIs_with_estimable_seasonal_branch_effect"], 0)
        self.assertEqual(a["preregistered_heldout_decision_UNCHANGED"], DECISION_FROZEN)

    def test_duplicate_study_taxon_or_undeclared_input_fails_closed(self):
        t = csv_read(ROOT, TIME)
        h = csv_read(ROOT, HEADING)
        duplicated = t + [t[0].copy()]
        with self.assertRaises(ValueError):
            choose_six_originals(duplicated, h)
        changed = deepcopy(h)
        changed[5]["joint_direction_solar_term"] = "-0.9"
        with self.assertRaises(ValueError):
            choose_six_originals(t, changed)

    def test_ineligible_source_bootstrap_sign_not_overclaimed(self):
        self.assertEqual(support_sign(-0.01, -0.02, 0.01), "UNRESOLVED_SIGN")
        self.assertEqual(support_sign(None, None, None), "UNAVAILABLE")
        with self.assertRaises(ValueError):
            support_sign(0.1, -0.2, -0.1)

    def test_src_file_hashes_are_original_and_output_contains_no_tracking(self):
        self.assertEqual(len(self.result["source_provenance_git_blob_sha"]), 8)
        rendered = json.dumps(self.result)
        for name in ('"animal_id"', '"individual_id"', '"latitude"', '"longitude"',
                     '"timestamp"', '"location_lat"', '"location_long"', '"gps_fix"'):
            self.assertNotIn(name, rendered)
        self.assertNotIn("NaN", rendered)
        self.assertEqual(self.result["status"],
                         "POST_OUTCOME_ORIGINAL_SOURCE_OBSERVATIONAL_CROSS_AXIS_ATLAS")

    def test_radiation_nonreplication_cannot_be_reclassified_by_adapter(self):
        d = json_read(ROOT, RAD_DECISION)
        self.assertEqual(d["prior_frozen_external_replication_decision_UNCHANGED"], DECISION_FROZEN)
        mutant = deepcopy(d)
        mutant["prior_frozen_external_replication_decision_UNCHANGED"] = "REPLICATED"
        with self.assertRaises(ValueError):
            radiation_time_overlap(json_read(ROOT, RAD_A), json_read(ROOT, RAD_B),
                                   mutant, csv_read(ROOT, TIME),
                                   json_read(ROOT, SEASON_A), json_read(ROOT, SEASON_B))


if __name__ == "__main__":
    unittest.main()
