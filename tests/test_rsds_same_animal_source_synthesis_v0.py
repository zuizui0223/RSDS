"""Original two-region same-animal ERA5 research result synthesis tests."""
from __future__ import annotations
import json
import unittest
from copy import deepcopy
from pathlib import Path

from scripts.analyze_rsds_same_animal_source_results_v0 import (
    FROZEN_NEGATIVE, ORIGINAL, PREVIOUS, RADIATION,
    analyze, check_complete_original_sources, signed_source
)


class OriginalSameIndividualSourceSynthesis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.first=json.loads(ORIGINAL.read_text(encoding="utf-8"))
        cls.previous=json.loads(PREVIOUS.read_text(encoding="utf-8"))
        cls.radiation=json.loads(RADIATION.read_text(encoding="utf-8"))
        cls.actual=analyze(cls.first,cls.previous,cls.radiation)

    def test_full_two_source_6_original_native_taxa_always_reported(self):
        self.assertEqual(self.actual["original_studies"],2)
        self.assertEqual(self.actual["native_taxon_study_rows"],6)
        self.assertEqual(len(self.actual["all_taxon_study_results"]),6)
        self.assertEqual(self.actual["first_original_full_source_state"],
                         {"Utah":"SAME_ORIGINAL_ANIMAL_WEATHER_FRAME_COMPLETE",
                          "Wyoming":"SAME_ORIGINAL_ANIMAL_WEATHER_FRAME_COMPLETE"})
        self.assertEqual(sum(x["source_original_q_replayed"] for x in self.actual["all_taxon_study_results"]),6)

    def test_observed_original_Utah_bighorn_direction_positive_and_q_unresolved(self):
        z=next(z for z in self.actual["all_taxon_study_results"] if
               z["region"]=="Utah" and z["taxon"]=="Ovis canadensis")
        c=z["measured_source_axis_coefficients"]["sun_cos_per100W"]
        self.assertEqual(c["eligible_original_animals"],117)
        self.assertAlmostEqual(c["source_mean"],0.01686335143164922,places=12)
        self.assertEqual(c["sign_support"],"POSITIVE_SOURCE_BOOTSTRAP")
        self.assertEqual(z["measured_source_axis_coefficients"]["sun_sin_per100W"]["sign_support"],
                         "POSITIVE_SOURCE_BOOTSTRAP")
        self.assertEqual(z["measured_source_axis_coefficients"]["q_daypart_per100W"]["sign_support"],
                         "UNRESOLVED_SIGN")

    def test_all_six_source_daypart_coefficients_have_descriptive_intervals_over_zero(self):
        for z in self.actual["all_taxon_study_results"]:
            q=z["measured_source_axis_coefficients"]["q_daypart_per100W"]
            self.assertEqual(q["sign_support"],"UNRESOLVED_SIGN")
            self.assertTrue(q["animal_bootstrap95"][0] <= 0 <= q["animal_bootstrap95"][1])

    def test_wyoming_bighorn_contextual_tangential_coefficients_reverse(self):
        z=next(z for z in self.actual["all_taxon_study_results"] if
               z["region"]=="Wyoming" and z["taxon"]=="Ovis canadensis")
        self.assertTrue(z["context_sign_sensitivity"]["sun_sin_per100W"]
                        ["both_directions_present_in_selected_contexts"])
        self.assertEqual(z["context_sign_sensitivity"]["sun_sin_per100W"]
                         ["positive_excluding_zero"],1)
        self.assertEqual(z["context_sign_sensitivity"]["sun_sin_per100W"]
                         ["negative_excluding_zero"],1)
        self.assertEqual(z["measured_source_axis_coefficients"]["sun_sin_per100W"]
                         ["sign_support"],"UNRESOLVED_SIGN")

    def test_original_mule_deer_direction_not_globally_replicated(self):
        comp={x["taxon"]:x for x in self.actual["two_study_native_taxon_comparisons"]}
        self.assertEqual(set(comp),{"Odocoileus hemionus","Ovis canadensis"})
        mule=comp["Odocoileus hemionus"]["weather_time_direction"]
        self.assertEqual(mule["sun_cos_per100W"]["Utah"]["sign_support"],"UNRESOLVED_SIGN")
        self.assertEqual(mule["sun_cos_per100W"]["Wyoming"]["sign_support"],"UNRESOLVED_SIGN")

    def test_existing_unreplicated_external_radiation_verdict_is_immutable(self):
        self.assertEqual(self.actual["frozen_external_radiation_nonreplication_UNCHANGED"],FROZEN_NEGATIVE)
        mutant=deepcopy(self.first)
        mutant["frozen_external_radiation_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            analyze(mutant,self.previous,self.radiation)

    def test_unknown_origin_or_omitted_source_taxon_is_fail_closed(self):
        mutant=deepcopy(self.first)
        mutant["sources"][0]["taxon_results"].pop()
        with self.assertRaises(ValueError):
            check_complete_original_sources(mutant,self.previous)
        mutant=deepcopy(self.first)
        mutant["sources"][0]["original_native_taxa"][0]="Changed native taxon"
        with self.assertRaises(ValueError):
            check_complete_original_sources(mutant,self.previous)

    def test_censored_unestimated_values_are_not_ecological_zeros(self):
        self.assertEqual(signed_source(None),"UNAVAILABLE")
        self.assertEqual(signed_source({"source_mean":None}),"UNAVAILABLE")
        with self.assertRaises(ValueError):
            signed_source({"source_mean":2,"source_animal_bootstrap95":[-1,1]})

    def test_public_ecology_output_contains_no_raw_telemetry_identifiers(self):
        rendered=json.dumps(self.actual)
        for field in ('"animal_id"', '"timestamp"', '"longitude"', '"latitude"',
                      '"source_gridkey"', '"individual_id"', '"gps_fix"'):
            self.assertNotIn(field,rendered)


if __name__=="__main__":
    unittest.main()
