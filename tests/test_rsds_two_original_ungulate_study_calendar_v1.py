"""Two ORIGINAL source DOI shared calendar native-taxon contrast known-truth tests."""
from __future__ import annotations
import copy
import json
import unittest
from pathlib import Path

from scripts.compare_rsds_two_ungulate_study_shared_date_v1 import (
    compare,validate_study,component,SOURCE,SCHEDULES,NEGATIVE,WYOMING,ROOT
)


def wyoming_fixture():
    return json.loads((ROOT/WYOMING).read_text(encoding="utf-8"))


def utah_synthetic_from_public_schema():
    obj=copy.deepcopy(wyoming_fixture())
    obj["study_Doi"]="10.5441/001/1.711"
    return obj


class TwoIndependentOriginalStudyCalendarComparison(unittest.TestCase):
    def test_wyoming_source_first_result_and_species_contract(self):
        z=wyoming_fixture()
        self.assertEqual(validate_study(z,"10.5441/001/1.703","Wyoming"),
                         "SHARED_DATE_SOURCE_EVALUATED")
        self.assertEqual(z["native_taxa"],list(SOURCE))
        self.assertEqual(z["shared_support"]["common_daylight_standard_reference_latitude_degrees"],43)
        for sch in SCHEDULES:
            r=component(z,sch)
            self.assertTrue(r["positive_mule_negative_bighorn_supported"])
            self.assertEqual(r["contrast"]["nominal_sign"],"NOMINAL_POSITIVE")

    def test_two_study_source_count_not_four_independent_taxa(self):
        utah=utah_synthetic_from_public_schema()
        wyoming=wyoming_fixture()
        out=compare(utah,wyoming)
        self.assertEqual(out["fixed_original_study_DOIs"],2)
        self.assertEqual(out["source_original_studies_full_calendar_comparison_identifiable"],2)
        self.assertEqual(len(out["all_study_source_original_native_photoperiod_sign_results"]),2)
        self.assertEqual(out["original_studies_with_same_opposite_native_daylength_sign_through_three_clocks"],
                         ["Wyoming","Utah"])
        self.assertTrue(out["both_original_studies_show_same_source_native_sign_if_supported"])

    def test_utah_insufficient_source_days_are_neither_confirmed_nor_negative_biology(self):
        utah=utah_synthetic_from_public_schema()
        utah["status"]="INSUFFICIENT_COOCCURRENT_SOURCE_DATE_SUPPORT"
        utah["source_taxon_results"]=[]
        utah["same_study_between_native_taxon"]=None
        utah["shared_support"].update({
            "original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules":0,
            "original_tags_with_10_days_on_both_taxa_shared_calendar":{
                SOURCE[0]:0,SOURCE[1]:0},
        })
        result=compare(utah,wyoming_fixture())
        self.assertEqual(result["source_original_studies_full_calendar_comparison_identifiable"],1)
        self.assertEqual(result["original_studies_with_same_opposite_native_daylength_sign_through_three_clocks"],
                         ["Wyoming"])
        self.assertFalse(result["both_original_studies_show_same_source_native_sign_if_supported"])
        self.assertTrue(all(
            x["status"]=="UNAVAILABLE" for x in
            result["all_study_source_original_native_photoperiod_sign_results"][1][
                "three_observation_schedules"].values()
        ))

    def test_bad_utah_original_MD5_or_source_replay_is_not_a_natural_null(self):
        utah=utah_synthetic_from_public_schema()
        utah["original_provider_GPS_MD5_verified"]=False
        r=compare(utah,wyoming_fixture())
        self.assertEqual(r["all_study_source_original_native_photoperiod_sign_results"][1][
            "shared_calendar_source_status"],"ORIGINAL_PROVIDER_OR_SOURCE_REPLAY_HOLD")
        self.assertFalse(r["both_original_studies_show_same_source_native_sign_if_supported"])

    def test_same_reference_daylength_and_fixed_taxa_fail_closed(self):
        utah=utah_synthetic_from_public_schema()
        utah["shared_support"]["common_daylight_standard_reference_latitude_degrees"]=40.0
        with self.assertRaisesRegex(ValueError,"reference daylength"):
            compare(utah,wyoming_fixture())
        utah=utah_synthetic_from_public_schema()
        utah["native_taxa"]=[SOURCE[0]]
        with self.assertRaisesRegex(ValueError,"Wrong original"):
            compare(utah,wyoming_fixture())

    def test_source_study_mule_bighorn_difference_must_reconstruct_the_tag_means(self):
        utah=utah_synthetic_from_public_schema()
        x=utah["same_study_between_native_taxon"][
            "q_daylength_slope_mule_minus_bighorn_by_schedule"]["schedules"]["original"]
        x["mule_minus_bighorn_slope"]+=0.0001
        with self.assertRaisesRegex(ValueError,"algebra"):
            compare(utah,wyoming_fixture())

    def test_prior_old_radiation_external_no_replication_not_modifiable(self):
        utah=utah_synthetic_from_public_schema()
        utah["historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED"]="REPLICATED"
        with self.assertRaisesRegex(ValueError,"reclassified"):
            compare(utah,wyoming_fixture())
        self.assertEqual(NEGATIVE,"NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")


if __name__=="__main__":
    unittest.main()
