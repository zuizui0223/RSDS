"""Original six-source-native first-response archive aggregation regression tests.

This test never downloads new original GPS/ERA5 or rewrites native source
effect values. Confirms exact q producer schema versus DAY/NIGHT producer
schema, a previously failed source-only aggregator field-name boundary.
"""
from __future__ import annotations

from copy import deepcopy
import unittest

from scripts import aggregate_rsds_two_source_six_native_month_heat_v1 as a
from scripts import rsds_two_source_all_native_ungulate_month_heat_v1 as source


def e_q(mean, n=6):
    return {
        "status":"ESTIMABLE",
        "original_tag_n":n,
        "estimate":mean,
        "tag_bootstrap95":[mean-.0001,mean+.0001],
        "nominal_CI_sign":"POSITIVE" if mean>.0001 else "NEGATIVE" if mean<-.0001 else "INCLUDES_ZERO",
    }


def e_part(mean,n=6):
    return {
        "status":"ESTIMABLE",
        "original_tag_n":n,
        "estimate":mean,
        "tag_bootstrap95":[mean-.0001,mean+.0001],
        "nominal_sign":"POSITIVE" if mean>.0001 else "NEGATIVE" if mean<-.0001 else "INCLUDES_ZERO",
    }


def fake_full_source(s,old,index):
    records=[]
    for taxon in s["taxa"]:
        baseline=old[s["doi"]][taxon]
        day=e_part(.0004)
        night=e_part(.0010)
        q=e_q(.0006)
        records.append({
            "original_native_taxon":taxon,
            "source_native_status":"SOURCE_NATIVE_MONTH_COMPONENTS_ESTIMABLE",
            "first_full_original_q":baseline["q"]["mean"],
            "first_full_original_q_tagged_animals":baseline["q"]["source_tags"],
            "first_original_daylength_slope":baseline["daylength"]["mean"],
            "first_original_daylength_tagged_animals":baseline["daylength"]["source_tags"],
            "month_joint_original_tag_n":6,
            "month_native_q_Tmax_given_daylength_antecedent_water":q,
            "month_native_DAY_Tmax_given_daylength_antecedent_water":day,
            "month_native_NIGHT_Tmax_given_daylength_antecedent_water":night,
            "original_native_q_NIGHT_minus_DAY_model_and_bootstrap_algebra_verified":True,
            "first_source_native_q_Tmax_new_ERA5_tag_bootstrap_replayed":True,
        })
    return {
        "version":source.SOURCE_VERSION,
        "source_index":index,
        "original_study_DOI":s["doi"],
        "original_repository_study_id":s["id"],
        "registered_native_taxa":s["taxa"],
        "source_status":"SOURCE_ALL_NATIVE_ERA5_MONTH_MODELS_EVALUATED",
        "source_ERA5_grid_count":15,
        "source_ERA5_grid_HOLDs":0,
        "source_ERA5_daily_physical_SHA256":"a"*64,
        "source_original_GPS_MD5_verified":True,
        "all_registered_native_q_and_solar_group_first_results_replayed":True,
        "original_q_NIGHT_minus_DAY_per_tag_date_reconstructed":True,
        "new_contemporaneous_source_native_ERA5_full_grid_coverage":True,
        "source_native_group_results":records,
        "original_source_mule_and_bighorn_shared_source_grid_yearmonth_count":4,
        "original_source_mule_and_bighorn_shared_grid_yearmonth_with_3_tags_each":2,
        "source_new_weather_is_not_previous_20261023_or_20261116_physical_replay":True,
        "no_original_animal_ID_tracked_coordinates_exact_date_GPS_weather_grid_or_individual_effect_exported":True,
        "previous_48_source_irradiance_heading_direction_NOT_REPLICATED":True,
        "previous_12_source_ERA5_raw_SHA_HOLD":True,
    }


class TestOriginalSourceOnlyArchiveRecovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.studies,cls.history=source.frozen()
        cls.data=[fake_full_source(s,cls.history,i) for i,s in enumerate(cls.studies)]

    def test_exact_q_and_component_producer_sign_keys(self):
        self.assertEqual(a.validated(e_q(.0006), "q"),"POSITIVE")
        self.assertEqual(a.validated(e_part(.0004), "day"),"POSITIVE")
        self.assertEqual(a.validated(e_part(-.0004), "night"),"NEGATIVE")
        self.assertEqual(a.validated(e_part(0), "night"),"INCLUDES_ZERO")

    def test_unexpected_or_ambiguous_keys_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"schema"):
            a.validated(e_q(.0006), "night")
        with self.assertRaisesRegex(ValueError,"schema"):
            a.validated(e_part(.0004), "q")
        both={**e_q(.0006),"nominal_sign":"POSITIVE"}
        with self.assertRaisesRegex(ValueError,"schema"):
            a.validated(both,"q")
        with self.assertRaisesRegex(ValueError,"class"):
            a.validated({**e_q(.0006),"nominal_CI_sign":"NEGATIVE"},"q")

    def test_source_q_missing_is_unavailable_not_zero(self):
        null={"status":"INSUFFICIENT_ORIGINAL_TAGS","original_tag_n":4,
              "estimate":None,"tag_bootstrap95":None,"nominal_CI_sign":"UNAVAILABLE"}
        self.assertEqual(a.validated(null,"q"),"UNAVAILABLE")
        with self.assertRaisesRegex(ValueError,"zero"):
            a.validated({**null,"estimate":0},"q")
        self.assertEqual(a.validated(None,"day"),"UNAVAILABLE")

    def test_all_six_natives_are_study_grain_not_six_original_DOIs(self):
        out=a.aggregate(deepcopy(self.data))
        self.assertEqual(out["original_source_studies"],2)
        self.assertEqual(out["original_native_taxon_study_rows"],6)
        self.assertEqual(out["original_source_complete_tag_month_native_rows"],6)
        self.assertEqual(out["source_native_q_Tmax_sign_counts"],{"POSITIVE":6})
        self.assertEqual(out["source_native_night_Tmax_sign_counts"],{"POSITIVE":6})
        self.assertEqual(out["source_native_day_Tmax_sign_counts"],{"POSITIVE":6})
        self.assertEqual(len(out["all_original_six_source_native_rows_including_HOLD"]),6)

    def test_baseline_replay_identity_mismatch_is_a_HARD_ERROR(self):
        rows=deepcopy(self.data)
        rows[0]["source_native_group_results"][0]["first_full_original_q"]+=.00001
        with self.assertRaisesRegex(ValueError,"replay"):
            a.aggregate(rows)

    def test_original_day_night_q_coefficient_algebra_mismatch_is_a_HARD_ERROR(self):
        rows=deepcopy(self.data)
        rows[1]["source_native_group_results"][0]["month_native_NIGHT_Tmax_given_daylength_antecedent_water"]["estimate"]+=.0005
        with self.assertRaisesRegex(ValueError,"algebra"):
            a.aggregate(rows)

    def test_full_source_HOLD_not_recast_as_negative(self):
        rows=deepcopy(self.data)
        rows[0]["source_status"]="SOURCE_HOLD"
        rows[0]["source_native_group_results"]=[]
        result=a.aggregate(rows)
        self.assertEqual(result["original_native_taxon_study_rows"],6)
        self.assertEqual(result["source_native_q_Tmax_sign_counts"],{"POSITIVE":2,"UNAVAILABLE":4})
        self.assertTrue(result["prior_12_ERA5_original_raw_SHA_HOLD"])


if __name__=="__main__":
    unittest.main()
