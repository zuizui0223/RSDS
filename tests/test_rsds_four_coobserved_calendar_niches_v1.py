"""Immutable original coobserved GPS calendars: daily niche ranks vs daylight slopes."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import scripts.build_rsds_four_coobserved_calendar_ecological_niches_v1 as m


class FourOriginalSourceStudyCoobservedCalendars(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=m.analyze()

    def test_four_independent_published_studies_ten_native_groups_eight_taxa(self):
        r=self.d
        self.assertEqual(r["original_source_studies_with_identical_native_source_calendar_days"],4)
        self.assertEqual(r["native_study_taxon_rows"],10)
        self.assertEqual(r["distinct_native_taxa_labels"],8)
        self.assertEqual(r["study_calendar_days_in_separate_study_frames_NOT_INDEPENDENT_DAILY_REPLICATES"],681)
        self.assertEqual([x["source_shared_calendar_days"] for x in r["original_study_calendar_q_daylength_comparisons"]],
                         [169,193,199,120])
        self.assertEqual([x["source_year_month_blocks"] for x in r["original_study_calendar_q_daylength_comparisons"]],
                         [9,8,10,5])

    def test_all_four_original_field_studies_have_a_supported_native_q_rank(self):
        r=self.d
        self.assertEqual(r["original_studies_with_at_least_one_resolved_native_relative_q_rank"],4)
        s=r["original_study_calendar_q_daylength_comparisons"]
        self.assertEqual([z["source_native_pairs_with_resolved_observed_relative_q_rank"] for z in s],
                         [3,1,1,2])
        self.assertEqual(len({z["original_study_DOI"] for z in s}),4)
        self.assertEqual([z["native_source_taxon_rows"] for z in s],[3,2,2,3])

    def test_opposing_source_native_daylight_slopes_not_a_global_test(self):
        r=self.d
        self.assertEqual(r["original_studies_with_opposite_supported_native_daylength_slope_signs"],3)
        s=r["original_study_calendar_q_daylength_comparisons"]
        self.assertEqual([z["source_contains_nominally_positive_AND_negative_native_photoperiod_slopes"]
                          for z in s],[True,True,False,True])
        goose=s[3]
        branta=next(x for x in goose["all_source_native_taxon_q_and_season"]
                    if x["taxon"]=="Branta leucopsis")
        anser=next(x for x in goose["all_source_native_taxon_q_and_season"]
                   if x["taxon"]=="Anser anser")
        self.assertEqual(branta["source_reference_daylength_slope_per_hour"]["nominal_sign"],"NEGATIVE")
        self.assertEqual(anser["source_reference_daylength_slope_per_hour"]["nominal_sign"],"POSITIVE")
        self.assertAlmostEqual(branta["source_q"]["mean"],-0.1840335077564767,places=13)

    def test_repeated_same_two_native_taxa_q_rank_two_DOIs_but_daylight_not_replicated(self):
        z=self.d["same_native_mule_deer_minus_bighorn_repeated_pair_across_two_DOIs"]
        self.assertEqual(z["source_studies"],2)
        self.assertEqual(set(z["two_published_source_DOIs"]),
                         {"10.5441/001/1.703","10.5441/001/1.711"})
        self.assertTrue(z["all_three_source_GPS_UTC_schedules_show_positive_q_pair_rank"])
        self.assertTrue(z["opposite_native_daylength_signs_supported_only_in_Wyoming"])
        self.assertTrue(z["Utah_native_daylength_signs_both_positive_and_pair_contrast_unresolved"])
        s=self.d["original_study_calendar_q_daylength_comparisons"]
        for study in s[1:3]:
            p=study["all_within_original_study_native_pair_q_rank_contrasts"][0]
            self.assertEqual(set(p["three_fixed_original_UTC_clock_q_contrast_signs"].values()),
                             {"POSITIVE"})
        self.assertEqual(s[1]["all_source_native_taxon_q_and_season"][0][
            "source_reference_daylength_slope_per_hour"]["nominal_sign"],"POSITIVE")
        self.assertEqual(s[1]["all_source_native_taxon_q_and_season"][1][
            "source_reference_daylength_slope_per_hour"]["nominal_sign"],"NEGATIVE")

    def test_months_are_source_conditional_with_Utah_NOT_TESTED(self):
        mth=self.d["native_original_month_conditional_pair_source_comparisons"]
        self.assertEqual(set(mth),{"Illinois","Wyoming","Geese","Utah"})
        self.assertEqual(len(mth["Illinois"]),3)
        self.assertEqual(len(mth["Wyoming"]),1)
        self.assertEqual(len(mth["Geese"]),3)
        self.assertEqual(mth["Utah"],"NOT_TESTED_WITHIN_MONTH_DO_NOT_INFER_ZERO")
        self.assertEqual([v["nominal_sign"] for v in mth["Geese"]],
                         ["INCLUDES_ZERO","INCLUDES_ZERO","NEGATIVE"])

    def test_one_of_the_four_source_DOIs_not_in_twelve_full_source_macro(self):
        d=self.d
        self.assertEqual(d["four_calendar_study_DOIs_already_in_twelve_study_full_source_macro"],3)
        self.assertTrue(d["Illinois_original_DOI_is_NOT_part_of_twelve_full_source_study_frame"])

    def test_previous_negative_radiation_and_old_ERA5_source_HOLDs_unchanged(self):
        d=self.d
        self.assertTrue(d["old_12_ERA5_source_version_raw_JSON_HOLDS_UNCHANGED"])
        self.assertEqual(d["old_48_source_negative_radiation_direction_external_replication_UNCHANGED"],
                         m.NEGATIVE)

    def test_prior_first_source_study_results_fail_closed_on_mutation(self):
        contract=json.loads(m.CONTRACT.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for x in contract["source_result_files"]:
                f=root/x["path"]
                f.parent.mkdir(parents=True,exist_ok=True)
                f.write_bytes((m.ROOT/x["path"]).read_bytes())
            first=root/contract["source_result_files"][0]["path"]
            d=json.loads(first.read_text(encoding="utf-8"))
            d["study_DOI"]="invented"
            first.write_text(json.dumps(d),encoding="utf-8")
            with patch.object(m,"ROOT",root):
                with self.assertRaisesRegex(ValueError,"Previously frozen source original calendar result Git blob changed"):
                    m.analyze()

    def test_no_individual_tag_GPS_event_location_or_date(self):
        s=json.dumps(self.d)
        for token in ('"individual_id"','"animal_id"','"tag_id"','"latitude"',
                      '"longitude"','"location_lat"','"location_long"',
                      '"source_private_gps"','"timestamp"','"source_event_date"'):
            self.assertNotIn(token,s)


if __name__=="__main__":
    unittest.main()
