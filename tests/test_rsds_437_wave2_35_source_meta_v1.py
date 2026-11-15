"""Known truth combined meta estimator with all 16+19 source HOLDS present."""
from __future__ import annotations
import copy
import unittest

from scripts import aggregate_rsds_437_wave2_35_source_meta_v1 as m

def synthetic():
    _,roster=m.freeze_original()
    rows=[]
    for i,s in enumerate(roster):
        names=s["registered_native_taxa"]
        pre=s["source_status"]
        status=("HUMAN_LOCATION_PRIVACY_HOLD" if pre=="HUMAN_LOCATION_PRIVACY_HOLD"
                else "SOURCE_FILE_METADATA_HOLD" if pre=="SOURCE_METADATA_HOLD"
                else "SOURCE_NO_NATIVE_Q_ESTIMABLE")
        rows.append({
          "version":"RSDS_WAVE2_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1",
          "wave2_index":i,
          "source_package_uuid":s["package_uuid"],
          "source_original_repository_study_id":s["original_repository_study_id"],
          "source_original_DOI_labels":s["original_DOI_labels"],
          "source_metadata_stratum":s["selection_stratum"],
          "original_registered_native_taxa":names,
          "original_file_source_status":pre,
          "original_file_MD5_bytes_verified":pre=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE",
          "source_science_status":status,
          "source_hold_stage":"HUMAN_SOURCE_PRIVACY_NO_ACCESS" if status=="HUMAN_LOCATION_PRIVACY_HOLD"
                  else "SOURCE_ONLY_ORIGINAL_FILE_METADATA" if status=="SOURCE_FILE_METADATA_HOLD" else None,
          "source_hold_reason":"NO_HUMAN_TRAJECTORY_SOURCE_ACCESS" if status=="HUMAN_LOCATION_PRIVACY_HOLD" else None,
          "native_group_results":[{"taxon":taxon,"q":None,"daylength":None,"status":"NATIVE_Q_UNAVAILABLE"} for taxon in names],
          "native_joint_calendar_covariance":{"status":"NOT_EVALUATED"},
          "original_source_GPS_event_rows_opened":0,
          "human_mobility_GPS_data_opened":False,
          "old_48_signed_solar_external_nonreplication_UNCHANGED":
              "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
          "old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED":True,
          "original_raw_GPS_dates_tag_IDs_coordinates_exported":False
        })
    return rows

def effect(y):
    return {"mean":y,"n_individuals":8,
            "animal_bootstrap_95_percentile":[y-.08,y+.08]}

class Total35SourceMeta(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=synthetic()

    def test_first_wave_unchanged_and_all_source_slots(self):
        result=m.summarize(self.rows)
        self.assertEqual(result["original_registry_source_packages"],437)
        self.assertEqual(result["original_GPS_GNSS_sensor_packages"],259)
        self.assertEqual(result["original_metadata_candidate_frame"],120)
        self.assertEqual(result["original_source_study_units"],35)
        self.assertEqual(result["wave1_source_study_units"],16)
        self.assertEqual(result["wave2_source_study_units"],19)
        self.assertEqual(len(result["combined_source_study_summary_rows"]),35)
        self.assertEqual(result["source_complete_effects_and_nonresponse_sensitivity"]["q"][
            "supported_complete_native_source_count"],4)
        self.assertTrue(result["human_location_source_no_access"])
        self.assertTrue(result["previous_48_source_sunlight_signed_external_not_replicated"])

    def test_two_stratum_measurability_denominators(self):
        result=m.summarize(self.rows)
        sample=result["finite_120_metadata_frame_source_measurability"]["q"]
        self.assertEqual(sample["sampled_fraction_single"],24/109)
        self.assertEqual(sample["sampled_fraction_multi"],1.0)
        self.assertEqual(sample["by_source_stratum"],{"single_native":2,"multi_native":2})
        self.assertEqual(sample["strict_finite_frame_identified_measurability_bounds"],
                         [4/120,89/120])
        self.assertAlmostEqual(sample["pseudo_random_SHA_frame_weighted_measurability"],
                               (109*2/24+2)/120)

    def test_success_q_and_daylight_on_registered_native_with_no_pseudo_replication(self):
        rows=copy.deepcopy(self.rows)
        i=next(i for i,v in enumerate(rows) if v["source_metadata_stratum"]=="single_native"
               and v["original_file_source_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        rows[i]["native_group_results"][0]["q"]=effect(-.3)
        rows[i]["native_group_results"][0]["daylength"]=effect(.012)
        rows[i]["source_science_status"]="SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE"
        out=m.summarize(rows)
        self.assertEqual(out["source_complete_effects_and_nonresponse_sensitivity"]["q"][
            "supported_complete_native_source_count"],5)
        self.assertEqual(out["finite_120_metadata_frame_source_measurability"]["q"][
            "by_source_stratum"],{"single_native":3,"multi_native":2})
        self.assertEqual(out["source_complete_effects_and_nonresponse_sensitivity"]["q"][
            "missing_complete_native_source_count"],30)
        self.assertIsNotNone(out["source_complete_effects_and_nonresponse_sensitivity"]["q"][
            "respondent_only_poststratified_mean_EXPLORATORY"])

    def test_partial_multi_native_not_passed_off_as_source_complete(self):
        rows=copy.deepcopy(self.rows)
        i=next(i for i,r in enumerate(rows) if r["source_metadata_stratum"]=="multi_native"
               and r["original_file_source_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        rows[i]["native_group_results"][0]["q"]=effect(1.3)
        rows[i]["source_science_status"]="SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE"
        result=m.summarize(rows)
        sample=result["new_wave2_source_study_summary_rows"][i]
        self.assertFalse(sample["complete_source_native_q"])
        self.assertEqual(sample["native_q_supported"],1)
        self.assertIsNone(sample["source_q_mean_equal_registered_native"])

    def test_human_privacy_reclassified_or_source_replaced_fails(self):
        rows=copy.deepcopy(self.rows)
        rows[10]["source_science_status"]="SOURCE_NO_NATIVE_Q_ESTIMABLE"
        with self.assertRaisesRegex(ValueError,"Human GPS private source"):
            m.summarize(rows)
        rows=copy.deepcopy(self.rows)
        rows[0]["source_package_uuid"]="post-result substitute"
        with self.assertRaisesRegex(ValueError,"Wave2 source ID"):
            m.summarize(rows)
        rows=copy.deepcopy(self.rows)
        rows[0]["old_48_signed_solar_external_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaisesRegex(ValueError,"previous biology"):
            m.summarize(rows)

    def test_missing_original_source_slot_not_discarded(self):
        with self.assertRaisesRegex(ValueError,"nineteen"):
            m.summarize(self.rows[:-1])

if __name__=="__main__":
    unittest.main()
