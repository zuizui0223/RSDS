"""Complete original six-study macro source-denominator + null tests."""
import unittest
from unittest.mock import patch

from scripts.aggregate_rsds_new_six_original_gps_ecology_v0 import (
    summarize, supported
)
from scripts.run_rsds_new_six_original_gps_ecology_v0 import NEGATIVE


def frame():
    ss=[]
    ee=[]
    for i in range(6):
        ss.append({"id":f"source{i}","doi":f"10.5441/001/1.fake{i}",
                   "taxon":f"Native species{i}","study_id":f"study{i}",
                   "guild":"guild","region":"region"})
        ee.append({"package_uuid":f"source{i}",
                   "source_metadata_status":"ORIGINAL_FILE_METADATA_ADMISSIBLE"})
    return ss,ee


def outcomes():
    sources,meta=frame()
    results={}
    for i,s in enumerate(sources):
        results[str(i)]={
            "version":"RSDS_NEW_SIX_ORIGINAL_GPS_FIRST_STUDY_RESULT_V0",
            "native_taxon":s["taxon"],
            "source_doI":s["doi"],
            "original_study_id":s["study_id"],
            "source_metadata_status":"ORIGINAL_FILE_METADATA_ADMISSIBLE",
            "source_status":"ORIGINAL_SOURCE_GPS_ECOLOGY_EVALUATED",
            "old_negative_external_shortwave_replication_UNCHANGED":NEGATIVE,
            "all_outcome_axes":{
                "day_night_q":None,"daylength_slope":None,"sunward_raw":None,
                "sunward_joint":None,"same_tag_two_contexts":None
            }
        }
    return results


class SixOriginalPublicSourceSynthesisKnownTruth(unittest.TestCase):
    def test_no_original_results_is_six_visible_UNAVAILABLE_sources(self):
        with patch("scripts.aggregate_rsds_new_six_original_gps_ecology_v0.load_frozen",
                   return_value=frame()):
            x=summarize(outcomes())
        self.assertEqual(len(x["six_full_source_study_taxon_results"]),6)
        self.assertEqual(x["source_sample_sign_support_counts"]["q"],{"UNAVAILABLE":6})
        self.assertEqual(x["source_sample_sign_support_counts"]["paired"],{"UNAVAILABLE":6})
        self.assertEqual(x["previous_external_48_study_solar_radiation_nonreplication_UNCHANGED"],NEGATIVE)

    def test_supported_source_mean_bootstrap_has_positive_and_null(self):
        self.assertEqual(supported({"n_individuals":6,"mean":.4,
                                    "animal_bootstrap_95_percentile":[.1,.8]},
                                   "mean","animal_bootstrap_95_percentile")["status"],
                         "POSITIVE_NOMINAL_BOOTSTRAP")
        self.assertEqual(supported({"n_individuals":5,"mean":-.05,
                                    "bootstrap95":[-.2,.3]},
                                   "mean","bootstrap95")["status"],
                         "INCLUDES_ZERO")
        with self.assertRaises(ValueError):
            supported({"n_individuals":4,"mean":.2,"bootstrap95":[.1,.3]},
                      "mean","bootstrap95")

    def test_no_replacement_or_duplicated_source_after_any_hold(self):
        r=outcomes()
        r["1"]["source_status"]="ORIGINAL_METADATA_HOLD"
        r["1"]["source_metadata_status"]="SOURCE_METADATA_HOLD"
        sources,meta=frame()
        meta[1]["source_metadata_status"]="SOURCE_METADATA_HOLD"
        with patch("scripts.aggregate_rsds_new_six_original_gps_ecology_v0.load_frozen",
                   return_value=(sources,meta)):
            z=summarize(r)
        self.assertEqual(z["number_of_studies_by_source_status"]["ORIGINAL_METADATA_HOLD"],1)
        self.assertEqual(z["source_sample_sign_support_counts"]["q"]["UNAVAILABLE"],6)
        del r["1"]
        with patch("scripts.aggregate_rsds_new_six_original_gps_ecology_v0.load_frozen",
                   return_value=(sources,meta)):
            with self.assertRaises(ValueError):
                summarize(r)

    def test_source_positive_time_does_not_invent_paired_same_tag_or_heading(self):
        r=outcomes()
        r["0"]["all_outcome_axes"]["day_night_q"]={
            "n_individuals":7,"mean":.05,
            "animal_bootstrap_95_percentile":[.01,.08]
        }
        with patch("scripts.aggregate_rsds_new_six_original_gps_ecology_v0.load_frozen",
                   return_value=frame()):
            z=summarize(r)
        self.assertEqual(z["source_sample_sign_support_counts"]["q"],{
            "UNAVAILABLE":5,"POSITIVE_NOMINAL_BOOTSTRAP":1
        })
        self.assertEqual(z["source_sample_sign_support_counts"]["sun_raw"],
                         {"UNAVAILABLE":6})
        self.assertEqual(z["source_sample_sign_support_counts"]["paired"],
                         {"UNAVAILABLE":6})


if __name__=="__main__":
    unittest.main()
