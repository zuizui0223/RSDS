"""Study DOI native species/pair HOLD-aware source atlas known-truth validation."""
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_four_new_multinative_temporal_v1 import aggregate
from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import fixed_sources,NEG,SOURCE_KEY


def held_studies():
    spec,metadata=fixed_sources()
    out={}
    for i,x in enumerate(spec["sources"]):
        out[i]={
            "version":SOURCE_KEY,"source_study_index":i,
            "original_study_DOI":x["doi"],
            "registered_native_taxa":x["registered_taxa"],
            "status":"ORIGINAL_SOURCE_OR_NATIVE_SHARED_CALENDAR_HOLD",
            "original_MD5_verified":False,
            "source_file_metadata_checked_again":False,
            "original_native_taxon_results":[{"taxon":t,"status":"SOURCE_HOLD"}
                                             for t in x["registered_taxa"]],
            "shared_source_calendar_support":None,
            "all_native_study_pair_shared_calendar_results":[],
            "previous_12_ERA5_source_version_HOLDS_UNCHANGED":True,
            "previous_48_source_negative_directional_radiation_replication_UNCHANGED":NEG,
        }
    return out


def e(n,mean):
    return {
        "qualified_original_tagged_animals":n,
        "mean":mean if n>=5 else None,
        "source_original_tag_bootstrap95":[mean-.02,mean+.02] if n>=5 else None,
        "nominal_interval_sign":"POSITIVE" if n>=5 and mean>.02 else
                               "NEGATIVE" if n>=5 and mean<-.02 else
                               "INCLUDES_ZERO" if n>=5 else "UNAVAILABLE",
    }


class FourNewNativeSourceStudyEcologyValidation(unittest.TestCase):
    def test_exact_source_denominator_with_four_true_HOLDS(self):
        d=aggregate(held_studies())
        self.assertEqual(d["fixed_original_published_studies"],4)
        self.assertEqual(d["registered_original_study_native_taxon_rows"],10)
        self.assertEqual(d["study_status_counts"],{"ORIGINAL_SOURCE_OR_NATIVE_SHARED_CALENDAR_HOLD":4})
        self.assertEqual(len(d["all_frozen_source_study_original_results"]),4)
        self.assertTrue(d["prior_12_ERA5_version_source_HOLDS_UNCHANGED"])
        self.assertEqual(d["supported_native_and_source_pair_nominal_sign_counts"],{})

    def test_complete_source_pair_and_native_results_need_MD5_and_full_comparison(self):
        s=held_studies()
        x=s[1]
        x["status"]="ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED"
        x["original_MD5_verified"]=True
        x["source_file_metadata_checked_again"]=True
        x["shared_source_calendar_support"]={
            "source_calendar_dates_shared_by_all_registered_taxa":95,
            "source_yearmonth_blocks":6,"registered_source_taxa":2,
            "source_registered_native_q_tag_counts":{t:6 for t in x["registered_native_taxa"]}
        }
        for t in x["original_native_taxon_results"]:
            t.update(status="ORIGINAL_NATIVE_Q_EVALUATED",
                     original_q_full_source_tag_effect={"n_individuals":6,"mean":.03},
                     original_q_daylength_full_source_tag_effect={"n_individuals":6,"mean":.02},
                     original_source_animals_with_diel_q=6,
                     original_source_animals_with_q_daylength=6,
                     same_study_all_native_common_calendar_q=e(6,.04),
                     same_study_all_native_common_calendar_daylength_slope=e(6,.025))
        t1,t2=x["registered_native_taxa"]
        x["all_native_study_pair_shared_calendar_results"]=[{
            "native_first":t1,"native_second":t2,
            "same_source_calendar_dates_across_ALL_taxa":95,
            "source_original_animal_q_pair":{
                "native_tag_support":[6,6],
                "first_minus_second_q":.1,
                "original_independent_tag_bootstrap95":[.03,.17],
                "nominal_interval_sign":"POSITIVE"
            },
            "source_common_day_q_difference_daylength":{
                "source_shared_calendar_days":95,
                "within_month_reference_daylength":{
                    "nominal_interval_sign":"INCLUDES_ZERO",
                }
            },
            "not_proven_spatial_sympatry":True,
        }]
        d=aggregate(s)
        self.assertEqual(d["study_status_counts"]["ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED"],1)
        self.assertEqual(d["supported_native_and_source_pair_nominal_sign_counts"]["native_q_POSITIVE"],2)
        self.assertEqual(d["supported_native_and_source_pair_nominal_sign_counts"]["pair_q_POSITIVE"],1)
        x["original_MD5_verified"]=False
        with self.assertRaises(ValueError):
            aggregate(s)

    def test_source_taxon_or_pub_DOI_replacement_fails(self):
        s=held_studies()
        s[0]["original_native_taxon_results"].pop()
        with self.assertRaises(ValueError):
            aggregate(s)
        s=held_studies()
        s[3]["original_study_DOI"]="10.5441/001/1.649"
        with self.assertRaises(ValueError):
            aggregate(s)

    def test_preserve_old_shortwave_negative_external_replication(self):
        s=held_studies()
        s[0]["previous_48_source_negative_directional_radiation_replication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate(s)
        self.assertEqual(NEG,"NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")


if __name__=="__main__":
    unittest.main()
