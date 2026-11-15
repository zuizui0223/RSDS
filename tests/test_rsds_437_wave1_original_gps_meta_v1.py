"""Freeze new source-first effects and meta denominators under known-truth."""
from __future__ import annotations

import copy
import unittest

from scripts.rsds_437_wave1_original_file_metadata_v1 import freeze
from scripts.aggregate_rsds_437_wave1_original_gps_v1 import meta_snapshot

def fake(w):
    a=[]
    for i,row in enumerate(w):
        is_pre=i in (0,1,13)
        a.append({
          "version":"RSDS_437_WAVE1_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1",
          "index":i,"original_ordinal":row["ordinal"],
          "source_package_uuid":row["package_uuid"],
          "original_repository_study_id":row["original_repository_study_id"],
          "source_DOI_label":row["source_DOI_label"],
          "source_selection_stratum":row["selection_stratum"],
          "registered_native_taxa":row["scientific_names"],
          "original_file_preflight_status":(
              "SOURCE_METADATA_HOLD" if is_pre else
              "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"),
          "original_provider_MD5_validated":not is_pre,
          "source_status":("PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD" if is_pre
                            else "ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO"),
          "source_hold_stage":("SOURCE_ONLY_ORIGINAL_FILE_METADATA" if is_pre else None),
          "source_hold_reason":("FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES"
                                if is_pre else None),
          "native_effects":[{"native_taxon":name,"q":None,"daylength":None}
                            for name in row["scientific_names"]],
          "native_common_calendar_covariance":{"status":"NOT_EVALUATED"},
          "observed_original_source_event_rows_if_opened":0,
          "old_48_signed_solar_shortwave_external_NONREPLICATED_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
          "old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED":True,
          "no_private_GPS_or_individual_data_in_output":True
        })
    return a

def effect(value):
    return {"n_individuals":6,"mean":value,
            "animal_bootstrap_95_percentile":[value-.05,value+.05]}

class Stage3SourceMeta(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,cls.sources=freeze()
        cls.data=fake(cls.sources)

    def test_all_sixteen_frozen_slots_with_no_response_no_false_positive(self):
        d=meta_snapshot(self.data)
        self.assertEqual(d["original_catalog_package_count"],437)
        self.assertEqual(d["original_GPS_GNSS_candidates"],259)
        self.assertEqual(d["metadata_eligible_finite_frame"],120)
        self.assertEqual(d["original_selection_wave_slots"],16)
        self.assertEqual(d["original_source_metadata_file_HOLD_slots"],3)
        self.assertEqual(sum(d["source_status_counts"].values()),16)
        for typ in ("q","daylength"):
            s=d["stratified_source_complete_effect_sensitivity"][typ]
            self.assertEqual(s["first_wave_source_complete_response"],0)
            self.assertIsNone(s["response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"])
            self.assertTrue(s["nonresponding_and_unsampled_sources_do_NOT_become_zero"])

    def test_poststratification_uses_109_11_not_8_8(self):
        data=copy.deepcopy(self.data)
        for i,v in [(2,.4),(8,-.3)]:
            row=data[i]
            row["source_status"]="ORIGINAL_SOURCE_GPS_Q_ESTIMABLE"
            for taxon in row["native_effects"]:
                taxon["q"]=effect(v)
                taxon["daylength"]=effect(.01)
        result=meta_snapshot(data)
        q=result["stratified_source_complete_effect_sensitivity"]["q"]
        self.assertEqual(q["first_wave_source_complete_response"],2)
        self.assertAlmostEqual(
          q["response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"],
          109/120*.4 + 11/120*(-.3),places=12
        )
        self.assertEqual(q["source_complete_response_counts"],
                         {"single_native":1,"multi_native":1})
        self.assertTrue(q["hypothetical_fill_missing_selected_sources_with_observed_range"][0]
                        <q["hypothetical_fill_missing_selected_sources_with_observed_range"][1])
        self.assertTrue(result["not_a_world_population_meta_analysis"])

    def test_partial_multi_native_is_not_complete_source(self):
        data=copy.deepcopy(self.data)
        a=data[8]
        a["source_status"]="ORIGINAL_SOURCE_GPS_Q_ESTIMABLE"
        a["native_effects"][0]["q"]=effect(.99)
        d=meta_snapshot(data)
        study=d["all_first_wave_source_study_summaries"][8]
        self.assertEqual(study["native_q_supported"],1)
        self.assertIsNone(study["source_q_equal_native_if_complete"])
        self.assertFalse(study["registered_source_complete_q"])
        self.assertEqual(d["stratified_source_complete_effect_sensitivity"]["q"][
          "first_wave_source_complete_response"],0)

    def test_no_best_species_or_source_substitution_and_old_science(self):
        data=copy.deepcopy(self.data)
        del data[8]["native_effects"][-1]
        with self.assertRaisesRegex(ValueError,"Source selection"):
            meta_snapshot(data)
        data=copy.deepcopy(self.data)
        data[2]["source_package_uuid"]="a-postoutcome-replacement"
        with self.assertRaisesRegex(ValueError,"selection"):
            meta_snapshot(data)
        data=copy.deepcopy(self.data)
        data[4]["old_48_signed_solar_shortwave_external_NONREPLICATED_UNCHANGED"]="REPLICATED"
        with self.assertRaisesRegex(ValueError,"reclassified"):
            meta_snapshot(data)
        data=copy.deepcopy(self.data)
        data[0]["source_status"]="ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO"
        with self.assertRaisesRegex(ValueError,"3 original"):
            meta_snapshot(data)

    def test_hypothetical_effects_do_not_make_earlier_three_HOLDS_zero(self):
        data=copy.deepcopy(self.data)
        for i in (0,1,13):
            self.assertIsNone(data[i]["native_effects"][0]["q"])
            self.assertEqual(data[i]["observed_original_source_event_rows_if_opened"],0)
        d=meta_snapshot(data)
        self.assertEqual(d["source_status_counts"]["PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD"],3)

if __name__=="__main__":
    unittest.main()
