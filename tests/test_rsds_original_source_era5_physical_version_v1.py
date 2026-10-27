"""Offline source GPS/ERA5 canonical provenance and unidentifiability tests."""
from __future__ import annotations
import unittest
from copy import deepcopy

from scripts.audit_rsds_era5_original_source_physical_v1 import (
    context_guard, freeze
)
from scripts.aggregate_rsds_era5_original_source_physical_v1 import aggregate


def fake_region(name):
    contract, original=freeze()
    s=original[name]
    taxa=[]
    for old in s["taxon_results"]:
        rr=[]
        for context in old["selected_contexts"]:
            rr.append({
                "status":"CURRENT_VERSION_TWO_REQUESTS_PHYSICAL_AUDITED",
                "original_source_q_coverage_rank":context["selected_gridyear_original_coverage_rank"],
                "original_q_qualified_source_animals":context["original_q_qualified_animals"],
                "old_raw_ERA5_source_SHA_still_unverified":True,
                "historical_old_weather_physical_array_not_archived":True,
                "raw_response_SHA_identical":False,
                "physical_hourly_numeric_SHA_identical":True,
                "returned_weather_grid_metadata_identical":True,
                "raw_SHA_disagreement_with_identical_physical_values_and_grid":True,
                "first_new_raw_SHA_equals_old_raw_SHA":False,
            })
        taxa.append({"taxon":old["taxon"],
                     "original_temporal_q_reproduced":True,
                     "original_selected_contexts":2,
                     "context_results":rr})
    return {
        "version":"RSDS_ERA5_ORIGINAL_GPS_SOURCE_CURRENT_PHYSICAL_CONTENT_AUDIT_V1",
        "region":name,"source_original_study_DOI":s["original_source_doi"],
        "status":"ORIGINAL_GPS_SOURCE_PHYSICAL_VERSION_AUDIT_COMPLETED",
        "source_original_GPS_MD5_replayed_before_weather":True,
        "all_native_taxon_q_and_daylength_replayed_before_weather":True,
        "native_source_taxa":len(taxa),"taxa":taxa,
        "previous_source_radiation_external_nonreplication_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
        "old_biological_weather_responses_reclassified":False,
    }


class ERA5OriginalStudyPhysicalAudit(unittest.TestCase):
    def test_historical_12_original_raw_SHA_HOLDS_remain_unchanged(self):
        spec,old=freeze()
        self.assertEqual(set(old),{"Utah","Wyoming"})
        self.assertEqual(len(old["Utah"]["taxon_results"]),4)
        self.assertEqual(len(old["Wyoming"]["taxon_results"]),2)
        self.assertEqual(spec["prior_old_raw_12"]["old_raw_JSON_SHA_mismatch_contexts"],12)

    def test_context_replay_checks_preweather_tag_coverage_and_steps(self):
        _,first=freeze()
        old=first["Utah"]["taxon_results"][0]["selected_contexts"][0]
        c={"rank":old["selected_gridyear_original_coverage_rank"],
           "source_q_days":old["original_qualified_q_days"],
           "source_q_animals":old["original_q_qualified_animals"],
           "preweather_solar_steps":old["original_moving_daylight_steps"]}
        context_guard(c,old)
        for k in ("source_q_days","source_q_animals","preweather_solar_steps"):
            changed=deepcopy(c)
            changed[k]+=1
            with self.assertRaises(ValueError):
                context_guard(changed,old)

    def test_six_original_taxa_twelve_contexts_two_studies(self):
        out=aggregate(fake_region("Utah"),fake_region("Wyoming"))
        self.assertEqual(out["original_studies"],2)
        self.assertEqual(out["original_study_taxon_groups"],6)
        self.assertEqual(out["original_selected_weather_contexts_reconstructed"],12)
        self.assertEqual(out["within_current_session_raw_JSON_sha_mismatch"],12)
        self.assertEqual(out["within_current_session_canonical_UTC_hourly_values_equal"],12)
        self.assertTrue(out["historical_original_12_source_version_HOLDs_UNCHANGED"])
        self.assertEqual(out["baseline_physical_content_OLD_vs_NEW"],
                         "NOT_IDENTIFIABLE_OLD_HOURLY_NUMERIC_VALUES_NOT_ARCHIVED")

    def test_missing_old_source_taxon_fails_instead_of_source_replacement(self):
        a,b=fake_region("Utah"),fake_region("Wyoming")
        a["taxa"].pop()
        with self.assertRaises(ValueError):
            aggregate(a,b)

    def test_a_new_provider_weather_hold_is_not_original_ERA5_zero(self):
        a,b=fake_region("Utah"),fake_region("Wyoming")
        z=a["taxa"][0]["context_results"][0]
        z.clear()
        z.update({"status":"NEW_PROVIDER_PHYSICAL_SOURCE_HOLD",
                  "old_raw_ERA5_source_SHA_still_unverified":True})
        out=aggregate(a,b)
        self.assertEqual(out["current_provider_contexts_two_response_comparisons"],11)
        self.assertEqual(out["context_provider_result_status_counts"]["NEW_PROVIDER_PHYSICAL_SOURCE_HOLD"],1)
        self.assertEqual(out["baseline_raw_JSON_SHA_historical_replay_HOLD_contexts"],12)

    def test_falsified_old_external_failure_fails_closed(self):
        a,b=fake_region("Utah"),fake_region("Wyoming")
        b["previous_source_radiation_external_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate(a,b)


if __name__=="__main__":
    unittest.main()
