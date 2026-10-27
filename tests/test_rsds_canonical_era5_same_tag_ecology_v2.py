"""Canonical NEW provider-version original tagged animal radiation effect tests."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from copy import deepcopy

from scripts.run_rsds_canonical_era5_same_tag_ecology_v2 import (
    frozen, compare_context_eligibility
)
from scripts.aggregate_rsds_canonical_era5_same_tag_ecology_v2 import (
    run, classify
)


def fake_source(region):
    spec, _, _ = frozen()
    key = "utah" if region == "Utah" else "wyoming"
    info = spec["study_coverage"][key]
    return {
        "version":"RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_ORIGINAL_SOURCE_RESULT",
        "region":region,
        "source_study_DOI":info["doi"],
        "original_GPS_provider_MD5_replayed":True,
        "original_q_and_daylength_replayed_before_weather":True,
        "source_native_taxa":info["native_taxa"],
        "taxa":[{"taxon":t,"status":"CANONICAL_WEATHER_HOLD",
                 "source_animal_aggregate":None,
                 "within_same_tag_rank2_minus_rank1":None,
                 "context_results":[]
                } for t in info["native_taxa"]],
        "status":"ORIGINAL_SOURCE_REANALYSIS_COMPLETED",
        "separate_new_physical_weather_version":True,
        "previous_ERA5_raw_SHA_12_context_HOLDS_UNCHANGED":True,
        "prior_48_source_radiation_external_NONREPLICATED_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
        "new_ERA5_shortwave_effect_is_NOT_replication_of_old_weather_dataset":True
    }


class RSDSCanonicalWeatherOriginalTaggedEcologyTests(unittest.TestCase):
    def test_new_contract_pins_first_actual_GPS_and_current_numeric_source_blobs(self):
        spec,first,audit=frozen()
        self.assertEqual(spec["study_coverage"]["selected_q_gridyear_contexts"],12)
        self.assertEqual(audit["current_provider_contexts_two_response_comparisons"],12)
        self.assertEqual(audit["within_current_session_canonical_UTC_hourly_values_equal"],12)
        self.assertEqual(len(first["Utah"]["taxon_results"]),4)
        self.assertEqual(len(first["Wyoming"]["taxon_results"]),2)

    def test_first_context_q_and_original_solar_steps_are_replayed_before_new_ERA5(self):
        _,first,audit=frozen()
        old=first["Utah"]["taxon_results"][0]["selected_contexts"][0]
        ref=next(x for x in audit["full_original_taxon_context_results"]
                 if x["region"]=="Utah" and x["taxon"]=="Antilocapra americana")["context_outcomes"][0]
        sample={"rank":old["selected_gridyear_original_coverage_rank"],
                "source_q_days":old["original_qualified_q_days"],
                "source_q_animals":old["original_q_qualified_animals"],
                "preweather_solar_steps":old["original_moving_daylight_steps"]}
        compare_context_eligibility(sample,old,ref)
        for key in ("source_q_days","source_q_animals","preweather_solar_steps"):
            fake=deepcopy(sample)
            fake[key]+=1
            with self.assertRaises(ValueError):
                compare_context_eligibility(fake,old,ref)

    def test_source_group_cannot_claim_supported_effect_with_under_five_tags(self):
        value={"n_source_animals":4,"source_mean":.2,
               "source_animal_bootstrap95":[.1,.3]}
        with self.assertRaises(ValueError):
            classify(value)
        self.assertEqual(classify({"n_source_animals":4,
                                   "source_mean":None,"source_animal_bootstrap95":None})["status"],
                         "UNSUPPORTED")
        self.assertEqual(classify({"n_source_animals":6,"source_mean":.3,
                                   "source_animal_bootstrap95":[.1,.5]})["status"],
                         "NOMINAL_POSITIVE_ANIMAL_BOOTSTRAP")
        self.assertEqual(classify({"n_source_animals":6,"source_mean":-.2,
                                   "source_animal_bootstrap95":[-.4,.1]})["status"],
                         "INCLUDES_ZERO")

    def test_six_original_taxa_and_null_holds_without_invented_biological_results(self):
        r=run(fake_source("Utah"),fake_source("Wyoming"))
        self.assertEqual(r["original_study_count"],2)
        self.assertEqual(r["original_source_taxon_rows"],6)
        self.assertEqual(r["selected_contexts_reported"],0)
        self.assertEqual(r["current_canonical_weather_contexts_exact_source_SHA_replayed"],0)
        self.assertEqual(r["source_animal_effect_nominal_sign_support"][
            "sun_cos_vs_irradiance"],{"UNSUPPORTED":6})
        self.assertEqual(r["historical_12_ERA5_raw_JSON_SHA_HOLDs_UNCHANGED"],True)

    def test_original_source_replay_must_precede_any_new_biological_effect(self):
        a,b=fake_source("Utah"),fake_source("Wyoming")
        a["original_q_and_daylength_replayed_before_weather"]=False
        with self.assertRaises(ValueError):
            run(a,b)

    def test_pinned_first_source_weather_SHA_change_fails_closed(self):
        import scripts.run_rsds_canonical_era5_same_tag_ecology_v2 as mod
        spec,first,audit=frozen()
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"fake_physical.json"
            p.write_text(json.dumps({**audit,"original_studies":42}),encoding="utf-8")
            with patch.object(mod,"AUDIT",p):
                with self.assertRaisesRegex(ValueError,"SHA changed"):
                    frozen()


if __name__=="__main__":
    unittest.main()
