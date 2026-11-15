"""Known-truth year-month FE, source-MD5 HOLD and original study-ID accounting."""
from __future__ import annotations
from datetime import date,timedelta
from unittest.mock import patch
import unittest
from scripts import rsds_five_source_native_tag_month_conditioned_v1 as m
from scripts import aggregate_rsds_five_source_native_month_v1 as a
from scripts.movebank_native_ecology_20261019 import daylength


def tagged_daily(beta=.2,ref=43.0):
    rows={}
    for month in (1,3,9):
        start=date(2020,month,5)
        for i in range(12):
            d=start+timedelta(days=i)
            L=daylength(ref,d)
            # Controlled year-month fixed intercept.
            q=2*month+beta*L
            rows[d]=(q,L)
    return {2020:rows}


def source_row(source,i,status="MONTH_NATIVE_ESTIMABLE"):
    p={
        "version":m.SOURCE_VERSION,"source_index":i,
        "original_study_DOI":source["doi"],"original_study_id":str(source["study_id"]),
        "original_native_target":source["taxon"],
        "all_original_registered_native_taxa":source["registered_taxa"],
        "old_original_source_GPS_access_state":source["source_file_prior_access"],
        "baseline_full_year_first_results_replayed":status!="SOURCE_HOLD",
        "original_file_MD5_bytes_verified":status!="SOURCE_HOLD",
        "status":status,
        "original_source_hold_stage":None,
        "original_source_hold_type":None,
        "month_conditioned_native":None,
        "old_solar_direction_external_replication_still_negative":True,
        "old_original_ERA5_raw_SHA_hold":True,
        "no_raw_animal_GPS_tracks_ID_exact_coordinates_dates_or_tag_slopes_exported":True
    }
    if status!="SOURCE_HOLD":
        same={"status":"ESTIMABLE","tagged_animals_support":6,
              "mean_slope_per_hour":.12,"original_tag_bootstrap95":[.11,.13],
              "nominal_interval_sign":"POSITIVE"}
        p["month_conditioned_native"]={
            "local_original_solar_daylength":same,
            "fixed_source_reference_latitude_daylength":same,
        }
    return p


class StudyMonthControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract=m.frozen()
        cls.sources=cls.contract["source_manifest"]

    def test_all_five_frozen_study_sources_and_target_effects(self):
        self.assertEqual(len(self.sources),5)
        self.assertEqual([s["taxon"] for s in self.sources],
                         ["Anser albifrons","Odocoileus hemionus","Odocoileus hemionus",
                          "Anser albifrons","Odocoileus hemionus"])
        self.assertEqual(len({z["doi"] for z in self.sources}),5)
        self.assertEqual(self.contract["baseline_replay_tolerance"],1e-10)

    def test_within_month_fixed_effect_removes_month_intercepts(self):
        year=tagged_daily(.2)
        s=m.group_month_stats(year[2020],43.0)
        self.assertTrue(s["supported"])
        self.assertEqual(s["eligible_months"],3)
        self.assertAlmostEqual(s["local_slope"],.2,places=9)
        self.assertAlmostEqual(s["reference_slope"],.2,places=9)

    def test_positive_and_negative_tag_equal_bootstrap(self):
        p={f"tag-{i}":tagged_daily(.2+0.001*i) for i in range(5)}
        z=m.same_source_month_effect(p,43.0,1)
        self.assertEqual(z["original_tagged_animals_with_year_month_conditions"],5)
        self.assertEqual(z["local_original_solar_daylength"]["nominal_interval_sign"],"POSITIVE")
        neg={f"tag-{i}":tagged_daily(-.3-.001*i) for i in range(5)}
        z=m.same_source_month_effect(neg,43.0,1)
        self.assertEqual(z["local_original_solar_daylength"]["nominal_interval_sign"],"NEGATIVE")

    def test_single_month_or_insufficient_independent_tags_hold(self):
        a0=tagged_daily(.2)[2020]
        one={d:v for d,v in a0.items() if d.month==1}
        x=m.group_month_stats(one,43.0)
        self.assertFalse(x["supported"])
        self.assertIsNone(x["local_slope"])
        tags={f"id{i}":tagged_daily() for i in range(4)}
        x=m.same_source_month_effect(tags,43.0,1)
        self.assertEqual(x["local_original_solar_daylength"]["status"],
                         "INSUFFICIENT_NATIVE_TAG_SUPPORT")
        self.assertIsNone(x["local_original_solar_daylength"]["mean_slope_per_hour"])

    def test_source_original_MD5_failure_does_not_create_month_effect(self):
        with patch.object(m,"original_files",side_effect=OSError("NO_NETWORK_OR_MD5")) as func:
            x=m.process(0)
            func.assert_called_once()
        self.assertEqual(x["status"],"SOURCE_HOLD")
        self.assertEqual(x["original_source_hold_type"],"OSError")
        self.assertFalse(x["original_file_MD5_bytes_verified"])
        self.assertIsNone(x["month_conditioned_native"])
        self.assertNotIn("NO_NETWORK_OR_MD5",str(x))

    def test_aggregation_retains_all_five_and_preserves_old_negative(self):
        good=[source_row(s,i) for i,s in enumerate(self.sources)]
        out=a.aggregate(good)
        self.assertEqual(out["five_frozen_original_field_source_DOI"],5)
        self.assertEqual(out["two_previously_targeted_species"],2)
        self.assertEqual(len(out["all_five_frozen_original_source_native_analyses_including_HOLD"]),5)
        self.assertTrue(out["old_48_solar_direction_external_nonreplication_UNCHANGED"])
        good[0]=source_row(self.sources[0],0,"SOURCE_HOLD")
        out=a.aggregate(good)
        self.assertEqual(out["source_status_count"]["SOURCE_HOLD"],1)
        with self.assertRaisesRegex(ValueError,"five"):
            a.aggregate(good[:4])


if __name__=="__main__":
    unittest.main()
