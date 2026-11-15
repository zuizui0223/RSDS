"""Known-truth tests for postoutcome same-source matched month decomposition."""
from __future__ import annotations
from datetime import date,timedelta
from unittest.mock import patch
import unittest

from scripts.movebank_native_ecology_20261019 import daylength
from scripts import rsds_five_source_same_tag_date_month_decomposition_v1 as m


def synthetic_year(beta=.2,lat=43):
    days={}
    for month,intercept in [(1,0),(3,15),(9,-10)]:
        first=date(2020,month,4)
        for i in range(13):
            d=first+timedelta(days=i)
            L=daylength(lat,d)
            days[d]=(intercept+beta*L,L)
    return {2020:days}


class SameOriginalTagDateMonth(unittest.TestCase):
    def test_same_date_unadjusted_and_FE_differ_due_to_month_intercepts(self):
        z=m.original_tag_year_selected_same_days(synthetic_year(.2)[2020],43)
        self.assertIsNotNone(z)
        self.assertEqual(z["same_original_dates"],39)
        self.assertAlmostEqual(z["local_month_FE_same_dates"],.2,places=8)
        self.assertAlmostEqual(z["reference_month_FE_same_dates"],.2,places=8)
        self.assertNotAlmostEqual(z["local_raw_same_dates"],.2,places=2)

    def test_original_tag_resampling_is_paired_and_same_support(self):
        animals={f"tag{i}":synthetic_year(.1+.003*i) for i in range(6)}
        r=m.source_matched_tag_effect(animals,43,20261115)
        self.assertTrue(r["all_contrasts_use_identical_original_tag_year_month_date_support"])
        self.assertEqual(r["source_original_eligible_tagged_animals"],6)
        fx=r["source_q_daylength_effects"]
        for key in ("local_raw_same_dates","local_month_FE_same_dates","local_month_minus_raw"):
            self.assertEqual(fx[key]["tag_n"],6)
            self.assertEqual(fx[key]["status"],"ESTIMABLE")
        self.assertAlmostEqual(
            fx["local_month_FE_same_dates"]["mean"]-fx["local_raw_same_dates"]["mean"],
            fx["local_month_minus_raw"]["mean"],places=9)
        self.assertEqual(fx["local_month_FE_same_dates"]["nominal_sign"],"POSITIVE")

    def test_month_support_unavailable_is_not_zero(self):
        one=synthetic_year()[2020]
        only={d:z for d,z in one.items() if d.month==1}
        self.assertIsNone(m.original_tag_year_selected_same_days(only,43))
        r=m.source_matched_tag_effect({f"tag{i}":synthetic_year() for i in range(4)},43,20261115)
        z=r["source_q_daylength_effects"]["local_month_minus_raw"]
        self.assertEqual(z["status"],"INSUFFICIENT_NATIVE_TAG_SUPPORT")
        self.assertIsNone(z["mean"])
        self.assertEqual(z["nominal_sign"],"UNAVAILABLE")

    def test_original_5_studies_and_prior_first_month_blob_pinned(self):
        c,old,prior=m.fixed()
        self.assertEqual(c["original_same_source_native_study_slots"],5)
        self.assertEqual(len(old["source_manifest"]),5)
        self.assertEqual(len(prior),5)
        self.assertTrue(all(x["source_gps_status"]=="MONTH_NATIVE_ESTIMABLE" for x in prior))
        self.assertEqual(c["prior_five_source_month_first_result_git_blob"],
                         "27643d9c69ca6683edde5ee75945f6a4fb0eb413")

    def test_original_provider_failure_stops_before_any_gps(self):
        with patch.object(m,"original_files",side_effect=OSError("URL_WITH_PRIVATE_INFO")) as fetch:
            z=m.source_process(0)
            fetch.assert_called_once()
        self.assertEqual(z["source_status"],"SOURCE_HOLD")
        self.assertFalse(z["old_q_and_daylength_original_baseline_replayed"])
        self.assertFalse(z["original_file_UUID_MD5_bytes_verified"])
        self.assertIsNone(z["matched_month_effect"])
        self.assertEqual(z["source_hold_stage"],"ORIGINAL_PROVIDER_FROZEN_SOURCE_UUID_BYTES_MD5")
        self.assertNotIn("URL_WITH_PRIVATE_INFO",str(z))


if __name__=="__main__":
    unittest.main()
