"""Synthetic source-common-calendar month fixed-effect ecology known truth tests."""
from __future__ import annotations

from datetime import date,timedelta
import math
import unittest

from scripts.movebank_native_ecology_20261019 import daylength
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import (
    compare_frozen_datewise_prior, interpret_month_adjustment,
    month_conditional_slope, sufficient_months, sign
)


def fake_sources(nmonth=7, daysper=14, f=None):
    out={}
    for m in range(1,nmonth+1):
        for j in range(daysper):
            d=date(2021,m,4)+timedelta(days=j)
            L=daylength(43,d)
            val=f(m,L,j) if f else 2*L+(m%2)*7
            out[(d.year,d)]=val
    return out


class RSDSMonthConditionalPhotoperiodKnownTruth(unittest.TestCase):
    def test_within_month_slope_identifies_synthetic_within_block_component(self):
        dates=fake_sources()
        a=month_conditional_slope(dates,43,daylength,seed=20261029)
        self.assertEqual(a["status"],"MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE")
        self.assertEqual(a["eligible_source_month_blocks"],7)
        self.assertEqual(a["eligible_source_calendar_days"],98)
        self.assertAlmostEqual(a["within_month_daylength_slope"],2.0,places=10)
        self.assertAlmostEqual(a["month_cluster_bootstrap95"][0],2,places=10)
        self.assertAlmostEqual(a["month_cluster_bootstrap95"][1],2,places=10)
        self.assertEqual(a["nominal_interval_sign"],"POSITIVE")
        self.assertEqual(a,month_conditional_slope(dates,43,daylength,seed=20261029))

    def test_strong_between_month_changes_do_not_become_within_month_daylength(self):
        dates=fake_sources(f=lambda m,L,j:100*m)
        a=month_conditional_slope(dates,43,daylength)
        self.assertEqual(a["status"],"MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE")
        self.assertAlmostEqual(a["within_month_daylength_slope"],0,places=12)
        self.assertEqual(a["nominal_interval_sign"],"INCLUDES_ZERO")

    def test_negative_within_month_component_recovers_against_positive_month_shifts(self):
        dates=fake_sources(f=lambda m,L,j:-3*L+40*m)
        a=month_conditional_slope(dates,43,daylength)
        self.assertAlmostEqual(a["within_month_daylength_slope"],-3,places=10)
        self.assertEqual(a["nominal_interval_sign"],"NEGATIVE")
        self.assertTrue(all(math.isfinite(z) for z in a["leave_one_month_block_out_range"]))

    def test_insufficient_source_months_are_HOLD_not_zero(self):
        dates=fake_sources(nmonth=3,daysper=13)
        x=month_conditional_slope(dates,43,daylength)
        self.assertEqual(x["status"],"INSUFFICIENT_COMMON_SOURCE_MONTH_SUPPORT")
        self.assertIsNone(x["within_month_daylength_slope"])
        self.assertIsNone(x["month_cluster_bootstrap95"])
        self.assertEqual(x["nominal_interval_sign"],"UNAVAILABLE")

    def test_comparison_to_frozen_first_result_is_nonnegotiable(self):
        self.assertEqual(compare_frozen_datewise_prior(.023,.023),"ARCHIVED_UNADJUSTED_SLOPE_REPLAYED")
        with self.assertRaisesRegex(ValueError,"modified"):
            compare_frozen_datewise_prior(.023,.024)

    def test_eligible_source_months_selected_on_date_support_not_response_values(self):
        a=fake_sources()
        b={d:-200*x for d,x in a.items()}
        blocks1,one=sufficient_months(a,43,daylength)
        blocks2,two=sufficient_months(b,43,daylength)
        for k in ("source_calendar_dates","source_year_month_blocks",
                  "eligible_source_month_blocks","eligible_source_calendar_days",
                  "within_month_daylength_centered_sum_squares_hours2"):
            self.assertEqual(one[k],two[k])
        self.assertEqual(one["status"],two["status"])

    def test_nominal_source_sign_claims_are_exploratory(self):
        self.assertIn("WITHIN_MONTH",interpret_month_adjustment("NOMINAL_POSITIVE","POSITIVE"))
        self.assertEqual(interpret_month_adjustment("POSITIVE","INCLUDES_ZERO"),
                         "WITHIN_MONTH_INTERVAL_INCLUDES_ZERO")
        self.assertEqual(interpret_month_adjustment("NEGATIVE","UNAVAILABLE"),
                         "WITHIN_MONTH_SUPPORT_UNAVAILABLE")


if __name__=="__main__":
    unittest.main()
