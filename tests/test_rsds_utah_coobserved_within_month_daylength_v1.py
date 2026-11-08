"""Known-truth 199-day Utah coobserved pair/month sensitivity, source SHA locked."""
from __future__ import annotations

import json
import math
import unittest
from datetime import date
from unittest.mock import patch

from scripts.movebank_native_ecology_20261019 import daylength
from scripts.movebank_utah_diel_photoperiod_20261019 import slope
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import month_conditional_slope
from scripts.run_rsds_utah_coobserved_within_month_daylength_v1 import (
    TAXA, SCHEDULES, REFERENCE_LATITUDE, audit_dates_and_months,
    compare_months, frozen
)


def same_utah_199_day_synthetic(known_slope=.025):
    """10 predeclared year×month blocks, 2 calendar years, exactly 199 dates."""
    dates={}
    for year in (2020,2021):
        for m in (4,5,6,7,8):
            ds=[date(year,m,i) for i in range(1,21)]
            Ls=[daylength(REFERENCE_LATITUDE,d) for d in ds]
            mid=sum(Ls)/len(Ls)
            for i,d in enumerate(ds):
                if year==2021 and m==8 and i==19:
                    continue
                sunlight=Ls[i]
                # Within each month: exactly fixed known slope; between
                # months: strong arbitrary environment baseline by month.
                offset=(m-6)*.11 + (.04 if year==2021 else 0)
                q=offset+known_slope*(sunlight-mid)
                dates[(year,d)]={s:q+(.001 if s=="even" else -.001 if s=="odd" else 0)
                                for s in SCHEDULES}
    assert len(dates)==199
    return dates


class UtahOriginalTwoSpeciesWithinMonthKnownTruth(unittest.TestCase):
    def test_original_study_and_month_sampling_contract_sha_lock(self):
        spec, prior_spec, previous_absolute, previous_delta, original=frozen()
        self.assertEqual(spec["version"],"RSDS_UTAH_TWO_NATIVE_SHARED_CALENDAR_WITHIN_MONTH_DAYLENGTH_V1")
        self.assertEqual(original["study_Doi"],"10.5441/001/1.711")
        self.assertEqual(original["native_taxa"],list(TAXA))
        self.assertEqual(original["shared_support"][
            "original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"],199)
        self.assertEqual(original["shared_support"][
            "original_tags_with_10_days_on_both_taxa_shared_calendar"],
            {TAXA[0]:10,TAXA[1]:42})
        self.assertEqual(original["calendar_date_paired_source_difference"]["shared_year_month_blocks"],10)
        self.assertEqual(spec["month_conditioning"]["min_blocks"],5)

    def test_synthetic_month_slope_recovers_known_true_within_month_not_between_month(self):
        dates=same_utah_199_day_synthetic(known_slope=.025)
        vals={key:x["original"] for key,x in dates.items()}
        g=month_conditional_slope(vals,REFERENCE_LATITUDE,daylength,seed=20261108)
        self.assertEqual(g["source_calendar_dates"],199)
        self.assertEqual(g["source_year_month_blocks"],10)
        self.assertGreaterEqual(g["eligible_source_month_blocks"],5)
        self.assertEqual(g["status"],"MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE")
        self.assertAlmostEqual(g["within_month_daylength_slope"],.025,places=10)
        self.assertTrue(g["month_cluster_bootstrap95"][0]>0)
        self.assertGreater(g["within_month_daylength_centered_sum_squares_hours2"],.25)
        naive=slope([(daylength(REFERENCE_LATITUDE,d),v) for (_,d),v in vals.items()])
        self.assertGreater(abs(naive-.025),.01)

    def test_replay_three_old_unadjusted_slopes_before_new_month_effects(self):
        data=same_utah_199_day_synthetic()
        ref={
            "calendar_date_paired_source_difference":{
                "coobserved_original_source_dates":199,
                "shared_year_month_blocks":10,
                "slopes":{
                    s:{"mule_q_minus_bighorn_q_vs_shared_reference_daylength":
                        slope([(daylength(REFERENCE_LATITUDE,d),v[s]) for (_,d),v in data.items()]),
                        "status":"INCLUDES_ZERO"}
                    for s in SCHEDULES
                }
            }
        }
        report=compare_months(data,ref,{
            "source_shared_original_calendar_days":199,"original_source_year_month_blocks":10
        })
        self.assertEqual(len(report),3)
        self.assertEqual([x["source_GPS_UTC_clock_schedule"] for x in report],list(SCHEDULES))
        self.assertTrue(all(x["all_months_prior_source_slope_replay"]==
                            "ARCHIVED_UNADJUSTED_SLOPE_REPLAYED" for x in report))
        self.assertTrue(all(abs(x["month_centered_source_daylength_response"][
            "within_month_daylength_slope"]-.025)<1e-10 for x in report))
        ref["calendar_date_paired_source_difference"]["slopes"]["odd"][
            "mule_q_minus_bighorn_q_vs_shared_reference_daylength"]+=.01
        with self.assertRaisesRegex(ValueError,"modified"):
            compare_months(data,ref,{"source_shared_original_calendar_days":199,
                                      "original_source_year_month_blocks":10})

    def test_month_support_cannot_rescue_a_different_date_or_taxon_sample(self):
        _, _, _, _, prior=frozen()
        aligned={taxon:{"private_synthetic_"+str(k):{}
                         for k in range(n)} for taxon,n in ((TAXA[0],10),(TAXA[1],42))}
        original_days=set(same_utah_199_day_synthetic())
        r=audit_dates_and_months(aligned,original_days,prior)
        self.assertEqual(r["source_shared_original_calendar_days"],199)
        self.assertEqual(r["original_source_year_month_blocks"],10)
        self.assertFalse(r["original_private_dates_or_tag_ids_exported"])
        altered=set(original_days)
        altered.pop()
        with self.assertRaisesRegex(ValueError,"changed"):
            audit_dates_and_months(aligned,altered,prior)
        del aligned[TAXA[0]]["private_synthetic_0"]
        with self.assertRaisesRegex(ValueError,"changed"):
            audit_dates_and_months(aligned,original_days,prior)

    def test_historic_negative_radiation_and_ERA5_hold_remain_immutable(self):
        _, _, _, _, prior=frozen()
        self.assertTrue(prior["historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED"])
        self.assertEqual(prior["historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")


if __name__=="__main__":
    unittest.main()
