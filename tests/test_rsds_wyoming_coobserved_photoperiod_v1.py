"""Known-truth two native taxa on exact same Wyoming source calendar dates.

No real GPS fix, tagged source IDs, dates or coordinates are ever output.
"""
from __future__ import annotations
import json
import math
import unittest
from datetime import date,timedelta
from scripts.run_rsds_wyoming_coobserved_photoperiod_v1 import (
    TAXA,SCHEDULES,frozen,daylength,as_three_original_calendar,
    matched_two_taxon_calendar,per_tag_effects,bootstrap_source,
    between_taxa_effect,same_calendar_date_difference,
    monthly_block_calendar_contrast,replay_original_triple
)


def synthetic(n1=7,n2=8,days=120):
    pool={sp:{} for sp in TAXA}
    date0=date(2020,1,1)
    for species,n,sign in ((TAXA[0],n1,1),(TAXA[1],n2,-1)):
        for i in range(n):
            tag=f"PRIVATE_{species}_{i}"
            rows={}
            for k in range(days):
                d=date0+timedelta(days=k)
                L=daylength(43.,d)
                base=(.02 if sign==1 else -.04)+sign*.02*L+i*.00001
                rows[d]={"original":base,"even":base+.003,
                         "odd":base-.004}
            pool[species][tag]={2020:rows}
    return pool


class WyomingWithinStudySharedDates(unittest.TestCase):
    def test_frozen_original_MD5_study_frame_and_prior_result_immutable(self):
        spec,prior,prior_delta=frozen()
        self.assertEqual(spec["original_study"]["doi"],"10.5441/001/1.703")
        self.assertEqual([x["taxon"] for x in prior["native_source_results"]],
                         list(TAXA))
        self.assertEqual({x["taxon"] for x in prior_delta["taxa"]},set(TAXA))
        self.assertTrue(all("new_strict_triple_original_tag_effects" in x
                            for x in prior_delta["taxa"]))
        self.assertTrue(all("original_same_date_q_source_estimates" in x
                            for x in prior["native_source_results"]))
        self.assertEqual(spec["date_filter"]["minimum_date_level"].split("shared dates")[0].strip(),
                         "For pooled calendar-date study comparison, >=20")

    def test_replay_uses_actual_distinct_archived_native_sign_schema(self):
        from unittest.mock import patch
        spec, absolute, delta=frozen()
        def fake_replay(taxon, original_gps, original_delta_row):
            self.assertIn("new_strict_triple_original_tag_effects",original_delta_row)
            historical=next(x for x in absolute["native_source_results"]
                            if x["taxon"]==taxon)
            return {
                "source_q_sign_recurrence":
                    historical["source_q_three_schedule_nominal_sign_stability"],
                "source_daylength_slope_sign_recurrence":
                    historical["source_daylength_three_schedule_nominal_sign_stability"],
                "q_original_tag_means_on_exact_same_dates":
                    historical["original_same_date_q_source_estimates"],
                "q_daylength_per_hour_on_exact_same_dates":
                    historical["original_same_date_daylength_slope_estimates"],
            }
        with patch("scripts.run_rsds_wyoming_coobserved_photoperiod_v1.replay_first_strict_three",
                   side_effect=fake_replay):
            self.assertTrue(replay_original_triple({s:{} for s in TAXA},absolute,delta))

    def test_all_synthetic_dates_shared_and_fixed_study_photoperiod_x(self):
        pool=synthetic()
        kept,days,n=matched_two_taxon_calendar(pool)
        self.assertEqual(len(days),120)
        self.assertEqual(n,1)
        self.assertEqual(len(kept[TAXA[0]]),7)
        self.assertEqual(len(kept[TAXA[1]]),8)
        a,s1=per_tag_effects(kept[TAXA[0]])
        b,s2=per_tag_effects(kept[TAXA[1]])
        self.assertEqual(len(s1),7)
        self.assertEqual(len(s2),8)
        signs=between_taxa_effect(s1,s2,23)
        for sch in SCHEDULES:
            v=signs["schedules"][sch]
            self.assertAlmostEqual(v["mule_minus_bighorn_slope"],.04,places=12)
            self.assertEqual(v["status"],"NOMINAL_POSITIVE")
        source=bootstrap_source(s1,3)
        self.assertTrue(all(source["schedules"][sch]["nominal_sign"]=="POSITIVE"
                            for sch in SCHEDULES))
        e=bootstrap_source(s2,4)
        self.assertTrue(all(e["schedules"][sch]["nominal_sign"]=="NEGATIVE"
                            for sch in SCHEDULES))

    def test_cross_taxon_daily_data_truly_matched_by_date_year_not_tag(self):
        pool=synthetic()
        kept,dates,_=matched_two_taxon_calendar(pool)
        pairs=same_calendar_date_difference(kept,dates)
        self.assertEqual(len(pairs),120)
        r=monthly_block_calendar_contrast(pairs,42)
        self.assertEqual(r["coobserved_original_source_dates"],120)
        self.assertGreaterEqual(r["shared_year_month_blocks"],4)
        # Study month-block minimum is >=5, and these 120 days cover 4
        # months only, so block-bootstrap is correctly unavailable.
        self.assertTrue(all(
            v["status"]=="INSUFFICIENT_CALENDAR_BLOCK_SUPPORT"
            for v in r["slopes"].values()
        ))

    def test_date_overlap_pruning_not_converted_to_animal_nonresponse(self):
        pool=synthetic()
        # Move one species' dates far enough away to prevent overlap.
        t=TAXA[1]
        for tag,years in pool[t].items():
            year=years[2020]
            years[2020]={d+timedelta(days=200):v for d,v in year.items()}
        filtered,shared,n=matched_two_taxon_calendar(pool)
        self.assertEqual(len(shared),0)
        self.assertEqual(filtered[TAXA[0]],{})
        self.assertEqual(filtered[TAXA[1]],{})
        self.assertGreaterEqual(n,2)
        a,b=per_tag_effects(filtered[TAXA[0]])
        self.assertEqual(bootstrap_source(b,12)["status"],"UNAVAILABLE")
        self.assertEqual(between_taxa_effect(a,b,42)["schedules"]["original"]["status"],
                         "INSUFFICIENT_COMMON_DATE_TAGGED_ANIMALS")

    def test_common_date_filter_ignores_response_signs(self):
        p1=synthetic()
        p2=synthetic()
        for sp,bytag in p2.items():
            for tag,years in bytag.items():
                for yr,daily in years.items():
                    for d in daily:
                        for sch in SCHEDULES:
                            daily[d][sch]=-1000*daily[d][sch]
        r1,dates1,_=matched_two_taxon_calendar(p1)
        r2,dates2,_=matched_two_taxon_calendar(p2)
        self.assertEqual(dates1,dates2)
        self.assertEqual({sp:set(r1[sp]) for sp in TAXA},
                         {sp:set(r2[sp]) for sp in TAXA})


if __name__=="__main__":
    unittest.main()
