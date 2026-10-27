"""Known-truth tests of source-original GPS duty cycle, matched dates and zero/HOLD boundaries."""
from __future__ import annotations
import math
import unittest
from datetime import datetime, timedelta
from unittest.mock import patch

from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    PAIRS, verify_frame, segment, per_tag_year_days, matched_per_tag,
    summary, analyze_taxon, source_baseline_check, run_source
)
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon


def synthetic_animals(n=6,days=90):
    """Only synthetic telemetry; no real animal GPS IDs or locations used."""
    out={}
    start=datetime(2020,1,1,0)
    for a in range(n):
        f={}
        la=40.0+0.001*a
        lon=-105.0
        for k in range(days*24):
            t=start+timedelta(hours=k)
            # Different observed straight-line speeds at daytime and nighttime.
            # Stable GPS timestamps at exact hour boundaries, no missing data.
            la+=.003 if t.hour<7 or t.hour>=18 else .012
            f[t]=(t,la,lon)
        out[f"synthetic{a}"]=f
    return {"Synthetic.test":out}


class FixedOriginalGPSCadenceMatchedSupport(unittest.TestCase):
    def test_original_eight_study_and_source_holds_are_frozen(self):
        old,new,frame,files=verify_frame()
        self.assertEqual(old["fixed_original_study_count"],2)
        self.assertEqual(new["frozen_source_original_study_count"],6)
        self.assertEqual(len(frame),6)
        self.assertEqual(len(files),6)
        self.assertEqual(sum(x["source_metadata_status"]=="SOURCE_METADATA_HOLD" for x in files),1)

    def test_hour_and_both_UTC_two_hour_phases_have_no_interpolation(self):
        records=synthetic_animals(6,90)["Synthetic.test"]
        full,base=per_tag_year_days(records,"ORIGINAL_1H")
        self.assertEqual(base["tags_with_some_q_observable"],6)
        for mode in PAIRS:
            y,stats=per_tag_year_days(records,mode)
            self.assertGreater(stats["qualified_same_solar_period_GPS_intervals"],0)
            self.assertLess(stats["candidate_adjacent_GPS_intervals"],
                            base["candidate_adjacent_GPS_intervals"])
            q,s,support=matched_per_tag(full,y)
            self.assertEqual(support["source_tagged_animals_with_ge10_matching_days"],6)
            self.assertEqual(len(q),6)
            self.assertTrue(all(math.isfinite(v) for v in q.values()))
            self.assertGreater(support["original_tagged_days_identical_between_1h_and_2h"],60)

    def test_original_full_hour_replays_native_source_q_and_daylength(self):
        d=synthetic_animals()
        obj,_=one_taxon("Synthetic.test",d["Synthetic.test"])
        fixed={"Synthetic.test":{
            "day_night_speed_contrast":obj["day_night_speed_contrast"],
            "daylength_association":obj["daylength_association"],
        }}
        self.assertTrue(source_baseline_check(d,fixed))
        altered={"Synthetic.test":{
            "day_night_speed_contrast":dict(obj["day_night_speed_contrast"]),
            "daylength_association":obj["daylength_association"],
        }}
        altered["Synthetic.test"]["day_night_speed_contrast"]["mean"]+=.2
        with self.assertRaises(ValueError):
            source_baseline_check(d,altered)

    def test_same_days_only_not_introduce_fictitious_animals(self):
        base={"tagA":{2020:{datetime(2020,1,1).date()+timedelta(days=5*i):
                            (-.1,.4+i*.12) for i in range(11)}},
              "tagB":{2020:{datetime(2020,1,1).date()+timedelta(days=i):
                            (-.1,9.0+i*.04) for i in range(11)}}}
        thinning={"tagA":{2020:{d:(q+.02,L) for d,(q,L) in base["tagA"][2020].items()}},
                  "tagOther":{2020:base["tagB"][2020]}}
        q,s,support=matched_per_tag(base,thinning)
        self.assertEqual(support["source_tagged_animals_with_ge10_matching_days"],1)
        self.assertEqual(set(q),{"tagA"})
        self.assertAlmostEqual(q["tagA"],.02)
        self.assertAlmostEqual(s["tagA"],0,places=12)
        self.assertEqual(summary(q,42)["status"],"INSUFFICIENT_MATCHED_TAGS")
        self.assertIsNone(summary(q,42)["paired_tag_mean_2h_minus_same_day_1h"])

    def test_original_q_and_2h_parity_effect_source_reports_only_group_counts(self):
        d=synthetic_animals()
        r=analyze_taxon("Synthetic.test",d["Synthetic.test"])
        self.assertEqual(set(r["source_original_even_odd_2h_phase_sensitivity"]),set(PAIRS))
        for p in r["source_original_even_odd_2h_phase_sensitivity"].values():
            self.assertEqual(p["matched_q_velocity_proxy_difference"]["status"],
                             "SOURCE_TAGGED_ANIMAL_DESCRIPTIVE")
            self.assertGreaterEqual(p["identical_original_tag_date_support"][
                "source_tagged_animals_with_ge10_matching_days"],5)
            self.assertEqual(len(p["matched_q_velocity_proxy_difference"]["original_tag_resample95"]),2)
        self.assertNotIn("synthetic0",str(r))
        self.assertNotIn('"longitude"',str(r))


if __name__=="__main__":
    unittest.main()
