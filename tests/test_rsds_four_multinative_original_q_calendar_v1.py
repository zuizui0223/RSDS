"""Known-truth source-native shared-calendar GPS q and photoperiod comparisons."""
from __future__ import annotations

import json
import math
import unittest
from datetime import date,timedelta

from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import (
    fixed_sources, common_original_calendar, tag_level_q_and_season,
    effect,pair_native_q,pair_season,source_reference_latitude
)


def synthetic_native(n_tags=6,days=120):
    start=date(2020,3,1)
    sources={}
    for taxon,baseline,slope_day in [
            ("Native species A",.08,.006),
            ("Native species B",-.07,-.002),
            ("Native species C",.025,.0)]:
        track={}
        for k in range(n_tags):
            records={}
            for i in range(days):
                d=start+timedelta(days=i)
                q=baseline+slope_day*i/30+.001*k
                records[d]=(q,12.0)
            track[f"PRIVATE-{taxon}-{k}"]={2020:records}
        sources[taxon]=track
    return sources


class FourNewOriginalMultitaxonEcologyKnownTruth(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen,cls.metadata=fixed_sources()

    def test_four_original_source_file_pinned_before_any_new_GPS_event(self):
        self.assertEqual(len(self.frozen["sources"]),4)
        self.assertEqual(len(self.metadata),4)
        self.assertEqual(sum(len(s["registered_taxa"]) for s in self.frozen["sources"]),10)
        self.assertTrue(all(s["source_status"]=="SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE"
                            for s in self.metadata))
        self.assertTrue(all(s["source_original_file"]["checksum"]["checkSumAlgorithm"]=="MD5"
                            for s in self.metadata))

    def test_all_three_source_taxa_same_dates_not_chosen_by_q_sign(self):
        native=synthetic_native()
        shared,days,rounds=common_original_calendar(native)
        self.assertEqual(len(days),120)
        self.assertTrue(all(len(v)==6 for v in shared.values()))
        self.assertEqual({t:len(v) for t,v in shared.items()},{
            "Native species A":6,"Native species B":6,"Native species C":6
        })
        a,grads=tag_level_q_and_season(shared["Native species A"],43)
        b,_=tag_level_q_and_season(shared["Native species B"],43)
        self.assertEqual(len(a),6)
        self.assertEqual(len(grads),6)
        q=pair_native_q(a,b,20261029)
        self.assertEqual(q["nominal_interval_sign"],"POSITIVE")
        self.assertGreater(q["first_minus_second_q"],.1)
        self.assertNotIn("PRIVATE",json.dumps(q))

    def test_exclusive_dates_yield_no_fabricated_species_niche(self):
        native=synthetic_native(days=40)
        c=native["Native species C"]
        shift={tag:{2020:{d+timedelta(days=200):q
                            for d,q in ys[2020].items()}}
               for tag,ys in c.items()}
        native["Native species C"]=shift
        shared,days,_=common_original_calendar(native)
        self.assertEqual(days,set())
        self.assertTrue(all(not sp for sp in shared.values()))
        q,_=tag_level_q_and_season(shared["Native species A"],43)
        r=effect(q,20261029)
        self.assertIsNone(r["mean"])
        self.assertEqual(r["nominal_interval_sign"],"UNAVAILABLE")

    def test_equatorial_reference_daylength_slope_unsupported_without_dropping_q(self):
        native=synthetic_native()
        share,days,_=common_original_calendar(native)
        q,grad=tag_level_q_and_season(share["Native species A"],0.0)
        self.assertEqual(len(q),6)
        self.assertEqual(len(grad),0)
        self.assertEqual(effect(grad,20261029)["nominal_interval_sign"],"UNAVAILABLE")
        self.assertNotEqual(effect(q,20261029)["nominal_interval_sign"],"UNAVAILABLE")

    def test_common_source_datewise_month_threshold_not_outcome_selected(self):
        native=synthetic_native(days=60)
        grouped,days,_=common_original_calendar(native)
        from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import daily_source_q
        a=daily_source_q(grouped["Native species A"],days)
        b=daily_source_q(grouped["Native species B"],days)
        result=pair_season(a,b,days,43,20261029)
        self.assertEqual(result["source_shared_calendar_days"],60)
        self.assertEqual(result["all_yearmonth_slope_status"],"UNAVAILABLE")
        self.assertEqual(result["within_month_reference_daylength"]["nominal_interval_sign"],"UNAVAILABLE")

    def test_original_study_reference_location_never_serialized(self):
        fake={
            "A":{"x":{0:(None,0.25, 10),1:(None,0.20,20)}},
            "B":{"y":{0:(None,-0.15,15)}}
        }
        lat=source_reference_latitude(fake)
        self.assertTrue(-1<lat<1)
        self.assertNotIn("longitude",str(lat))


if __name__=="__main__":
    unittest.main()
