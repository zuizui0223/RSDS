"""Known truth: identical original tagged animal/calendar date q ecology over 3 GPS schedules."""
from __future__ import annotations

import json
import math
import unittest
from datetime import date, timedelta
from unittest.mock import patch

from scripts.run_rsds_triple_matched_ecological_q_v1 import (
    Q_KEYS, SLOPE_KEYS, freeze, private_tag_q_and_daylength,
    paired_bootstrap, classify_stability, compare_previous_source_delta
)
from scripts.run_rsds_gps_three_phase_exact_dates_v2 import (
    triple_from_private_original_source, aggregate_private_triple
)


def artificial(n=7,days=24,change=None):
    full,even,odd={},{},{}
    start=date(2020,4,1)
    for j in range(n):
        tag="PRIVATE_TAG_"+str(j)
        tables=[{},{},{}]
        for i in range(days):
            d=start+timedelta(days=i)
            L=10.0+.075*i
            base=.025+.012*L+.002*j
            qs=(base,base+.003+j*.0002,base-.004+j*.0001)
            for k in range(3):
                tables[k][d]=(qs[k],L+(100 if k==1 else -100 if k==2 else 0))
        if change is not None and j in change:
            k,deln=change[j]
            tables[k]={d:val for d,val in list(tables[k].items())[:deln]}
        for target,table in zip((full,even,odd),tables):
            target[tag]={2020:table}
    return full,even,odd


class ThreeScheduleEcologicalQKnownTruth(unittest.TestCase):
    def test_frozen_8_studies_12_original_taxa_and_old_negative(self):
        c,old=freeze()
        self.assertEqual(old["fixed_original_source_studies"],8)
        self.assertEqual(old["native_source_study_taxon_rows"],12)
        self.assertEqual(c["version"],"RSDS_SAME_TAG_TRIPLE_DATE_ECOLOGICAL_Q_AND_DAYLENGTH_INVARIANCE_V1")

    def test_absolute_q_levels_and_slope_all_share_same_original_days(self):
        original,even,odd=artificial()
        q,s,den=private_tag_q_and_daylength(original,even,odd)
        self.assertEqual(len(q),7)
        self.assertEqual(len(s),7)
        self.assertEqual(den["matched_original_tagged_animal_calendar_days"],7*24)
        self.assertEqual(den["source_tag_years_with_10_or_more_identical_dates"],7)
        group=paired_bootstrap(q,Q_KEYS)
        slopes=paired_bootstrap(s,SLOPE_KEYS)
        self.assertEqual(classify_stability(group,Q_KEYS),"THREE_SCHEDULES_NOMINAL_POSITIVE")
        self.assertEqual(classify_stability(slopes,SLOPE_KEYS),"THREE_SCHEDULES_NOMINAL_POSITIVE")
        self.assertAlmostEqual(group["q_even"]["group_original_tag_mean"]-
                               group["q_odd"]["group_original_tag_mean"],.0073,places=12)
        self.assertAlmostEqual(slopes["daylength_even"]["group_original_tag_mean"],.012,places=12)
        self.assertAlmostEqual(slopes["daylength_odd"]["group_original_tag_mean"],.012,places=12)
        self.assertEqual(group,paired_bootstrap(q,Q_KEYS))
        self.assertNotIn("PRIVATE_TAG_",json.dumps(group))

    def test_group_absolute_levels_reconstruct_previous_original_source_phase_deltas(self):
        original,even,odd=artificial()
        q,s,_=private_tag_q_and_daylength(original,even,odd)
        pair,_=triple_from_private_original_source(original,even,odd)
        prior={"new_strict_triple_original_tag_effects":aggregate_private_triple(pair)}
        compare_previous_source_delta(paired_bootstrap(q,Q_KEYS),
                                      paired_bootstrap(s,SLOPE_KEYS),prior)
        bad=json.loads(json.dumps(prior))
        bad["new_strict_triple_original_tag_effects"]["q_even_minus_odd"][
            "group_source_mean_delta"]+=.001
        with self.assertRaisesRegex(ValueError,"fail to reconstruct"):
            compare_previous_source_delta(paired_bootstrap(q,Q_KEYS),
                                          paired_bootstrap(s,SLOPE_KEYS),bad)

    def test_disjoint_three_dates_or_missing_tag_does_not_make_a_zero(self):
        f,e,o=artificial(n=7,change={0:(2,6),1:(1,5),2:(2,4)})
        q,s,den=private_tag_q_and_daylength(f,e,o)
        self.assertEqual(len(q),4)
        a=paired_bootstrap(q,Q_KEYS)
        self.assertTrue(all(x["status"]=="UNAVAILABLE" and
                            x["group_original_tag_mean"] is None and
                            x["nominal_source_tag_bootstrap95"] is None
                            for x in a.values()))
        self.assertEqual(classify_stability(a,Q_KEYS),"UNAVAILABLE")

    def test_support_on_different_days_is_rejected_from_original_delta(self):
        f,e,o=artificial()
        q,s,den=private_tag_q_and_daylength(f,e,o)
        group=paired_bootstrap(q,Q_KEYS)
        slopes=paired_bootstrap(s,SLOPE_KEYS)
        _,_ = triple_from_private_original_source(f,e,o)
        fake={"new_strict_triple_original_tag_effects":{
            key:{"source_original_same_tags_and_triple_dates":7,
                 "group_source_mean_delta":0}
            for key in ("q_even_minus_original","q_odd_minus_original","q_even_minus_odd",
                        "daylength_even_minus_original","daylength_odd_minus_original",
                        "daylength_even_minus_odd")}}
        with self.assertRaisesRegex(ValueError,"fail to reconstruct"):
            compare_previous_source_delta(group,slopes,fake)

    def test_opposite_signs_are_not_reported_as_stable(self):
        d={
            "q_original":{"status":"POSITIVE"},
            "q_even":{"status":"POSITIVE"},
            "q_odd":{"status":"NEGATIVE"}
        }
        self.assertEqual(classify_stability(d,Q_KEYS),"OPPOSITE_NOMINAL_SOURCE_TAG_SIGN")
        d["q_odd"]["status"]="INCLUDES_ZERO"
        self.assertEqual(classify_stability(d,Q_KEYS),
                         "UNRESOLVED_AT_LEAST_ONE_ZERO_SPANNING_INTERVAL")


if __name__=="__main__":
    unittest.main()
