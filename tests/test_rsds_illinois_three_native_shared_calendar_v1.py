"""Known-truth three-native Illinois original study day-matched q ecological tests."""
from __future__ import annotations
import json
import unittest
from copy import deepcopy
from datetime import date,timedelta

from scripts.run_rsds_illinois_three_native_shared_calendar_v1 import (
    TAXA,PAIRS,freeze,triple_common_calendar,tag_effects,pair_delta,datewise_pair_slope
)


def synthetic(n=(8,7,10),n_days=155):
    result={}
    start=date(2020,2,1)
    for k,taxon in enumerate(TAXA):
        original={}
        for i in range(n[k]):
            data={}
            for j in range(n_days):
                day=start+timedelta(days=j)
                # Known source true ordering coyote > bobcat > deer.
                level=[.20,.09,-.02][k]+j*.0001+i*.0002
                data[day]={"original":level}
            original["synthetic_private_"+str(k)+"_"+str(i)]={2020:data}
        result[taxon]=original
    return result


class OriginalIllinoisThreeNativeCalendarSourceTest(unittest.TestCase):
    def test_prior_source_first_result_and_all_three_native_taxa_are_frozen(self):
        spec,first,solar=freeze()
        self.assertEqual(spec["source_original_study"]["doi"],"10.5441/001/1.649")
        self.assertEqual(set(first["source_native_taxa"]),set(TAXA))
        self.assertEqual(set(x["taxon"] for x in solar["source_taxa"]),set(TAXA))
        self.assertTrue(all(first["source_native_taxa"][t]["q"]["mean"] is not None for t in TAXA))

    def test_identical_source_calendar_day_support_selected_without_looking_at_q_sign(self):
        p=synthetic()
        data,days,rounds=triple_common_calendar(p)
        self.assertEqual(len(days),155)
        self.assertEqual({t:len(data[t]) for t in TAXA},
                         dict(zip(TAXA,[8,7,10])))
        reversed_sign=deepcopy(p)
        for sp in TAXA:
            for tagged in reversed_sign[sp].values():
                for year in tagged:
                    tagged[year]={d:{"original":-v["original"]} for d,v in tagged[year].items()}
        s,d2,_=triple_common_calendar(reversed_sign)
        self.assertEqual(d2,days)
        self.assertEqual({k:set(v) for k,v in s.items()},
                         {k:set(v) for k,v in data.items()})

    def test_three_native_source_positive_pairwise_ranking(self):
        p=synthetic()
        shared,days,_=triple_common_calendar(p)
        all_scores={}
        for sp in TAXA:
            q,s=tag_effects(shared[sp],("original",))
            all_scores[sp]=(q,s)
        for i,(a,b) in enumerate(PAIRS):
            z=pair_delta(all_scores[a][0],all_scores[b][0],20261029+i)
            self.assertEqual(z["status"],"POSITIVE")
            self.assertGreater(z["first_minus_second_mean"],.10)
            self.assertGreater(z["original_tag_bootstrap95"][0],0)
            pred=datewise_pair_slope(shared[a],shared[b],days,20261029+i)
            self.assertEqual(pred["source_calendar_days"],155)
            self.assertIsNotNone(pred["slope"])
            self.assertAlmostEqual(pred["slope"],0,places=10)

    def test_single_native_source_without_tag_or_qualified_date_is_hold_not_effect_zero(self):
        p=synthetic(n=(8,4,10))
        shared,dates,_=triple_common_calendar(p)
        a,b=tag_effects(shared[TAXA[0]],("original",))
        c,d=tag_effects(shared[TAXA[1]],("original",))
        contrast=pair_delta(a,c,20261029)
        self.assertEqual(contrast["status"],"UNAVAILABLE")
        self.assertIsNone(contrast["first_minus_second_mean"])
        self.assertIsNone(contrast["original_tag_bootstrap95"])
        p=synthetic()
        for tag in p[TAXA[2]].values():
            tag[2020]={d+timedelta(days=330):v for d,v in tag[2020].items()}
        _,shared_dates,_=triple_common_calendar(p)
        self.assertEqual(len(shared_dates),0)


if __name__=="__main__":
    unittest.main()
