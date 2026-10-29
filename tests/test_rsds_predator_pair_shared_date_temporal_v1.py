"""Known-truth original published coyote/bobcat shared-calendar q and daylight tests."""
from __future__ import annotations
import json
import unittest
from datetime import date,timedelta
from copy import deepcopy

from scripts.run_rsds_predator_pair_shared_date_temporal_v1 import (
    TAXA,SCHEDULES,frozen,same_native_calendar,tag_effects,
    summarize_tags,pair_contrast,calendar_date_slope
)


def pool(n1=7,n2=8,days=80,negative=False):
    """Private test-only tag-keyed shared astronomy dates, no real GPS fixes."""
    c={}
    origin=date(2020,2,1)
    for i,taxon in enumerate(TAXA):
        members={}
        for tag in range(n1 if i==0 else n2):
            all_days={}
            for k in range(days):
                d=origin+timedelta(days=k)
                q=(.1 if i==0 else -.03)+.002*k+tag*.001
                if negative:q=-q
                all_days[d]={"original":q,"even":q-.01,"odd":q+.01}
            members[f"private_{i}_{tag}"]={2020:all_days}
        c[taxon]=members
    return c


class CoyotesBobcatsNewOriginalSourceCalendarTests(unittest.TestCase):
    def test_two_public_source_UUID_MD5_license_and_species_are_fixed_before_new_q(self):
        c,earlier=frozen()
        self.assertEqual([x["region"] for x in c["source_studies"]],["Washington","Illinois"])
        self.assertEqual([x["doi"] for x in c["source_studies"]],
                         ["10.5441/001/1.gm93267b","10.5441/001/1.649"])
        self.assertEqual(c["source_studies"][1]["native_taxa"],
                         ["Canis latrans","Lynx rufus","Odocoileus virginianus"])
        self.assertTrue(all(earlier[x["package_uuid"]]["source_state"]==
                            "original_MD5_verified_source_analyzed"
                            for x in c["source_studies"]))
        self.assertEqual(len({x["study_id"] for x in c["source_studies"]}),2)

    def test_same_original_calendar_days_selected_without_q_sign(self):
        x=pool()
        a,days,n=same_native_calendar(x[TAXA[0]],x[TAXA[1]],SCHEDULES)
        self.assertEqual(len(a[TAXA[0]]),7)
        self.assertEqual(len(a[TAXA[1]]),8)
        self.assertEqual(len(days),80)
        for t in TAXA:
            self.assertEqual({(year,d) for z in a[t].values()
                             for year,ys in z.items() for d in ys},days)
        y=pool(negative=True)
        b,ds,_=same_native_calendar(y[TAXA[0]],y[TAXA[1]],SCHEDULES)
        self.assertEqual({t:set(z) for t,z in a.items()},
                         {t:set(z) for t,z in b.items()})
        self.assertEqual(ds,days)

    def test_native_source_daily_q_and_same_day_photo_astro_positive_contrast(self):
        # Five or more year-month blocks are required for the date-block slope.
        x=pool(days=160)
        shared,dates,_=same_native_calendar(x[TAXA[0]],x[TAXA[1]],SCHEDULES)
        scores={}
        for sp in TAXA:
            q,g=tag_effects(shared[sp],SCHEDULES)
            scores[sp]=(q,g)
            self.assertEqual(len(q),7 if sp==TAXA[0] else 8)
            self.assertEqual(len(g),len(q))
        a,b=scores[TAXA[0]][0],scores[TAXA[1]][0]
        pair=pair_contrast(a,b,SCHEDULES,123)
        for axis in SCHEDULES:
            self.assertEqual(pair[axis]["nominal_sign"],"POSITIVE")
            self.assertAlmostEqual(pair[axis]["coyote_minus_bobcat_mean"],.1295,places=12)
        day=calendar_date_slope(shared,dates,SCHEDULES,999)
        self.assertEqual(day["coobserved_source_calendar_dates"],160)
        self.assertEqual(set(day["slopes"]),set(SCHEDULES))
        self.assertTrue(all(abs(v["same_calendar_coyote_minus_bobcat_q_daylength_slope"])<1e-12
                            for v in day["slopes"].values()))

    def test_original_source_tag_bootstrap_needs_five_individuals(self):
        d=pool(n1=3,n2=7)
        chosen,days,_=same_native_calendar(d[TAXA[0]],d[TAXA[1]],SCHEDULES)
        a,_=tag_effects(chosen[TAXA[0]],SCHEDULES)
        b,_=tag_effects(chosen[TAXA[1]],SCHEDULES)
        q=summarize_tags(a,SCHEDULES,123)
        self.assertTrue(all(v["nominal_sign"]=="UNAVAILABLE" and
                            v["source_animal_mean"] is None for v in q.values()))
        comp=pair_contrast(a,b,SCHEDULES,123)
        self.assertTrue(all(v["nominal_sign"]=="UNAVAILABLE" and
                            v["coyote_minus_bobcat_mean"] is None for v in comp.values()))

    def test_no_common_dates_cannot_be_mistaken_for_zero_ecological_response(self):
        x=pool(n1=7,n2=8,days=80)
        y=deepcopy(x[TAXA[1]])
        for z in y.values():
            for year in list(z):
                z[year]={d+timedelta(days=350):v for d,v in z[year].items()}
        result,dates,iters=same_native_calendar(x[TAXA[0]],y,SCHEDULES)
        self.assertEqual(len(dates),0)
        self.assertEqual(len(result[TAXA[0]]),0)
        self.assertEqual(len(result[TAXA[1]]),0)


if __name__=="__main__":
    unittest.main()
