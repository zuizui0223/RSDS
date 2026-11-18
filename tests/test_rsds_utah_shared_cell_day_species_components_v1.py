"""Known truth: same original GPS location AND date, source-cell/month test.

No original provider activity in ordinary tests. Maintain original 95
opportunity denominator, reject map/date manipulation and SOURCE_HOLD with
no published GPS positions/animal IDs.
"""
from __future__ import annotations
from collections import defaultdict
from datetime import date,timedelta
from math import sin,cos
import unittest
from unittest.mock import patch

from scripts import rsds_utah_shared_cell_same_day_mule_bighorn_components_v1 as m
from scripts.movebank_native_ecology_20261019 import daylength


def one_block(cell,month,coef=.002):
    rows=[]
    for k in range(15):
        d=date(2020,month,2)+timedelta(days=k)
        L=daylength(cell[0]+.5,d)
        T=10+.31*k+2.3*sin(1.11*k)+.08*month+.12*cell[0]
        W=.2+.03*k+.19*cos(.87*k)+.02*month
        # Shared source date and SAME cell/week/month weather
        dq=coef*T+.007*L+.005*W+2.0*month
        dd=-.0007*T+.004*L+.005*W+.1*month
        dn=dq+dd
        rows.append((T,L,W,dq,dd,dn))
    return rows


class OriginalUtahCellDayPair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c,cls.utah,cls.legacy,cls.native=m.frozen()

    def test_fixed_utah_original_source_and_doi(self):
        self.assertEqual(len(self.native),4)
        self.assertEqual(self.utah["doi"],"10.5441/001/1.711")
        self.assertEqual(self.c["prior_source_shared_grid_yearmonth_pairs_any_gps_tag"],118)
        self.assertEqual(self.c["prior_source_shared_grid_yearmonth_with_at_least_3_distinct_tags_each_taxon"],95)
        self.assertEqual(self.c["source_first_full_ERA5_grid_count"],28)
        self.assertTrue(self.c["selected_after_prior_source_positive_temperature_outcomes"])
        self.assertEqual(self.c["previous_48_signed_solar_heading_external_verdict"],"NOT_REPLICATED")
        self.assertEqual(self.c["previous_12_ERA5_raw_SHA_source_status"],"HOLD")

    def test_within_source_cell_month_joint_q_night_minus_day_estimator(self):
        blocks=[one_block((42,-113),2),one_block((42,-113),4)]
        result=m.cell_month_demeaned_temperature_slopes(blocks,self.c)
        self.assertIsNotNone(result)
        v,n,b=result
        self.assertEqual(n,30)
        self.assertEqual(b,2)
        self.assertAlmostEqual(v["q"],.002,places=9)
        self.assertAlmostEqual(v["day"],-.0007,places=9)
        self.assertAlmostEqual(v["night"],.0013,places=9)
        self.assertAlmostEqual(v["night"]-v["day"],v["q"],places=12)

    def test_geographic_cell_bootstrap_equal_weight_not_pseudo_day_events(self):
        blocks={}
        for lat in range(40,46):
            cell=(lat,-112)
            blocks[cell]=[one_block(cell,2),one_block(cell,4)]
        out=m.original_source_equal_cell_bootstrap(blocks,self.c)
        self.assertEqual(out["source_status"],"SHARED_ORIGINAL_CELL_DAY_MULE_BIGHORN_EFFECT_ESTIMABLE")
        self.assertEqual(out["original_geographic_cells_with_identifiable_joint_Tmax_L_water"],6)
        self.assertEqual(out["original_common_cellmonth_blocks_in_final_cell_level_regression"],12)
        self.assertEqual(out["original_common_solar_dates_in_final_regressions_not_independent"],180)
        fit=out["source_shared_cell_day_species_components_estimates"]
        self.assertAlmostEqual(fit["q"]["mean"],.002,places=9)
        self.assertAlmostEqual(fit["night"]["mean"]-fit["day"]["mean"],fit["q"]["mean"],places=10)
        self.assertEqual(fit["q"]["nominal_sign"],"POSITIVE")
        self.assertTrue(out["same_geographic_cell_bootstrap_draws_for_q_DAY_NIGHT"])

    def test_underpowered_shared_cell_support_is_a_null_HOLD_not_zero(self):
        blocks={}
        for lat in range(42,46):
            cell=(lat,-113)
            blocks[cell]=[one_block(cell,2),one_block(cell,4)]
        out=m.original_source_equal_cell_bootstrap(blocks,self.c)
        self.assertEqual(out["source_status"],"INSUFFICIENT_SHARED_CELL_MONTH_GEOSPATIAL_SUPPORT")
        self.assertEqual(out["original_geographic_cells_with_identifiable_joint_Tmax_L_water"],4)
        self.assertIsNone(out["source_shared_cell_day_species_components_estimates"]["q"]["mean"])
        self.assertEqual(out["source_shared_cell_day_species_components_estimates"]["q"]["nominal_sign"],"UNAVAILABLE")

    def test_only_a_single_source_month_is_insufficient(self):
        self.assertIsNone(m.cell_month_demeaned_temperature_slopes(
            [one_block((42,-113),2)],self.c))

    def test_source_grid_month_original_tag_count_and_shared_dates(self):
        # Synthetic original GPS accepted q-only triples, with three tags
        # per taxon in each of 2 cell/month; weather already in final daily units.
        cell=(42,-113)
        source={m.MULE:defaultdict(dict),m.SHEEP:defaultdict(dict)}
        weather={cell:{}}
        for i in range(150):
            d=date(2019,12,1)+timedelta(days=i)
            weather[cell][d]=(15.,17.,1.,.2)
        for month in (2,4):
            ym=(2020,month);key=(cell,ym)
            for taxon in (m.MULE,m.SHEEP):
                days={}
                for k in range(10):
                    d=date(2020,month,2)+timedelta(days=k)
                    days[d]=[(f"{taxon}{j}",.3,.1,.4) for j in range(3)]
                source[taxon][key]=days
        c=dict(self.c)
        c["prior_source_shared_grid_yearmonth_pairs_any_gps_tag"]=2
        c["prior_source_shared_grid_yearmonth_with_at_least_3_distinct_tags_each_taxon"]=2
        with patch.object(m,"climate_at_gps_day",
                          side_effect=lambda _,d:(15.,17.+2*sin(d.day),.2)):
            summaries,blocks=m.source_shared_day_blocks(source,weather,c)
        self.assertEqual(summaries["source_both_native_original_shared_cell_yearmonth_total"],2)
        self.assertEqual(summaries["source_both_native_original_at_least_three_tags_each_shared_cell_yearmonth"],2)
        self.assertEqual(summaries["original_exact_shared_tag_calendar_dates_across_all_candidate_cell_months"],20)
        self.assertEqual(summaries["source_cell_yearmonth_with_min_same_original_GPS_dates"],2)
        self.assertEqual(len(blocks),1)

    def test_original_provider_failure_is_private_SOURCE_HOLD(self):
        with patch.object(m,"original_files",side_effect=OSError("a private original URL")):
            r=m.process()
        self.assertEqual(r["status"],"SOURCE_HOLD")
        self.assertEqual(r["source_hold_type"],"OSError")
        self.assertFalse(r["original_published_source_GPS_MD5_verified"])
        self.assertIsNone(r["original_source_shared_cell_day_block_model"])
        self.assertNotIn("a private original URL",str(r))


if __name__=="__main__":
    unittest.main()
