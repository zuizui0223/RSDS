"""Original q, deterministic equal 3 intervals, exact source-tag support rules."""
from __future__ import annotations
from datetime import datetime,date,timedelta
import math
from unittest.mock import patch
import unittest

from scripts import rsds_three_mule_equal_gps_day_night_opportunity_v1 as m
from scripts import rsds_three_mule_day_night_heat_components_v1 as previous
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as original
from scripts import aggregate_rsds_three_mule_equal_gps_day_night_opportunity_v1 as agg
from scripts.movebank_native_ecology_20261019 import daylength

class ThreeMuleOriginalEqualGPSSampling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.sources,cls.prev,cls.prev_components=m.frozen()

    def test_three_frozen_original_sources_and_old_q_heat(self):
        self.assertEqual(len(self.sources),3)
        self.assertEqual([x["source_key"] for x in self.sources],
                         ["mule_old_utah","mule_old_wyoming","mule_new_sublette"])
        self.assertTrue(all(z["source_status"]=="SOURCE_COMPONENTS_ESTIMABLE"
                            for z in self.prev_components))
        self.assertTrue(self.contract["species_sources_selected_after_seeing_previous_positive_Tmax"])
        self.assertEqual(self.contract["old_48_source_irradiance_heading_external_direction_signed_result"],"NOT_REPLICATED")

    def test_equal_three_from_original_timestamp_order_not_sorted_speed(self):
        self.assertEqual(m.choose_three_by_original_interval_order([9,2,8]),[9,2,8])
        self.assertEqual(m.choose_three_by_original_interval_order([9,2,8,3]),[9,2,3])
        self.assertEqual(m.choose_three_by_original_interval_order([10,4,20,3,8,1,6]),[10,3,6])
        with self.assertRaisesRegex(ValueError,"DAY/NIGHT"):
            m.choose_three_by_original_interval_order([1,2])

    def test_original_q_date_speed_components_exact_reproduction(self):
        obs={}
        for k in range(15):
            d=date(2020,3,1)+timedelta(days=k)
            for partstart,step in ((0,.002),(9,.001)):
                for hour in range(partstart,partstart+5):
                    t=datetime(d.year,d.month,d.day,hour)
                    obs[t]=(t,20+(.002*k)+(hour-partstart)*step,0.)
        prior,_=original.source_days({"T0":obs})
        balanced=m.original_balanced_date_rows({"T0":obs})
        audit=previous.verify_original_daily_q_reconstruction(
            prior,[(tag,year,[r[:7] for r in rows])
                   for tag,year,rows in balanced])
        self.assertGreater(audit["original_q_reconstructed_tag_solar_dates_not_independent"],9)
        for tag,year,rows in balanced:
            for d,L,q,cell,ym,day,night,bd,bn,nd,nn in rows:
                self.assertAlmostEqual(night-day,q,places=10)
                self.assertGreaterEqual(nd,3)
                self.assertGreaterEqual(nn,3)
                self.assertTrue(math.isfinite(bn-bd))

    def test_known_true_balanced_calendar_part_regression(self):
        rows=[]
        for month in (2,4,10):
            for k in range(15):
                d=date(2020,month,3)+timedelta(days=k)
                L=daylength(42,d)
                T=10+2.1*math.sin(k*1.17)+.22*k+.17*month
                water=.55+.22*math.cos(k*.91)+.035*k+.08*month
                day=.52+.002*T+.04*L+.02*water+.03*month
                night=.72+.01*T+.03*L+.04*water-.01*month
                bd=day-.003*T+.02
                bn=night+.004*T-.04
                rows.append((d,T,water,L,night-day,day,night,bd,bn,5,9))
        r=m.balanced_month_model(rows)
        self.assertIsNotNone(r)
        coeff,nDates,imbalance=r
        self.assertEqual(nDates,45)
        self.assertEqual(imbalance,45)
        self.assertAlmostEqual(coeff["old_month_Q_Tmax"],.008,places=7)
        self.assertAlmostEqual(coeff["balanced_month_Q_Tmax"],.015,places=7)
        self.assertAlmostEqual(coeff["balanced_month_DAY_Tmax"],-.001,places=7)
        self.assertAlmostEqual(coeff["balanced_month_NIGHT_Tmax"],.014,places=7)
        self.assertAlmostEqual(coeff["balanced_minus_old_month_Q_Tmax"],.007,places=7)

    def test_error_original_provider_cannot_claim_robustness(self):
        with patch.object(m,"original_files",side_effect=OSError("private location token")) as v:
            result=m.process(0)
            v.assert_called_once()
        self.assertEqual(result["status"],"SOURCE_HOLD")
        self.assertFalse(result["provider_source_MD5_verified"])
        self.assertIsNone(result["result"])
        self.assertEqual(result["hold_stage"],"SOURCE_ORIGINAL_MD5")
        self.assertNotIn("private location token",str(result))

    def test_source_original_all_three_hold_not_substituted(self):
        rows=[]
        for i,s in enumerate(self.sources):
            rows.append({
                "version":m.SOURCE_VERSION,
                "source_index":i,"source_key":s["source_key"],
                "original_study_DOI":s["doi"],
                "original_study_ID":str(s["study_id"]),
                "source_native":"Odocoileus hemionus",
                "original_registered_native_taxa":s["registered_taxa"],
                "status":"SOURCE_HOLD","result":None,
                "hold_stage":"SOURCE_ORIGINAL_MD5","hold_type":"OSError",
                "source_48_irradiance_heading_direction_external_sign_NOT_REPLICATED":True,
                "source_12_old_ERA5_raw_SHA_HOLD":True,
                "no_original_tag_ID_GPS_exact_date_climate_cell_or_individual_speeds_exported":True,
            })
        final=agg.aggregate(rows)
        self.assertEqual(final["original_source_study_DOI_count"],3)
        self.assertEqual(final["original_source_studies_with_equal_step_model"],0)
        self.assertEqual(final["source_status_counts"],{"SOURCE_HOLD":3})
        with self.assertRaisesRegex(ValueError,"3"):
            agg.aggregate(rows[:2])

if __name__=="__main__":
    unittest.main()
