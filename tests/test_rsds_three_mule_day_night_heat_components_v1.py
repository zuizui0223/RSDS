"""Known-truth tests for Q = NIGHT - DAY original GPS speed regressions."""
from __future__ import annotations
from datetime import datetime,date,timedelta
import math
from unittest.mock import patch
import unittest

from scripts.movebank_native_ecology_20261019 import daylength
from scripts import rsds_three_mule_day_night_heat_components_v1 as m
from scripts import rsds_three_mule_within_month_joint_weather_v1 as old
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as original
from scripts import aggregate_rsds_three_mule_day_night_heat_components_v1 as agg


def synthetic_tag(tag):
    """Original private DAY/NIGHT component log speeds with three month intercepts."""
    components=[];weather={};prior=[]
    for month in (2,4,10):
        for k in range(15):
            d=date(2020,month,3)+timedelta(days=k)
            L=daylength(42,d)
            Tmax=10+2.1*math.sin(k*1.17)+.22*k+.17*month
            water=.55+.22*math.cos(k*.91)+.035*k+.08*month
            day=.63-.004*Tmax+.04*L+.018*water+.012*month
            night=.73+.007*Tmax+.03*L+.020*water-.009*month
            q=night-day
            row=(d,L,q,(40,-100),(2020,month),day,night)
            components.append(row)
            prior.append(row[:5])
            weather[d]=(Tmax-1.5,Tmax,water)
    return (tag,2020,components),(tag,2020,prior),weather


class OriginalMuleDayNightDecomposition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c,cls.studies,cls.old,cls.prior=m.frozen()

    def test_frozen_three_studies_previous_q_temperature_first(self):
        self.assertEqual(len(self.studies),3)
        self.assertEqual(len(self.prior),3)
        self.assertEqual([x["source_key"] for x in self.studies],
                         ["mule_old_utah","mule_old_wyoming","mule_new_sublette"])
        self.assertEqual([z["source_status"] for z in self.prior],
                         ["SOURCE_MONTH_JOINT_ESTIMABLE"]*3)
        self.assertEqual(self.c["previous_12_ERA5_raw_SHA_source_status"],"HOLD")

    def test_one_tag_year_q_day_night_temperature_identity(self):
        comp,prior,weather=synthetic_tag("T0")
        rows=[(d,weather[d][1],weather[d][2],L,q,day,night)
              for d,L,q,cell,ym,day,night in comp[2]]
        z=m.matched_year_component_effect(rows)
        self.assertIsNotNone(z)
        self.assertAlmostEqual(z["month_day_tmax"],-.004,places=6)
        self.assertAlmostEqual(z["month_night_tmax"],.007,places=6)
        self.assertAlmostEqual(z["month_q_tmax"],.011,places=6)
        self.assertAlmostEqual(z["month_night_tmax"]-z["month_day_tmax"],
                               z["month_q_tmax"],places=12)
        self.assertAlmostEqual(z["raw_night_tmax"]-z["raw_day_tmax"],
                               z["raw_q_tmax"],places=12)

    def test_original_eligible_tag_and_date_source_q_components_exact(self):
        events={}
        start=date(2020,3,1)
        for ix in range(15):
            d=start+timedelta(days=ix)
            for span,step in ((0,.002),(9,.001)):
                for h in range(span,span+4):
                    t=datetime(d.year,d.month,d.day,h,0,0)
                    lat=20.+(0.002*ix)+(h-span)*step
                    events[t]=(t,lat,0.)
        original_frame,geography=original.source_days({"T0":events})
        comp=m.original_source_components({"T0":events})
        result=m.verify_original_daily_q_reconstruction(original_frame,comp)
        self.assertGreaterEqual(result["original_q_reconstructed_tag_solar_dates_not_independent"],10)
        self.assertTrue(result["original_source_tag_day_grid_yearmonth_coverage_exactly_identical"])

    def test_shared_tag_resampling_matches_previous_source_q_bootstrap_exact(self):
        original_rows=[];comp_rows=[];weather={}
        for j in range(6):
            comp,original_i,values=synthetic_tag("tag"+str(j))
            original_rows.append(original_i)
            comp_rows.append(comp)
            weather.update(values)
        with patch.object(old,"climate_at_gps_day",side_effect=lambda _,d:weather[d]):
            previous=old.per_source_same_sample_models(original_rows,{(40,-100):{}},20261116)
        fake={"new_month_same_tag_model":previous}
        with patch.object(m,"climate_at_gps_day",side_effect=lambda _,d:weather[d]):
            res=m.summarize_original_tag_components(comp_rows,{(40,-100):{}},fake,20261116)
        self.assertEqual(res["source_status"],"SOURCE_COMPONENTS_ESTIMABLE")
        self.assertEqual(res["eligible_original_tagged_animals"],6)
        self.assertTrue(res["source_q_temperature_exact_mean_and_95_CI_replayed"])
        self.assertTrue(res["q_night_minus_day_algebra_exact_in_all_tag_years_and_bootstrap"])
        self.assertAlmostEqual(res["components"]["month_q_tmax"]["estimate"],.011,places=6)
        self.assertEqual(res["components"]["month_day_tmax"]["nominal_sign"],"NEGATIVE")
        self.assertEqual(res["components"]["month_night_tmax"]["nominal_sign"],"POSITIVE")

    def test_original_source_failure_is_HOLD_not_zero(self):
        with patch.object(m,"original_files",side_effect=OSError("private GPS path")) as ff:
            result=m.process(0)
            ff.assert_called_once()
        self.assertEqual(result["status"],"SOURCE_HOLD")
        self.assertIsNone(result["source_component_summary"])
        self.assertFalse(result["source_original_GPS_MD5_bytes_verified"])
        self.assertEqual(result["hold_type"],"OSError")
        self.assertNotIn("private GPS path",str(result))

    def test_all_three_source_HOLDs_cannot_count_as_temperature_signs(self):
        a=[]
        for i,source in enumerate(self.studies):
            a.append({
                "version":m.SOURCE_VERSION,
                "source_index":i,
                "source_key":source["source_key"],
                "original_study_DOI":source["doi"],
                "original_study_id":str(source["study_id"]),
                "original_registered_taxa":source["registered_taxa"],
                "status":"SOURCE_HOLD",
                "source_component_summary":None,
                "source_weather_grid_count":0,"source_weather_grid_HOLD_count":0,
                "no_original_GPS_tag_ids_positions_exact_dates_climate_grid_component_coefficients_exported":True,
                "old_48_irradiance_heading_direction_signed_test_NOT_REPLICATED":True,
                "old_12_ERA5_raw_SHA_HOLD":True,
                "hold_stage":"ORIGINAL_PROVIDER_MD5_UUID_BYTES","hold_type":"OSError",
                "original_full_native_q_daylength_replayed":False,
                "original_source_daily_q_with_day_night_components_reconstructed":False,
                "same_20261116_daily_ERA5_physical_SHA_replayed":False,
                "original_six_weather_source_model_means_replayed":False,
                "original_q_Tmax_month_first_bootstrap_replayed":False,
            })
        result=agg.aggregate(a)
        self.assertEqual(result["registered_original_source_DOI_studies"],3)
        self.assertEqual(result["original_source_status_counts"],{"SOURCE_HOLD":3})
        self.assertEqual(result["original_tag_source_component_nominal_CI_sign_counts"]["source_Q_Tmax_nominal_positive"],0)
        with self.assertRaisesRegex(ValueError,"All 3"):
            agg.aggregate(a[:2])

if __name__=="__main__":unittest.main()
