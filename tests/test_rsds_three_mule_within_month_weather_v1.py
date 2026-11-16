"""Known truth and fail-closed tests: same-tag month weather partialing."""
import math
import unittest
from datetime import date,timedelta
from unittest.mock import patch

from scripts.movebank_native_ecology_20261019 import daylength
from scripts import rsds_three_mule_within_month_joint_weather_v1 as m
from scripts import aggregate_rsds_three_mule_within_month_weather_v1 as agg

def synthetic_weather_tag(tag="T0"):
    """Three month groups, linear model q=0.13 Tmax+0.08 L-0.06 water+month FE."""
    days=[]
    climate={}
    for mo in (2,4,10):
        for i in range(15):
            d=date(2020,mo,3)+timedelta(days=i)
            L=daylength(42.0,d)
            T=10.0+2.1*math.sin(i*1.17)+.22*i+.17*mo
            water=.55+.22*math.cos(i*.91)+.035*i+.08*mo
            q=.13*T+.08*L-.06*water+2.5*mo
            days.append((d,L,q,(40,-100),(2020,mo)))
            climate[d]=(T-1.5,T,water)
    return (tag,2020,days),climate

class ThreeNativeMonthWeatherTest(unittest.TestCase):
    def test_freeze_source_IDs_and_old_positive_sun_and_weather_completeness(self):
        c,sources,prior=m.freeze()
        self.assertEqual(len(sources),3)
        self.assertEqual([s["source_key"] for s in sources],
                         ["mule_old_utah","mule_old_wyoming","mule_new_sublette"])
        self.assertEqual(len(prior),3)
        self.assertTrue(all(p["Tmax_given_daylength_water_nominal_CI_signs_by_fully_covered_source"]["mule_old_utah"]=="POSITIVE" for p in []))
        self.assertTrue(all(len(z["current_provider_daily_physical_fingerprint"])==64 for z in prior))
        self.assertTrue(c["species_selection_after_prior_positive_temperature_findings"])
        self.assertEqual(c["historic_12_ERA5_raw_SHA_source_status"],"HOLD")

    def test_known_linear_within_month_joint_three_predictors(self):
        (tag,year,rows),weather=synthetic_weather_tag()
        data=[(d,weather[d][1],weather[d][2],L,q) for d,L,q,cell,ym in rows]
        fit=m.within_year_month(data)
        self.assertTrue(fit["supported"],str(fit))
        self.assertEqual(fit["eligible_month_blocks"],3)
        self.assertEqual(fit["eligible_original_dates"],45)
        self.assertAlmostEqual(fit["estimates"]["month_Tmax_given_L_water"],.13,places=7)
        self.assertAlmostEqual(fit["estimates"]["month_L_given_Tmax_water"],.08,places=7)
        self.assertAlmostEqual(fit["estimates"]["month_water_given_Tmax_L"],-.06,places=7)

    def test_six_original_tags_share_9_estimand_support(self):
        recs={};years=[]
        for n in range(6):
            (tag,year,rows),weather=synthetic_weather_tag(f"T{n}")
            years.append((tag,year,rows))
            recs.update(weather)
        with patch.object(m,"climate_at_gps_day",side_effect=lambda _,d:recs[d]):
            result=m.per_source_same_sample_models(years,{(40,-100):{}},20261116)
        self.assertEqual(result["eligible_original_tagged_animals_joint_all_three"],6)
        self.assertEqual(result["original_tag_qualified_month_joint_years"],6)
        self.assertEqual(len(result["month_partial_joint_estimands"]),9)
        self.assertTrue(result["all_nine_estimands_share_identical_original_tag_year_month_date_support"])
        self.assertEqual(result["month_partial_joint_estimands"]["month_Tmax_given_L_water"]["nominal_CI_sign"],"POSITIVE")
        self.assertEqual(result["month_partial_joint_estimands"]["month_water_given_Tmax_L"]["nominal_CI_sign"],"NEGATIVE")
        self.assertAlmostEqual(result["month_partial_joint_estimands"]["month_Tmax_given_L_water"]["estimate"],.13,places=7)

    def test_single_month_and_low_tag_support_do_not_create_zero(self):
        (_,_,rows),weather=synthetic_weather_tag()
        data=[(d,weather[d][1],weather[d][2],L,q)
              for d,L,q,cell,ym in rows if d.month==2]
        z=m.within_year_month(data)
        self.assertFalse(z["supported"])
        self.assertIsNone(z["estimates"])
        records=[synthetic_weather_tag(f"T{n}")[0] for n in range(4)]
        with patch.object(m,"climate_at_gps_day",side_effect=lambda _,d:weather[d]):
            res=m.per_source_same_sample_models(records,{(40,-100):{}},20261116)
        self.assertEqual(res["month_partial_joint_estimands"]["month_Tmax_given_L_water"]["status"],"INSUFFICIENT_ORIGINAL_TAGS")
        self.assertIsNone(res["month_partial_joint_estimands"]["month_Tmax_given_L_water"]["estimate"])

    def test_original_provider_failure_cannot_promote_biology(self):
        with patch.object(m,"original_files",side_effect=OSError("secret original URL")) as mocked:
            r=m.run(0)
            mocked.assert_called_once()
        self.assertEqual(r["status"],"SOURCE_HOLD")
        self.assertEqual(r["hold_stage"],"ORIGINAL_SOURCE_MD5")
        self.assertIsNone(r["source_month_joint"])
        self.assertFalse(r["original_weather_version_exact_20261116_physical_sha_replayed"])
        self.assertNotIn("secret original URL",str(r))

    def test_source_first_weather_missing_cannot_be_counted_as_positive(self):
        c,s,p=m.freeze()
        row=[]
        for i,source in enumerate(s):
            row.append({
                "version":m.SOURCE_VERSION,"source_index":i,
                "source_key":source["source_key"],"source_original_DOI":source["doi"],
                "source_original_study_id":str(source["study_id"]),
                "original_source_registered_native_taxa":source["registered_taxa"],
                "status":"SOURCE_HOLD","hold_stage":"SAME_VERSION_ALL_SOURCE_ERA5_DAILY_PHYSICAL_SHA",
                "hold_type":"ValueError",
                "source_month_joint":None,"provider_original_MD5_verified":True,
                "original_native_q_and_daylength_mean_n_replayed_before_ERA5":True,
                "original_weather_version_exact_20261116_physical_sha_replayed":False,
                "prior_six_climate_model_mean_n_replayed":False,
                "weather_grid_count":p[i]["new_weather_source_grid_count"],
                "weather_grid_HOLD_count":1,
                "original_48_solar_heading_signed_external_not_replicated":True,
                "original_12_ERA5_old_raw_sha_hold":True,
                "GPS_animal_ID_coordinates_dates_weather_grid_key_and_individual_slopes_not_exported":True
            })
        result=agg.aggregate(row)
        self.assertEqual(result["original_source_study_DOI_denominator"],3)
        self.assertEqual(result["nominal_within_month_Tmax_positive_source_count"],0)
        self.assertEqual(result["source_status_counts"],{"SOURCE_HOLD":3})
        with self.assertRaisesRegex(ValueError,"Three"):
            agg.aggregate(row[:2])

if __name__=="__main__":unittest.main()
