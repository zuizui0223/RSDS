"""Three original source weather data study unit and known-truth synthetic tests."""
import unittest
from datetime import date,timedelta
from unittest.mock import patch
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as m
from scripts import aggregate_rsds_three_mule_original_daily_hydrothermal_v1 as agg
from scripts.movebank_wyoming_daily_climate_20261023 import adjusted_daily_effects


def faux_endpoints(mean=.004,n=10):
    return {name:{
        "admitted_animals":n,"admitted_animal_years":n,
        "mean":mean,"animal_bootstrap_95":[mean-.001,mean+.001],
        "status":"descriptive_within_original_animal_year"}
        for name in agg.ENDPOINTS}


def sample_result(source,index,complete=True):
    status=("ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED" if complete
            else "SOURCE_HOLD")
    d={**faux_endpoints(),
       "animal_years_with_at_least_ten_climate_complete_source_days_and_solar_span":10,
       "weather_complete_original_animal_days":100,
       "all_weather_grid_inputs_retrieved":True}
    return {
        "version":m.SOURCE_VERSION,
        "source_index":index,"source_key":source["source_key"],
        "original_study_DOI":source["doi"],"original_study_id":str(source["study_id"]),
        "original_native_taxon":source["taxon"],
        "registered_original_native_taxa":source["registered_taxa"],
        "original_source_file_MD5_verified":complete,
        "original_source_native_q_and_daylength_replayed_before_weather":complete,
        "original_GPS_source_event_rows_private_only":1000 if complete else 0,
        "source_status":status,
        "source_hold_stage":"PREVIOUS_PROVIDER_ORIGINAL_UUID_SIZE_MD5" if not complete else None,
        "source_hold_type":"OSError" if not complete else None,
        "source_weather_context_version":"NEW_20261116_ERA5_DAILY_NOT_OLD_RAW_SHA_REPLAY",
        "source_ERA5_cell_count":2 if complete else 0,
        "source_ERA5_grid_HOLD_count":0,
        "new_weather_physical_SHA256_not_coordinates":"9"*64 if complete else None,
        "source_weather_grid_results_differ_from_historical_not_compared":True,
        "original_q_eligible_tag_years":10 if complete else 0,
        "all_source_native_weather_adjusted_effects":d if complete else None,
        "original_48_sunlight_heading_external_NOT_REPLICATED":True,
        "original_12_ERA5_raw_SHA_HOLD_unchanged":True,
        "no_original_tag_ID_exact_dates_coordinates_GPS_weather_events_exported":True
    }


class OriginalMuleWeatherTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c,cls.source=m.frozen()

    def test_exact_all_three_mule_source_IDs(self):
        self.assertEqual(len(self.source),3)
        self.assertEqual([s["source_key"] for s in self.source],
             ["mule_old_utah","mule_old_wyoming","mule_new_sublette"])
        self.assertEqual(len({s["study_id"] for s in self.source}),3)
        self.assertTrue(all(s["taxon"]=="Odocoileus hemionus" for s in self.source))
        self.assertTrue(self.c["no_original_trajectory_tag_coordinates_dates_or_individual_slopes_exported"])

    def test_synthetic_weather_grid_original_day_counts(self):
        t=date(2020,2,1)
        days=[(t+timedelta(days=i),12.0,.04,(40,-100),(2020,2)) for i in range(15)]
        # Source admission is performed only by the frozen helper in production;
        # test received original selected tag-year tuples at the private boundary.
        dummy=[("tagA",2020,days)]
        with patch.object(m,"source_admitted_daily",return_value=dummy):
            a,c=m.source_days({"tagA":{}})
        self.assertEqual(a,dummy)
        self.assertEqual(len(c),1)
        self.assertEqual(len(c[(40,-100)]),15)

    def test_weather_physical_fingerprint_ignores_dict_insertion_order(self):
        a={(40,-100):{date(2020,1,1):(1,2,3,4),date(2020,1,2):None},
           (41,-99):{date(2020,1,1):(0,2,0,2)}}
        b={(41,-99):a[(41,-99)],(40,-100):dict(reversed(list(a[(40,-100)].items())))}
        self.assertEqual(m.physical_fingerprint(a),m.physical_fingerprint(b))
        self.assertEqual(len(m.physical_fingerprint(a)),64)

    def test_source_provider_error_is_HOLD_not_a_weather_or_biology_effect(self):
        with patch.object(m,"original_files",side_effect=OSError("PRIVATE_LOCATION_OR_KEY")) as p:
            r=m.source(0)
            p.assert_called_once()
        self.assertEqual(r["source_status"],"SOURCE_HOLD")
        self.assertFalse(r["original_source_native_q_and_daylength_replayed_before_weather"])
        self.assertEqual(r["source_hold_type"],"OSError")
        self.assertIsNone(r["all_source_native_weather_adjusted_effects"])
        self.assertNotIn("PRIVATE_LOCATION_OR_KEY",str(r))

    def test_all_studies_kept_missingness_is_explicit(self):
        rows=[sample_result(s,i) for i,s in enumerate(self.source)]
        result=agg.aggregate(rows)
        self.assertEqual(result["source_study_DOI_denominator"],3)
        self.assertEqual(result["full_weather_coverage_source_studies"],3)
        self.assertEqual(len(result["all_three_original_source_rows_including_HOLD"]),3)
        self.assertEqual(result["primary_daylength_joint_adjusted_nominal_CI_signs_by_fully_covered_source"],{
            "mule_old_utah":"POSITIVE","mule_old_wyoming":"POSITIVE",
            "mule_new_sublette":"POSITIVE"})
        rows[2]=sample_result(self.source[2],2,complete=False)
        result=agg.aggregate(rows)
        self.assertEqual(result["source_study_DOI_denominator"],3)
        self.assertEqual(result["full_weather_coverage_source_studies"],2)
        self.assertEqual(result["original_source_status_counts"]["SOURCE_HOLD"],1)
        with self.assertRaisesRegex(ValueError,"3"):
            agg.aggregate(rows[:2])

    def test_source_bootstrap_interval_unavailable_not_zero(self):
        self.assertEqual(agg.categorize({"mean":None}),"UNAVAILABLE")
        self.assertEqual(agg.categorize(faux_endpoints(mean=-.01)["daylight_given_tmax_and_water"]),"NEGATIVE")
        self.assertEqual(agg.categorize(faux_endpoints(mean=.01)["daylight_given_tmax_and_water"]),"POSITIVE")


if __name__=="__main__":
    unittest.main()
