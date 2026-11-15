"""Known-truth source-group aggregation and fail-closed first biology."""
from __future__ import annotations

import copy
import unittest

from scripts import aggregate_rsds_35_original_gps_observation_funnel_v1 as a

FIRST=a.pinned(*a.EARLIER)["combined_source_study_summary_rows"]

def synthetic_funnel():
    result=[]
    for i,record in enumerate(FIRST):
        wave,index=("wave1",i) if i<16 else ("wave2",i-16)
        prior=record["original_source_status"]
        is_human=(prior=="HUMAN_LOCATION_PRIVACY_HOLD")
        held=prior in (
            "PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD","SOURCE_FILE_METADATA_HOLD",
            "HUMAN_LOCATION_PRIVACY_HOLD","ORIGINAL_SOURCE_GPS_HOLD","SOURCE_GPS_HOLD")
        x={
          "version":"RSDS_35_ORIGINAL_GPS_SOURCE_OBSERVATION_FUNNEL_V1",
          "wave":wave,"index":index,
          "source_original_DOI":record["source_DOI_labels"][0],
          "source_repository_study_id":record["source_original_repository_study_id"],
          "original_registered_native_taxa":["Synthetic_species_"+str(j)
                                             for j in range(record["registered_native_count"])],
          "prior_first_source_q_status":prior,
          "original_first_q_study_complete":record["complete_source_native_q"],
          "original_first_daylength_study_complete":record["complete_source_daylength"],
          "original_MD5_replayed":not held,
          "source_status":("HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS" if is_human
                           else "ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED" if held
                           else "REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL"),
          "original_native_funnel":[],
          "source_UTC_observation_cadence_replay_status":("NO_HUMAN_GPS_ACCESS" if is_human else
              "NO_ADMISSIBLE_ARCHIVED_NATIVE_SOURCE" if held else "REPLAY_VERIFIED"),
          "original_raw_original_event_count":(
             None if held else record["original_source_event_row_count"]),
          "source_hold_reason":record["source_hold_error_type"] if held else None,
          "original_48_signed_sunlight_external_NONREPLICATION_UNCHANGED":a.NEG,
          "original_12_raw_ERA5_source_HOLD_UNCHANGED":True,
          "no_individual_ID_GPS_location_date_timestamp_or_q_values_exported":True,
        }
        if not held:
            source_native=record["registered_native_count"]
            for n in range(source_native):
                elig=n<record["native_q_supported"]
                elig_day=n<record["native_daylength_supported"]
                x["original_native_funnel"].append({
                  "native_taxon":"Synthetic_species_"+str(n),
                  "source_original_support_counts":{
                    "source_original_tagged_animals_in_hourly_parsed_track":10,
                    "original_selected_UTC_hour_fix_count":2400,
                    "UTC_even_hour_selected_fix_count":1200,
                    "UTC_odd_hour_selected_fix_count":1200,
                    "original_adjacent_source_fix_intervals":2300,
                    "accepted_original_0p5_to_2h_time_intervals":1800,
                    "accepted_same_solar_date_and_daypart_intervals":1600,
                    "accepted_finite_speed_max60kmh_intervals":1550,
                    "accepted_native_tag_solar_dates_with_ge3_each_daypart":150,
                    "accepted_native_tag_years_with_ge10_daypart_dates":8 if elig else 0,
                    "accepted_native_tags_with_ge10_daypart_dates":8 if elig else 0,
                    "accepted_native_tag_years_with_ge10_dates_and_daylength_span_ge0p5h":7 if elig_day else 0,
                    "accepted_native_tags_with_within_year_daylength_slope":7 if elig_day else 0
                  },
                  "original_observation_gate_ratios":{
                    "original_time_0p5_to_2h_over_adjacent":1800/2300,
                    "same_solar_date_daypart_over_time_eligible":1600/1800,
                    "plausible_speed_over_solar_eligible":1550/1600,
                    "total_accepted_over_adjacent":1550/2300,
                    "UTC_even_hour_fraction":.5
                  },
                  "replayed_archived_native_q_support":{
                    "n_tag_q":8 if elig else 0,
                    "n_tag_daylength":7 if elig_day else 0,
                    "q_estimable":elig,"daylength_estimable":elig_day},
                  "no_animal_ID_date_location_GPS_row_or_movement_values_exported":True
                })
        result.append(x)
    return result


class SourceEqualObservationFunnel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=synthetic_funnel()
        cls.result=a.summarize(cls.rows)

    def test_frozen35_sources_and_original_HOLDs(self):
        r=self.result
        self.assertEqual(r["source_original_registered_native_rows"],50)
        self.assertEqual(len(r["individual_original_source_diagnostic_rows"]),35)
        self.assertEqual(sum(r["original_source_diagnostic_status_counts"].values()),35)
        self.assertEqual(r["original_source_diagnostic_status_counts"],{
            "REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL":25,
            "ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED":9,
            "HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS":1,
        })
        self.assertEqual(r["original_human_GPS_source_opened"],False)
        self.assertEqual(r["original_preaccess_file_or_privacy_HOLDS"],7)
        self.assertEqual(r["original_sourcewide_GPS_HOLDS"],3)
        self.assertEqual(r["first_source_study_q_complete"],10)

    def test_source_equal_ratios_do_not_use_event_counts_as_independent_studies(self):
        r=self.result
        complete=r["source_equal_observation_gate_contrasts_by_archived_q_availability"][
            "complete_native_q"]
        incomplete=r["source_equal_observation_gate_contrasts_by_archived_q_availability"][
            "not_complete_native_q"]
        self.assertEqual(complete["parsed_source_study_count"],10)
        self.assertEqual(incomplete["parsed_source_study_count"],15)
        self.assertAlmostEqual(complete["source_equal_median_original_time_0p5_to_2h_over_adjacent"],
                               1800/2300)
        self.assertAlmostEqual(incomplete["source_equal_median_UTC_even_hour_fraction"],.5)
        # Nested GPS intervals descriptive only; cannot set source N to segment count.
        self.assertGreater(r["original_per_fix_gate_census_NESTED_but_NOT_independent_samples"][
            "original_adjacent_source_fix_intervals"],25000)
        self.assertEqual(r["source_original_registered_native_rows"],50)

    def test_no_replacing_original_HOLD_by_zero(self):
        s=copy.deepcopy(self.rows)
        s[0]["source_status"]="REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL"
        with self.assertRaisesRegex(ValueError,"promoted to HOLD or vice versa"):
            a.summarize(s)
        s=copy.deepcopy(self.rows)
        s[26]["source_status"]="REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL" # human Wave2 index10
        with self.assertRaisesRegex(ValueError,"promoted to HOLD or vice versa"):
            a.summarize(s)

    def test_original_tag_group_not_dropped_on_multi_native(self):
        s=copy.deepcopy(self.rows)
        i=next(j for j,x in enumerate(FIRST)
               if x["registered_native_count"]>1 and x["source_stratum"]=="multi_native"
               and s[j]["source_status"]=="REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL")
        s[i]["original_native_funnel"].pop()
        with self.assertRaisesRegex(ValueError,"source-native measured support not replayed"):
            a.summarize(s)

    def test_observation_time_and_solar_gate_are_distinct(self):
        s=copy.deepcopy(self.rows)
        i=2 # first parsed Wave1 source
        item=s[i]["original_native_funnel"][0]["source_original_support_counts"]
        item["accepted_original_0p5_to_2h_time_intervals"]=1500
        with self.assertRaisesRegex(ValueError,"observation gate not nested"):
            a.summarize(s)
        s=copy.deepcopy(self.rows)
        item=s[i]["original_native_funnel"][0]["original_observation_gate_ratios"]
        item["UTC_even_hour_fraction"]=1.2
        with self.assertRaisesRegex(ValueError,"cadence gate proportion"):
            a.summarize(s)

    def test_old_negative_external_and_prior_source_blob_immutable(self):
        self.assertEqual(self.result["original_48_solar_negative_external_NONREPLICATION_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertTrue(self.result["old_12_ERA5_raw_source_versions_HOLD_UNCHANGED"])


if __name__=="__main__":
    unittest.main()
