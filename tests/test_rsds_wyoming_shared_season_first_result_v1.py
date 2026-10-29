"""Evidence lock and known values for one-study coobserved Wyoming photoperiod contrast."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import scripts.validate_rsds_wyoming_shared_season_ecology_v1 as model


class WyomingCoObservedCalendarEcologyValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=model.validate()

    def test_all_original_source_provider_and_two_taxa_denominator(self):
        a=self.d
        self.assertEqual(a["original_study_DOI"],"10.5441/001/1.703")
        self.assertEqual(a["original_field_study_replications"],1)
        self.assertEqual(set(a["original_native_taxa"]),
                         {"Odocoileus hemionus","Ovis canadensis"})
        self.assertEqual(a["source_coobserved_calendar_dates"],193)
        self.assertEqual(a["distinct_original_calendar_years"],2)
        self.assertEqual(a["distinct_year_month_blocks_for_date_resampling"],8)
        self.assertEqual(a["astronomical_daylength_common_reference_latitude_degrees"],43)

    def test_mule_positive_bighorn_negative_each_same_study_reference_daylight(self):
        e=self.d["original_taxon_source_tagged_95pct_results"]
        self.assertEqual(e["Odocoileus hemionus"]["source_original_tag_n"],108)
        self.assertEqual(e["Ovis canadensis"]["source_original_tag_n"],41)
        for phase in ("original","even","odd"):
            positive=e["Odocoileus hemionus"]["common_date_reference_photoperiod_coefficients"][phase]
            negative=e["Ovis canadensis"]["common_date_reference_photoperiod_coefficients"][phase]
            self.assertEqual(positive["nominal_sign"],"POSITIVE")
            self.assertEqual(negative["nominal_sign"],"NEGATIVE")
            self.assertGreater(positive["animal_bootstrap95"][0],0)
            self.assertLess(negative["animal_bootstrap95"][1],0)

    def test_between_native_tag_diff_and_equal_calendar_block_diff_both_positive(self):
        a=self.d
        b=a["within_study_mule_minus_bighorn_tag_bootstrap_slopes"]
        c=a["paired_same_calendar_date_year_month_bootstrap_slopes"]
        reference={
           "original":(.012347378623912636,.011673653876371936),
           "even":(.01280098790008331,.013358456372272055),
           "odd":(.015842745850824742,.011558629015356208)
        }
        for phase,(tag_delta,day_delta) in reference.items():
            self.assertAlmostEqual(b[phase]["mule_minus_bighorn_slope"],tag_delta,places=12)
            self.assertAlmostEqual(c[phase]["mule_q_minus_bighorn_q_vs_shared_reference_daylength"],
                                   day_delta,places=12)
            self.assertGreater(b[phase]["descriptive_original_tag_bootstrap95"][0],0)
            self.assertGreater(c[phase]["year_month_block_descriptive95"][0],0)

    def test_original_raw_era5_HOLDs_and_external_negative_not_changed(self):
        r=self.d
        self.assertTrue(r["previous_twelve_ERA5_raw_source_version_HOLDS_UNCHANGED"])
        self.assertEqual(r["previous_forty_eight_source_negative_direction_verdict_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertFalse(r["previous_forty_eight_source_negative_radiation_direction_replicated"])

    def test_source_publication_never_leaks_original_tag_or_calendar_year_id(self):
        j=json.dumps(self.d)
        for exact in ('"individual_id"','"animal_id"','"tag_id"','"source_year"',
                      '"calendar_date"','"latitude"','"longitude"','"UTC_timestamp"',
                      '"GPS_event"','"private_fix"'):
            self.assertNotIn(exact,j)

    def test_prior_source_sha_changing_even_one_byte_blocks_new_interpretation(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for path in (model.SOURCE,model.PREVIOUS,model.CONTRACT):
                t=root/path
                t.parent.mkdir(parents=True,exist_ok=True)
                t.write_bytes((model.ROOT/path).read_bytes())
            src=root/model.SOURCE
            payload=json.loads(src.read_text())
            payload["shared_support"]["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]=194
            src.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,"frozen source"):
                model.validate(root)


if __name__=="__main__":
    unittest.main()
