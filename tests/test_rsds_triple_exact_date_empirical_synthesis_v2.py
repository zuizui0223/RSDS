"""Frozen source original 3-schedule identical tagged animal/day GPS result tests."""
import json
import tempfile
import unittest
from pathlib import Path

from scripts.summarize_rsds_three_phase_original_ecology_v2 import ROOT,PIN,SOURCE,OLDER,run


class EightOriginalStudyTripleDateEcologyKnownTruth(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.o=run()

    def test_eight_DOIs_twelve_native_taxa_and_all_original_HOLDS(self):
        z=self.o
        self.assertEqual(z["original_field_studies"],8)
        self.assertEqual(z["original_study_taxon_rows"],12)
        self.assertEqual(z["qualifying_same_tag_original_three_phase_source_taxon_rows"],8)
        self.assertEqual(z["source_state_counts"],{
            "TRIPLE_SOURCE_SUPPORT_EVALUATED":7,
            "PREDECLARED_FILE_METADATA_HOLD":1
        })
        self.assertTrue(z["previous_12_ERA5_source_json_version_HOLDS_UNCHANGED"])

    def test_both_original_Wyoming_taxa_same_tag_dates_retain_observed_phase_bias(self):
        z=self.o["original_Wyoming_two_native_taxa_exact_three_phase_recurrence"]
        self.assertEqual(set(z),{"Odocoileus hemionus","Ovis canadensis"})
        expected={
            "Odocoileus hemionus":(113,.005393070134477692),
            "Ovis canadensis":(41,.013338853585317417)
        }
        for taxon,(n,mean) in expected.items():
            r=z[taxon]
            self.assertEqual(r["original_triply_matched_tag_count"],n)
            e,o,d=(r[x]["original_tag_group_mean"] for x in (
                "q_even_minus_original","q_odd_minus_original","q_even_minus_odd"))
            self.assertAlmostEqual(d,mean,places=14)
            self.assertAlmostEqual(d,e-o,places=14)
            self.assertEqual(r["q_even_minus_original"]["nominal_bootstrap_interval_sign"],
                             "NOMINAL_POSITIVE")
            self.assertEqual(r["q_odd_minus_original"]["nominal_bootstrap_interval_sign"],
                             "NOMINAL_NEGATIVE")
            self.assertEqual(r["q_even_minus_odd"]["nominal_bootstrap_interval_sign"],
                             "NOMINAL_POSITIVE")

    def test_positive_triple_date_UTC_phase_bias_only_three_taxa_two_study_DOIs(self):
        r=self.o["qualified_original_source_taxa_with_positive_even_minus_odd_q"]
        self.assertEqual(r["n"],3)
        self.assertEqual({(x["study_doi"],x["taxon"]) for x in r["source_taxa"]},{
            ("10.5441/001/1.711","Cervus elaphus"),
            ("10.5441/001/1.703","Odocoileus hemionus"),
            ("10.5441/001/1.703","Ovis canadensis")
        })
        signs=self.o["source_axis_nominal_interval_counts"]
        self.assertEqual(signs["q_even_minus_odd"],{"NOMINAL_POSITIVE":3,
                 "INCLUDES_ZERO":5,"UNAVAILABLE":3})

    def test_ineligible_original_studies_not_false_evidence_of_no_parity_bias(self):
        d={(x["study_doi"],x["native_taxon"]):x for x in
           self.o["all_eight_source_study_native_taxa_group_effects"]}
        for key in (("10.5441/001/1.671","Meleagris gallopavo"),
                    ("10.5441/001/1.644","Prionailurus viverrinus")):
            z=d[key]["strict_three_phase_same_tag_date_effects"]
            self.assertEqual(z["q_even_minus_odd"]["nominal_bootstrap_interval_sign"],
                             "UNAVAILABLE")
            self.assertIsNone(z["q_even_minus_odd"]["original_tag_group_mean"])
        egret=d[("10.5441/001/1.641","Ardea alba")]
        self.assertEqual(egret["first_source_status"],"PREDECLARED_FILE_METADATA_HOLD")
        self.assertIsNone(egret["strict_three_phase_same_tag_date_effects"])

    def test_old_ERA5_and_48_source_negative_radiation_results_remain_immutable(self):
        z=self.o
        self.assertTrue(z["previous_12_ERA5_source_json_version_HOLDS_UNCHANGED"])
        self.assertEqual(z["previous_48_source_signed_sunlight_external_NONREPLICATED_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_changing_previous_source_original_blob_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            for path in PIN:
                f=p/path
                f.parent.mkdir(parents=True,exist_ok=True)
                f.write_bytes((ROOT/path).read_bytes())
            first=p/SOURCE
            z=json.loads(first.read_text(encoding="utf-8"))
            z["native_source_study_taxon_rows"]=0
            first.write_text(json.dumps(z),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"frozen original GPS source"):
                run(p)

    def test_no_original_source_tags_fix_coords_or_UTC_times_in_audit(self):
        s=json.dumps(self.o)
        for exact in ('"animal_id"','"individual_id"','"tag_id"','"latitude"',
                      '"longitude"','"event_date"','"utc_hour"','"location-lat"',
                      '"location-long"','"source_year"','"private_gps"'):
            self.assertNotIn(exact,s)


if __name__=="__main__":
    unittest.main()
