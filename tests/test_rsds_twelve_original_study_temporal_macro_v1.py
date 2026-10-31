"""Source-locked 12 original Movebank-study 22 native q/daylength macro ecology tests."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import scripts.build_rsds_twelve_original_study_temporal_macro_v1 as m


class OriginalTwelveStudyTemporalMacro(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=m.main_analysis()

    def test_twelve_original_studies_twenty_two_native_rows_ten_new_taxa(self):
        o=self.data
        self.assertEqual(o["all_original_studies"],12)
        self.assertEqual(o["all_source_native_taxon_rows"],22)
        self.assertEqual(o["distinct_native_species_labels"],20)
        self.assertEqual(o["new_original_studies"],4)
        self.assertEqual(len(o["all_original_study_native_taxon_rows_including_HOLD"]),22)
        self.assertEqual(len({r["original_study_doi"] for r in o["all_original_study_native_taxon_rows_including_HOLD"]}),12)

    def test_9_source_and_15_native_q_and_photoperiod_support(self):
        d=self.data
        for k in ("temporal_q_original_study_equal_inference",
                  "astronomical_daylength_original_study_equal_inference"):
            s=d[k]
            self.assertEqual(s["original_field_study_DOIs_including_HOLD"],12)
            self.assertEqual(s["original_studies_with_estimable_native_source_mean"],9)
            self.assertEqual(s["native_source_taxon_rows_with_estimate"],15)
            self.assertEqual(len(s["all_original_DOI_source_means_including_unavailable"]),12)

    def test_signed_observational_nonuniversality_with_all_HOLDs(self):
        d=self.data
        self.assertEqual(d["native_source_q_nominal_animal_CI_sign_counts"],
                         {"POSITIVE":4,"NEGATIVE":11,"UNAVAILABLE":7})
        self.assertEqual(d["native_source_daylength_nominal_animal_CI_sign_counts"],
                         {"POSITIVE":5,"NEGATIVE":2,"INCLUDES_ZERO":8,"UNAVAILABLE":7})
        source=d["temporal_q_original_study_equal_inference"]
        self.assertEqual(source["original_DOI_source_mean_point_signs"],
                         {"POSITIVE":4,"NEGATIVE":5,"EXACT_ZERO":0})

    def test_source_equal_and_row_equal_macro_means_diverge_by_units(self):
        d=self.data
        recent=d["temporal_q_original_study_equal_inference"]
        prior=d["first_eight_study_original_q_DOI_equal_reference"]
        self.assertAlmostEqual(prior["source_DOI_equal_mean"],0.03464040880669044,places=12)
        self.assertAlmostEqual(recent["source_DOI_equal_mean"],0.010262345757767848,places=12)
        self.assertAlmostEqual(recent["native_taxon_row_equal_mean_OVERWEIGHTS_MULTI_TAXON_STUDIES"],
                               -0.022509718658262047,places=12)
        self.assertLess(recent["native_taxon_row_equal_mean_OVERWEIGHTS_MULTI_TAXON_STUDIES"],0)
        self.assertGreater(recent["source_DOI_equal_mean"],0)
        self.assertEqual(len(recent["source_DOI_bootstrap95_DESCRIPTIVE_NONPROBABILITY"]),2)
        self.assertTrue(recent["not_global_CI"])

    def test_new_goose_120_shared_days_are_not_the_full_source_annual_effect(self):
        d=self.data
        self.assertEqual(d["new_source_native_goose_original_shared_calendar_days"],120)
        self.assertEqual(d["new_source_native_goose_original_year_month_blocks"],5)
        self.assertEqual(len(d["one_source_coobserved_three_goose_native_q_and_reference_daylength"]),3)
        self.assertEqual(len(d["one_source_three_goose_pair_q_and_month_reference_slopes"]),3)
        branta=next(r for r in d["all_original_study_native_taxon_rows_including_HOLD"]
                    if r["native_taxon"]=="Branta leucopsis")
        common=next(r for r in d["one_source_coobserved_three_goose_native_q_and_reference_daylength"]
                    if r["native_taxon"]=="Branta leucopsis")
        self.assertAlmostEqual(branta["q"]["mean"],-0.1654734745042472,places=12)
        self.assertAlmostEqual(common["shared_120_date_q"],-0.1840335077564767,places=12)
        self.assertNotEqual(branta["daylength"]["mean"],
                            common["shared_120_date_reference_daylength_slope_per_hour"])

    def test_original_schema_and_cadence_unavailable_are_not_zeros(self):
        d=self.data
        rows=d["all_original_study_native_taxon_rows_including_HOLD"]
        for taxon in ("Chelonoidis hoodensis","Chelonoidis porteri",
                      "Chelonoidis donfaustoi","Torgos tracheliotus",
                      "Gyps africanus","Ardea alba","Anas acuta"):
            native=next(x for x in rows if x["native_taxon"]==taxon)
            self.assertEqual(native["q"]["status"],"UNAVAILABLE")
            self.assertIsNone(native["q"]["mean"])
            self.assertEqual(native["daylength"]["status"],"UNAVAILABLE")
        self.assertEqual(
            len([x for x in rows if x["original_source_stratum"]=="FOUR_NEW_SOURCE_FROZEN_STUDIES"]),10
        )

    def test_deterministic_bootstrap_and_DOI_not_taxon_unit(self):
        first=m.main_analysis()
        second=m.main_analysis()
        self.assertEqual(first,second)
        studies=first["temporal_q_original_study_equal_inference"][
            "all_original_DOI_source_means_including_unavailable"]
        ut=next(x for x in studies if x["source_study_doi"]=="10.5441/001/1.711")
        goose=next(x for x in studies if x["source_study_doi"]=="10.5441/001/1.fk899541")
        self.assertEqual(ut["total_native_taxon_rows"],4)
        self.assertEqual(goose["total_native_taxon_rows"],3)
        self.assertEqual(ut["estimable_native_taxa"],4)
        self.assertEqual(goose["estimable_native_taxa"],3)

    def test_old_radiation_negative_and_weather_sha_HOLD_immutable(self):
        d=self.data
        self.assertEqual(d["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"],
                         m.NEGATIVE)
        self.assertTrue(d["original_12_ERA5_raw_JSON_source_HOLDS_UNCHANGED"])

    def test_changed_prior_first_source_blob_fails_closed(self):
        spec=json.loads(m.CONTRACT.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for f in spec["source_blobs"]:
                path=f["path"]
                target=root/path
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes((m.ROOT/path).read_bytes())
            bad=root/spec["source_blobs"][1]["path"]
            obj=json.loads(bad.read_text(encoding="utf-8"))
            obj["fixed_original_published_studies"]=2
            bad.write_text(json.dumps(obj),encoding="utf-8")
            with patch.object(m,"ROOT",root):
                with self.assertRaisesRegex(ValueError,"Frozen original group scientific result changed"):
                    m.main_analysis()


if __name__=="__main__":
    unittest.main()
