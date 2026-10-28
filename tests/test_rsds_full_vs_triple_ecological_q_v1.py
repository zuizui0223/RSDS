"""Source-bloc-pinned full versus same-tag same-date time ecology sign audit."""
import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_rsds_full_vs_matched_temporal_ecology_v1 import (
    ROOT,SHA,OLD,MATCHED,run,compare_sign
)


class FullVersusMatchedOriginalTimeNicheEcology(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=run()

    def test_eight_field_studies_twelve_taxon_rows_plus_strict_missing(self):
        a=self.report
        self.assertEqual(a["original_field_study_doi_count"],8)
        self.assertEqual(a["source_taxon_rows"],12)
        self.assertEqual(len(a["full_source_and_strict_same_original_tag_date_results"]),12)
        by={r["native_taxon"]:r for r in a["full_source_and_strict_same_original_tag_date_results"]}
        self.assertEqual(by["Meleagris gallopavo"]["original_identical_triple_matched_q_n_tagged_animals"],0)
        self.assertEqual(by["Prionailurus viverrinus"]["full_vs_matched_source_q_inference_class"],
                         "UNAVAILABLE_TRIPLE_SAME_TAG_DATE")

    def test_five_source_q_contrasts_survive_full_and_all_three_schedules(self):
        a=self.report
        q=a["full_vs_triple_q_support_sign_count"]
        self.assertEqual(q["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_POSITIVE"],3)
        self.assertEqual(q["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_NEGATIVE"],2)
        self.assertEqual(a["n_distinct_original_field_studies_with_q_sign_supported_in_full_and_all_three_schedules"],5)
        keep={(x["original_study_DOI"],x["native_taxon"]) for x in a[
            "full_source_and_strict_same_original_tag_date_results"]
            if x["full_vs_matched_source_q_inference_class"].startswith("SAME_SIGN_SUPPORTED")}
        self.assertEqual(keep,{
            ("10.5441/001/1.711","Ovis canadensis"),
            ("10.5441/001/1.703","Ovis canadensis"),
            ("10.5441/001/1.323","Lynx rufus"),
            ("10.5441/001/1.304","Vulpes vulpes"),
            ("10.5441/001/1.740","Phascolarctos cinereus")
        })

    def test_only_two_original_studies_retain_same_q_daylength_sign_full_and_three(self):
        d=self.report
        x=d["full_vs_triple_daylength_slope_sign_count"]
        self.assertEqual(x["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_POSITIVE"],2)
        self.assertEqual(d["n_distinct_original_field_studies_with_daylength_sign_supported_in_full_and_all_three_schedules"],2)
        taxa={(z["original_study_DOI"],z["native_taxon"]) for z in
              d["full_source_and_strict_same_original_tag_date_results"]
              if z["full_vs_matched_source_daylength_inference_class"].startswith("SAME_SIGN_SUPPORTED")}
        self.assertEqual(taxa,{
            ("10.5441/001/1.703","Odocoileus hemionus"),
            ("10.5441/001/1.323","Lynx rufus")
        })

    def test_original_full_vs_strict_animal_tag_source_ratio_is_reported_not_misread_as_ecology(self):
        rows={(x["original_study_DOI"],x["native_taxon"]):x for x in self.report[
            "full_source_and_strict_same_original_tag_date_results"]}
        elk=rows[("10.5441/001/1.711","Cervus elaphus")]
        self.assertEqual(elk["full_original_source_q_n_tagged_animals"],337)
        self.assertEqual(elk["original_identical_triple_matched_q_n_tagged_animals"],6)
        self.assertAlmostEqual(elk["original_same_tag_day_retention_fraction_of_full_q_tags"],6/337)
        mule=rows[("10.5441/001/1.703","Odocoileus hemionus")]
        self.assertEqual(mule["full_original_source_q_n_tagged_animals"],119)
        self.assertEqual(mule["original_identical_triple_matched_q_n_tagged_animals"],113)

    def test_strict_stability_with_full_source_ci_zero_is_not_full_original_study_replication(self):
        rows={(x["original_study_DOI"],x["native_taxon"]):x for x in self.report[
            "full_source_and_strict_same_original_tag_date_results"]}
        koala=rows[("10.5441/001/1.740","Phascolarctos cinereus")]
        self.assertEqual(koala["same_tag_three_schedule_daylength_sign"],"THREE_SCHEDULES_NOMINAL_POSITIVE")
        self.assertEqual(koala["full_original_source_q_daylength_slope"]["status"],"INCLUDES_ZERO")
        self.assertEqual(koala["full_vs_matched_source_daylength_inference_class"],
                         "MATCHED_ONLY_SIGN_SUPPORTED_NOT_FULL_NATIVE_SOURCE_POSITIVE")

    def test_historic_negative_external_and_12_weather_source_HOLDs_untouched(self):
        self.assertTrue(self.report["old_12_ERA5_serialized_raw_SHA_HOLDS_UNCHANGED"])
        self.assertEqual(self.report["old_48_source_external_signed_radiation_nonreplication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_changed_original_source_first_result_sha_prevents_new_ecological_claim(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in SHA:
                p=root/name
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_bytes((ROOT/name).read_bytes())
            raw=json.loads((root/MATCHED).read_text())
            raw["original_study_taxon_rows"]=7
            (root/MATCHED).write_text(json.dumps(raw))
            with self.assertRaisesRegex(ValueError,"Original source first scientific result changed"):
                run(root)

    def test_old_source_signed_vs_new_same_sample_labels(self):
        self.assertEqual(compare_sign(
            "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL","THREE_SCHEDULES_NOMINAL_POSITIVE"),
            "SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_POSITIVE")
        self.assertEqual(compare_sign(
            "INCLUDES_ZERO","THREE_SCHEDULES_NOMINAL_NEGATIVE"),
            "MATCHED_ONLY_SIGN_SUPPORTED_NOT_FULL_NATIVE_SOURCE_NEGATIVE")
        self.assertEqual(compare_sign(
            "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL",
            "UNRESOLVED_AT_LEAST_ONE_ZERO_SPANNING_INTERVAL"),
            "UNRESOLVED_IN_ONE_OR_MORE_MATCHED_GPS_SCHEDULES")


if __name__=="__main__":
    unittest.main()
