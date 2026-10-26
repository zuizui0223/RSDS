"""Eight original source-study equal-weight macro ecology known truth tests."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.source_equal_rsds_eight_study_temporal_solar_v1 import (
    ROOT, INPUT, CONTRACT, SEED, equal_source_bootstrap, run
)


class EightOriginalSourceStudyEqualMacroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=run()

    def test_original_study_unit_prevents_four_utah_taxa_becoming_four_ecosystems(self):
        r=self.result
        self.assertEqual(r["original_source_studies"],8)
        self.assertEqual(r["original_source_taxon_rows"],12)
        frame=r["fixed_original_source_sample_frame"]["study_native_taxon_counts"]
        self.assertEqual(sorted(frame.values()),[1,1,1,1,1,1,2,4])
        self.assertEqual(frame["10.5441/001/1.711"],4)
        self.assertEqual(frame["10.5441/001/1.703"],2)

    def test_seven_studies_with_temporal_q_and_one_egret_metadata_hold(self):
        r=self.result
        self.assertEqual(r["fixed_original_source_sample_frame"]["source_with_only_metadata_HOLD"],
                         ["10.5441/001/1.641"])
        q=r["source_equal_axes"]["night_day_q"]
        self.assertEqual(q["original_studies_with_supported_axis"],7)
        self.assertEqual(q["original_source_taxon_rows_with_supported_axis"],11)
        self.assertEqual(q["supported_original_source_signs"],{
            "positive":4,"negative":3,"exact_zero":0
        })
        self.assertAlmostEqual(q["equal_source_mean"],0.03464040880669044,places=12)
        self.assertAlmostEqual(q["equal_native_taxon_row_mean_NOT_INDEPENDENT_STUDIES"],
                               0.004796829023221451,places=12)

    def test_source_reweighting_is_not_a_new_world_population_parameter(self):
        q=self.result["source_equal_axes"]["night_day_q"]
        self.assertAlmostEqual(q["equal_source_minus_equal_taxon_row"],
                               0.02984357978346899,places=12)
        self.assertEqual(len(q["equal_source_95percentile_study_resampling_CONDITIONAL_NONPROBABILITY"]),2)
        low,high=q["leave_one_original_study_out_equal_source_mean_range"]
        self.assertTrue(low <= q["equal_source_mean"] <= high)
        self.assertTrue(any("purposively" in x.lower() for x in self.result["scientific_ceiling"]))

    def test_photoperiod_and_sun_source_support_remains_not_imputed(self):
        a=self.result["source_equal_axes"]
        self.assertEqual(a["q_daylength_per_hour"]["original_studies_with_supported_axis"],7)
        self.assertEqual(a["raw_solar_relative_cos"]["original_studies_with_supported_axis"],6)
        self.assertEqual(a["joint_circular_covariance"]["original_studies_with_supported_axis"],6)
        for axis in a:
            first=a[axis]["first_original_study_QA"]
            self.assertEqual(len(first),8)
            egret=next(x for x in first if x["study_doi"]=="10.5441/001/1.641")
            self.assertIsNone(egret["study_equal_taxon_mean"])
            self.assertEqual(egret["study_mean_sign"],"UNAVAILABLE")

    def test_old_joint_circular_point_only_does_not_gain_fictional_ci(self):
        r=self.result["source_equal_axes"]["joint_circular_covariance"]
        self.assertEqual(r["uncertainty_ceiling"],
                         "OLD_SIX_JOINT_ESTIMATES_ARE_POINT_ONLY_NO_TAG_INTERVAL")
        self.assertFalse(r["source_resampling_interval_available"])
        self.assertIsNone(r["equal_source_95percentile_study_resampling_CONDITIONAL_NONPROBABILITY"])

    def test_bootstrap_is_deterministic_by_study_id_not_original_gps_events(self):
        self.assertEqual(equal_source_bootstrap([-.2,.1,.15,.5],seed=SEED,draws=250),
                         equal_source_bootstrap([-.2,.1,.15,.5],seed=SEED,draws=250))

    def test_altered_pinned_original_source_result_blocks_new_claims(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for p in (INPUT,CONTRACT):
                target=root/p.relative_to(ROOT)
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes(p.read_bytes())
            old=json.loads((root/INPUT.relative_to(ROOT)).read_text(encoding="utf-8"))
            old["source_taxon_rows"]=99
            (root/INPUT.relative_to(ROOT)).write_text(json.dumps(old),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"changed or contract"):
                run(root)

    def test_fixed_prior_negative_radiation_and_ERA5_versions_unchanged(self):
        self.assertEqual(self.result["old_external_shortwave_direction_nonreplication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertTrue(self.result["old_12_ERA5_provider_source_SHA_HOLDS_UNCHANGED"])


if __name__=="__main__":
    unittest.main()
