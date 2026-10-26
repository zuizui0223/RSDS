"""Fixed source-blob original GPS-only ecology synthesis regression tests."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.summarize_rsds_same_tag_gps_geometry_v2 import (
    ALL_AXES, EXPECTED_BLOBS, GEOMETRY, WEATHER, FIRST,
    ROOT, PRESERVED_VERDICT, run,
)


class SameTagGPSGeometryEcologySynthesisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = run()

    def test_2_studies_6_source_taxa_269_within_original_tag_pairs(self):
        d = self.d
        self.assertEqual(d["study_units"], 2)
        self.assertEqual(d["study_taxon_rows"], 6)
        self.assertEqual(d["summed_source_paired_tag_counts_NOT_independent_biomes"], 269)
        self.assertEqual({x["region"] for x in d["source_taxon_ecological_results"]},
                         {"Utah","Wyoming"})
        self.assertEqual(sorted(x["original_paired_tagged_animals"]
                                for x in d["source_taxon_ecological_results"]),
                         [21,27,33,37,50,101])

    def test_5_temporal_q_intervals_exclude_zero_but_no_joint_circular_delta_does(self):
        d = self.d["each_axis_nominal_tag_bootstrap95_classification_counts"]
        self.assertEqual(d["q"]["POSITIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"], 2)
        self.assertEqual(d["q"]["NEGATIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"], 3)
        self.assertEqual(d["q"]["INCLUDES_ZERO"], 1)
        self.assertEqual(d["sun_cos"]["POSITIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"], 2)
        self.assertEqual(d["sun_cos"]["INCLUDES_ZERO"], 4)
        self.assertEqual(d["joint_circular"], {"INCLUDES_ZERO": 6})
        self.assertEqual(d["sun_sin"]["POSITIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"], 2)
        self.assertEqual(d["sun_sin"]["NEGATIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"], 1)
        self.assertEqual(d["sun_sin"]["INCLUDES_ZERO"], 3)

    def test_bighorn_two_study_q_signed_same_but_contexts_not_equivalent_exposure(self):
        r = self.d["repeated_native_taxa_with_two_original_studies"]
        self.assertEqual(len(r), 2)
        sheep = next(x for x in r if x["taxon"] == "Ovis canadensis")
        q = sheep["q_coverage_rank2_minus_rank1_in_each_original_study"]
        self.assertAlmostEqual(q["Utah"], -.02979419171930651, places=12)
        self.assertAlmostEqual(q["Wyoming"], -.011200213694510438, places=12)
        self.assertEqual(q["Utah_interval_classification"],
                         "NEGATIVE_UNADJUSTED_SOURCE_TAG_INTERVAL")
        self.assertEqual(q["Wyoming_interval_classification"],
                         "NEGATIVE_UNADJUSTED_SOURCE_TAG_INTERVAL")
        self.assertTrue(sheep["no_common_high_low_environmental_exposure_assigned"])

    def test_geographic_marginal_plus_joint_reconstructs_raw_solar_cos_difference(self):
        for item in self.d["source_taxon_ecological_results"]:
            z = item["results"]
            self.assertEqual(set(z), set(ALL_AXES))
            self.assertAlmostEqual(
                z["sun_cos"]["context2_minus_context1_mean"],
                z["geographic_marginal"]["context2_minus_context1_mean"] +
                z["joint_circular"]["context2_minus_context1_mean"],
                places=12
            )

    def test_same_tag_time_direction_rho_not_independence_proof(self):
        for r in self.d["source_taxon_ecological_results"]:
            rho = r["same_tag_rank_rho_delta_q_vs_delta_joint_solar_circular"]
            self.assertIsNotNone(rho)
            self.assertTrue(-1 <= rho <= 1)
        self.assertTrue(any("DOES NOT" in x for x in self.d["interpretation_ceiling"]))

    def test_12_ERA5_HOLDS_and_external_negative_direction_preserved(self):
        d = self.d
        self.assertEqual(d["earlier_weather_version_holds_preserved"], 12)
        self.assertEqual(d["earlier_weather_response_comparison_supported_taxa"], 0)
        self.assertEqual(d["old_external_shortwave_nonreplication_UNCHANGED"],
                         PRESERVED_VERDICT)

    def test_first_original_source_blob_replacement_fails_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for p in EXPECTED_BLOBS:
                dst = root / p
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_bytes((ROOT / p).read_bytes())
            d = json.loads((root / WEATHER).read_text())
            d["native_study_taxon_rows"] = 0
            (root / WEATHER).write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError, "replaced or updated"):
                run(root)

    def test_no_private_tag_coordinates_or_fix_times_in_public_synthesis(self):
        serialized = json.dumps(self.d)
        for key in ('"individual_id"', '"animal_id"', '"location_lat"',
                    '"location_long"', '"timestamp"', '"UTC_observation_time"',
                    '"original_tag_ids"', '"original_grid_latlon"'):
            self.assertNotIn(key, serialized)


if __name__ == "__main__":
    unittest.main()
