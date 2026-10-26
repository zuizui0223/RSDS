"""Known-truth tests for original tagged individual paired gridyear sensitivity."""
from __future__ import annotations

import json
import unittest
from copy import deepcopy

from scripts.run_rsds_same_tag_context_pair_v1 import (
    BASELINE, BOOT, CONTRACT, check_context_before_paired, first_sources,
    paired_context_effects, paired_percentiles, ranked_spearman,
    verify_first_context_effects,
)


class OriginalTaggedAnimalContextPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract, cls.first = first_sources()

    def test_frozen_first_study_frame_retains_two_studies_six_taxa(self):
        self.assertEqual(self.contract["version"], "RSDS_SAME_TAG_CONTEXT_PAIR_V1")
        self.assertEqual(set(self.first), {"Utah", "Wyoming"})
        n = sum(len(v["taxon_results"]) for v in self.first.values())
        self.assertEqual(n, 6)
        for source in self.first.values():
            self.assertTrue(source["source_original_replayed_before_NEW_weather"])
            for taxon in source["taxon_results"]:
                self.assertEqual(len(taxon["selected_contexts"]), 2)

    def test_paired_original_animal_not_equal_context_mean_difference(self):
        a = {
            "tagA": {"sun_cos": -0.2, "sun_sin": 0.1, "q": -0.01},
            "tagB": {"sun_cos": 0.1, "sun_sin": -0.2, "q": 0.01},
            "tagC": {"sun_cos": 0.2, "sun_sin": 0.3, "q": -0.02},
            "tagD": {"sun_cos": -0.1, "sun_sin": 0.4, "q": 0.02},
            "tagE": {"sun_cos": 0.0, "sun_sin": -0.4, "q": 0.03},
            "only_initial": {"sun_cos": 99, "sun_sin": 99, "q": 99}
        }
        b = {tag: {"sun_cos": z["sun_cos"] + .05,
                   "sun_sin": z["sun_sin"] - .01,
                   "q": z["q"] + .004} for tag, z in a.items()
             if tag != "only_initial"}
        b["only_later"] = {"sun_cos": -999, "sun_sin": -999, "q": -999}
        r = paired_context_effects(a, b)
        self.assertEqual(r["identical_tag_ids_in_both_original_contexts"], 5)
        self.assertEqual(r["all_three_axis_paired_original_tags"], 5)
        self.assertEqual(r["per_axis"]["q"]["status"], "PAIRED_ESTIMABLE")
        self.assertAlmostEqual(r["per_axis"]["sun_cos"]["rank2_minus_rank1_source_tag_mean"], .05)
        self.assertAlmostEqual(r["per_axis"]["sun_sin"]["rank2_minus_rank1_source_tag_mean"], -.01)
        self.assertAlmostEqual(r["per_axis"]["q"]["rank2_minus_rank1_source_tag_mean"], .004)
        self.assertTrue(r["release_contains_individual_keys"] is False)
        self.assertNotIn("tagA", json.dumps(r))
        self.assertNotIn("only_initial", json.dumps(r))
        self.assertEqual(r, paired_context_effects(a, b))

    def test_pair_order_not_high_versus_low_weather(self):
        x = {"t"+str(i): {"sun_cos": i/10, "sun_sin": -i/10, "q": i/100}
             for i in range(6)}
        y = {"t"+str(i): {"sun_cos": i/10 + .25, "sun_sin": -i/10 + .75,
                          "q": i/100 - .02} for i in range(6)}
        a = paired_context_effects(x, y)
        b = paired_context_effects(y, x)
        for k in ("sun_cos", "sun_sin", "q"):
            av = a["per_axis"][k]["rank2_minus_rank1_source_tag_mean"]
            bv = b["per_axis"][k]["rank2_minus_rank1_source_tag_mean"]
            self.assertAlmostEqual(av, -bv)

    def test_zero_overlap_is_UNSUPPORTED_not_zero_effect(self):
        left = {"A": {"sun_cos": .2, "sun_sin": .3, "q": -.1}}
        right = {"B": {"sun_cos": -.2, "sun_sin": -.3, "q": .1}}
        r = paired_context_effects(left, right)
        self.assertEqual(r["identical_tag_ids_in_both_original_contexts"], 0)
        self.assertIsNone(r["per_axis"]["sun_cos"]["rank2_minus_rank1_source_tag_mean"])
        self.assertEqual(r["per_axis"]["q"]["status"], "INSUFFICIENT_PAIRED_ANIMALS")
        self.assertIsNone(r["within_paired_tag_rho_delta_q_vs_delta_sun_sin"])

    def test_joint_three_axis_requires_the_same_tags_all_effects(self):
        first = {str(i): {"sun_cos": i/11, "sun_sin": i/12,
                          "q": i/13 if i % 2 else None} for i in range(10)}
        second = {str(i): {"sun_cos": i/11 + .1, "sun_sin": i/12 - .1,
                           "q": i/13 + .2 if i % 2 else None} for i in range(10)}
        r = paired_context_effects(first, second)
        self.assertEqual(r["identical_tag_ids_in_both_original_contexts"], 10)
        self.assertEqual(r["per_axis"]["sun_cos"]["paired_source_animals"], 10)
        self.assertEqual(r["per_axis"]["q"]["paired_source_animals"], 5)
        self.assertEqual(r["all_three_axis_paired_original_tags"], 5)

    def test_original_first_context_coverage_and_number_replay(self):
        orig = self.first["Utah"]["taxon_results"][0]["selected_contexts"][0]
        copy = {
            "rank": orig["selected_gridyear_original_coverage_rank"],
            "source_q_days": orig["original_qualified_q_days"],
            "source_q_animals": orig["original_q_qualified_animals"],
            "preweather_solar_steps": orig["original_moving_daylight_steps"]
        }
        check_context_before_paired(copy, orig)
        verify_first_context_effects(orig, orig)
        copy["source_q_animals"] += 1
        with self.assertRaises(ValueError):
            check_context_before_paired(copy, orig)
        altered = deepcopy(orig)
        altered["direction_sun_sin"]["source_mean"] += .15
        with self.assertRaises(ValueError):
            verify_first_context_effects(altered, orig)

    def test_pair_bootstrap_requires_individual_replicates_and_is_deterministic(self):
        self.assertIsNone(paired_percentiles([1, 2, 3, 4]))
        self.assertEqual(paired_percentiles([0.2] * 6),
                         paired_percentiles([0.2] * 6))
        a, b = paired_percentiles([.2] * 6)
        self.assertAlmostEqual(a, .2)
        self.assertAlmostEqual(b, .2)

    def test_source_legacy_external_negative_direction_verdict_is_locked(self):
        ref = json.loads(BASELINE.read_text(encoding="utf-8"))
        self.assertEqual(ref["frozen_external_radiation_nonreplication_UNCHANGED"],
                         self.contract["external_nonreplication_immutable"])

    def test_ranked_spearman_known_truth(self):
        self.assertAlmostEqual(ranked_spearman([1, 2, 3, 4, 5],
                                               [5, 4, 3, 2, 1]), -1)
        self.assertAlmostEqual(ranked_spearman([1, 2, 3, 4, 5],
                                               [1, 2, 3, 4, 5]), 1)
        self.assertIsNone(ranked_spearman([1, 1, 1, 1, 1],
                                          [1, 2, 3, 4, 5]))


if __name__ == "__main__":
    unittest.main()
