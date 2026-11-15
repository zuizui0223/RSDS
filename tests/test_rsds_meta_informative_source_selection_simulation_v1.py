"""Known truth controls for source selection, not effects on actual animals."""
from __future__ import annotations

import copy
import math
import unittest

from scripts import rsds_meta_informative_source_selection_simulation_v1 as m


class TrueSourceFrameSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=m.freeze()

    def test_real_first_source_evidence_guard_and_frame(self):
        c=self.c
        self.assertEqual(c["original_source_frame"],{
            "single_native":109,"multi_native":11})
        self.assertEqual(c["hypothetical_combined_sample"],{
            "single_native":24,"multi_native":11})
        self.assertEqual(c["repetitions_per_scenario"],1000)
        self.assertEqual(c["true_mean_log1p_night_minus_day_speed_q"],.04)
        self.assertEqual(
            [r["beta_q"] for r in c["source_q_response_dependence_scaled_betas"]],
            [0.0,-1.5,1.5])

    def test_stable_original_world_truth_and_sequential_wave(self):
        w=m.generate_world(self.c,7)
        self.assertEqual(len(w["single_native"]),109)
        self.assertEqual(len(w["multi_native"]),11)
        self.assertEqual(w,m.generate_world(self.c,7))
        first=m.sample(self.c,w,1)
        second=m.sample(self.c,w,2)
        for g in m.GROUPS:
            self.assertEqual(first[g],
                second[g][:self.c["hypothetical_wave1_sample"][g]])
        self.assertTrue(all(.35<=z["cadence"]<=.95
                            for units in w.values() for z in units))
        self.assertTrue(math.isfinite(m.target(self.c,w)))

    def test_MAR_response_has_no_q_dependence_conditional_on_cadence(self):
        c=self.c
        q0=m.probability(c,"single_native",.7,-.2,0)
        q1=m.probability(c,"single_native",.7,+.3,0)
        self.assertEqual(q0,q1)
        self.assertLess(m.probability(c,"single_native",.7,+.2,-1.5),
                        m.probability(c,"single_native",.7,-.2,-1.5))
        self.assertGreater(m.probability(c,"single_native",.7,+.2,+1.5),
                           m.probability(c,"single_native",.7,-.2,+1.5))
        self.assertAlmostEqual(
            m.probability(c,"single_native",.7,+.1,0),
            m.marginal_probability(c,"single_native",.7,0),
            places=12)

    def test_oracle_propensity_only_or_cadence_only_both_known_truth(self):
        c=self.c
        artificial={
            "single_native":[(.4,True,.5,.25),(-.2,False,.5,.25)],
            "multi_native":[(.1,True,.5,.25),(.0,False,.5,.25)]
        }
        oracle=m.HT(c,artificial,2)
        marg=m.HT(c,artificial,3)
        self.assertAlmostEqual(oracle,(109*.4+11*.1)/120,places=12)
        self.assertAlmostEqual(marg,2*oracle,places=12)
        self.assertIsNone(m.poststratified(c,{
            "single_native":[(.1,False,.5,.5)],
            "multi_native":[(.1,True,.5,.5)]
        },lambda a:a[1]))

    def test_same_truth_three_different_measurability_mechanisms(self):
        c=copy.deepcopy(self.c)
        c["repetitions_per_scenario"]=220
        result=m.analyze(c)
        a=result["synthetic_scenarios"]
        self.assertEqual(len(a),3)
        neutral=a["cadence_only_MAR_in_synthetic_world"]
        nmar_low=a["outcome_linked_negative_MNAR_synthetic_world"]
        nmar_high=a["outcome_linked_positive_MNAR_synthetic_world"]
        assert_name="respondent_only_35"
        self.assertLess(abs(neutral["estimators"][assert_name][
            "bias_vs_actual_finite_120_source_truth"]),.035)
        self.assertLess(nmar_low["estimators"][assert_name][
            "bias_vs_actual_finite_120_source_truth"],-.03)
        self.assertGreater(nmar_high["estimators"][assert_name][
            "bias_vs_actual_finite_120_source_truth"],+.03)
        self.assertTrue(all(v["estimators"][assert_name][
            "worlds_missing_estimator"]>0 for v in a.values()))
        self.assertTrue(all(v["estimators"][
            "true_source_q_response_probability_oracle_HT"][
                "worlds_with_defined_estimator"]==220 for v in a.values()))
        self.assertTrue(result["source_truth_and_cadence_shared_across_all_three_scenarios"])
        self.assertTrue(result["source_file_MD5_ORIGINAL_35_GPS_NOT_ACCESSED_BY_THIS_SIMULATION"])
        self.assertTrue(result["old_48_source_shortwave_direction_external_NOT_REPLICATED"])

    def test_no_meta_claim_from_zero_response_stratum(self):
        c=self.c
        a={"single_native":[(.02,False,.5,.5)],
           "multi_native":[(-.03,False,.25,.25)]}
        self.assertIsNone(m.poststratified(c,a,lambda v:v[1]))
        self.assertEqual(m.HT(c,a,2),0.0)
        # 0 in oracle HT is unbiased contribution of the missing responses,
        # NOT imputed biological q: it is a design estimator with known p.


if __name__=="__main__":
    unittest.main()
