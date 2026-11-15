"""Known-truth original native q and DOI-grain safety checks, no GPS source events."""
import unittest

from scripts import rsds_original_source_native_same_species_repeatability_v1 as m


def effect(v,n,low,high):
    return {"mean":v, "n_individuals":n, "animal_bootstrap_95_percentile":[low,high]}


class OriginalSourceNativeRepeatabilityRules(unittest.TestCase):
    def test_normalize_original_source_DOIs(self):
        self.assertEqual(m.norm_doi("doi:10.5441/001/1.304"), "10.5441/001/1.304")
        self.assertEqual(m.norm_doi("https://doi.org/10.5441/001/1.304"),
                         "10.5441/001/1.304")

    def test_source_native_not_multi_native_source_mean(self):
        # Two distinct species in one original source can point in opposite directions;
        # source-average may be zero, but must not be assigned to either species.
        fox=m.q_effect(effect(.25,6,.1,.4))
        goose=m.q_effect(effect(-.25,9,-.4,-.1))
        a=m.make_row({"doi":"doi:10.5441/001/1.291","study_id":"t1",
                      "wave":"wave1","n_native":2},"Vulpes vulpes",fox,m.q_effect(None),True)
        b=m.make_row({"doi":"doi:10.5441/001/1.291","study_id":"t1",
                      "wave":"wave1","n_native":2},"Vulpes lagopus",goose,m.q_effect(None),True)
        self.assertEqual(a["original_native_q"]["q"],.25)
        self.assertEqual(b["original_native_q"]["q"],-.25)
        self.assertEqual(a["original_DOI"],b["original_DOI"])
        self.assertEqual(a["registered_native_count"],2)

    def test_original_tag_bootstrap_sign_vs_mean_sign(self):
        e=m.q_effect(effect(-0.03,8,-0.09,0.02))
        self.assertEqual(e["point_sign"],"NEGATIVE")
        self.assertEqual(e["nominal_source_CI_sign"],"INCLUDES_ZERO")
        self.assertEqual(m.q_effect(effect(.1,5,.01,.2))["nominal_source_CI_sign"],"POSITIVE")

    def test_legacy_and_new_numerical_shapes(self):
        old=m.q_effect({"mean":-.1,"source_tags":20,
                        "tag_bootstrap95":[-.2,-.01]},old=True)
        new=m.q_effect(effect(-.1,20,-.2,-.01))
        self.assertEqual(old["nominal_source_CI_sign"],"NEGATIVE")
        self.assertEqual(old["q"],new["q"])
        self.assertTrue(old["available"])
        self.assertEqual(old["original_source_tagged_animals"],20)

    def test_missing_not_zero_and_cannot_infer_ci(self):
        self.assertFalse(m.q_effect(None)["available"])
        self.assertIsNone(m.q_effect(None)["q"])
        self.assertEqual(m.q_effect(None)["nominal_source_CI_sign"],"UNAVAILABLE")

    def test_invalid_tag_support_ci_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"Invalid original"):
            m.q_effect(effect(1,4,.2,1.5))
        with self.assertRaisesRegex(ValueError,"Invalid original"):
            m.q_effect(effect(1,6,1.5,.2))

    def test_source_wave3_original_ci_replay(self):
        obj={"available":True,"mean":-.02,"n_individuals":6,
             "bootstrap_95_percentile":[-.03,.001]}
        self.assertEqual(m.from_wave3_ci(obj)["nominal_source_CI_sign"],"INCLUDES_ZERO")
        self.assertFalse(m.from_wave3_ci({"available":False})["available"])


if __name__=="__main__":
    unittest.main()
