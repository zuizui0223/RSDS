"""Known bootstrap sign-support rules for frozen source-level Wave3 post-outcome audit."""
import unittest

from scripts import rsds_wave3_first_source_q_ci_audit_v1 as m


def effect(n, mu, lo, hi):
    return {
        "n_individuals": n,
        "mean": mu,
        "animal_bootstrap_95_percentile": [lo,hi],
    }


class SourceIntervalAuditRules(unittest.TestCase):
    def test_unavailable_is_not_zero(self):
        a=m.bootstrap_interval_status(None)
        self.assertFalse(a["available"])
        self.assertIsNone(a["mean"])
        self.assertEqual(a["nominal_bootstrap_sign"],"UNAVAILABLE")

    def test_negative_and_positive_bootstrap(self):
        self.assertEqual(m.bootstrap_interval_status(
            effect(6,-0.1,-0.2,-0.01))["nominal_bootstrap_sign"],"NOMINAL_NEGATIVE")
        self.assertEqual(m.bootstrap_interval_status(
            effect(5,0.1,0.01,0.2))["nominal_bootstrap_sign"],"NOMINAL_POSITIVE")

    def test_point_negative_not_supported_when_bootstrap_crosses_zero(self):
        a=m.bootstrap_interval_status(effect(6,-0.001,-0.02,0.03))
        self.assertEqual(a["nominal_point_sign"],"NEGATIVE")
        self.assertEqual(a["nominal_bootstrap_sign"],"INCLUDES_ZERO")

    def test_source_tag_count_and_interval_integrity_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"Incorrect original"):
            m.bootstrap_interval_status(effect(4,0.1,0.01,0.2))
        with self.assertRaisesRegex(ValueError,"Incorrect original"):
            m.bootstrap_interval_status(effect(6,0.1,0.2,0.01))

    def test_signed_slopes_cannot_imply_solar_causal_response(self):
        self.assertEqual(m.sign(-1),"NEGATIVE")
        self.assertEqual(m.sign(1),"POSITIVE")
        self.assertEqual(m.sign(0),"EXACT_ZERO")


if __name__=="__main__":
    unittest.main()
