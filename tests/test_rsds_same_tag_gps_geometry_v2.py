"""Known-truth GPS-only time vs Sun-direction paired-tag biology tests."""
from __future__ import annotations

import json
import math
import unittest
from datetime import date, timedelta

from scripts.run_rsds_same_tag_gps_geometry_v2 import (
    AXES, check_contract, select_matched_tag_contexts, individual_geometry,
    compare_private_tags
)


class RSDSOriginalTagGPSGeometryTests(unittest.TestCase):
    def test_source_studies_and_historical_weather_holds_are_unchanged(self):
        f = check_contract()
        self.assertEqual(set(f), {"Utah", "Wyoming"})
        self.assertEqual([f[x]["doi"] for x in ("Utah", "Wyoming")],
                         ["10.5441/001/1.711", "10.5441/001/1.703"])
        self.assertEqual(len(f["Utah"]["taxa"]) + len(f["Wyoming"]["taxa"]), 6)

    def test_context_pair_maximizes_shared_original_tags_not_q_values(self):
        def c(key, names, n=10):
            return {"_key":key,
                    "_rows":[{}] * (len(names)*n),
                    "_n":len(names), "_days":len(names)*n,
                    "_count":{animal:n for animal in names}}
        a=c((40,20,2020),list("ABCDE"))
        b=c((41,20,2020),list("ABCDE"))
        d=c((42,20,2020),list("XYZ"),200)
        pair,report=select_matched_tag_contexts([a,d,b])
        self.assertEqual(set(item[0] for item in pair), {a["_key"],b["_key"]})
        self.assertEqual(report["original_q_tags_in_both_selected_contexts"],5)
        self.assertEqual(report["eligible_context_pairs"],3)
        self.assertTrue(report["selection_is_coverage_only"])

    def test_single_context_retains_unavailability(self):
        pair,report=select_matched_tag_contexts([{"_key":(1,1,2020),"_count":{},
                                                "_rows":[],"_n":0,"_days":0}])
        self.assertEqual(pair,[])
        self.assertEqual(report["original_q_tags_in_both_selected_contexts"],0)

    def test_fixed_compass_bearing_is_not_a_joint_solar_response(self):
        # A constant geographic bearing B with changing Sun A produces a
        # varying raw cos(B-A) yet 0 circular covariance by exact identity.
        d0=date(2020,6,1)
        qdays=[{"animal":"private-test-id","date":d0+timedelta(days=i),
                "q":-.05,"daylength":14} for i in range(8)]
        steps=[]
        for i,row in enumerate(qdays):
            for j in range(3):
                steps.append({"animal":"private-test-id","date":row["date"],
                              "bearing":90.0,"sun_azimuth":10.0+i*5+j*8,
                              "solar_hour":9.0+j,
                              "elevation":30.0+j,
                              "daylength":14.0})
        tagged=individual_geometry({"_qdays":qdays,"_events":steps})
        self.assertEqual(len(tagged),1)
        z=next(iter(tagged.values()))
        self.assertEqual(z["_qdays"],8)
        self.assertEqual(z["_sun_steps"],24)
        self.assertAlmostEqual(z["sun_cos"],
                               z["geographic_marginal"]+z["joint_circular"],places=12)
        self.assertAlmostEqual(z["joint_circular"],0,places=12)
        self.assertAlmostEqual(z["q"],-.05)
        self.assertNotIn("private-test-id",json.dumps(z))

    def test_true_paired_tag_not_study_level_mean(self):
        first={f"secret-{i}":{k:float(i) for k in AXES} for i in range(5)}
        second={f"secret-{i}":{k:float(i)+.125 for k in AXES} for i in range(5)}
        first["only-first"]={k:10000.0 for k in AXES}
        second["only-second"]={k:-10000.0 for k in AXES}
        r=compare_private_tags(first,second)
        self.assertEqual(r["same_original_tagged_animals_all_axes"],5)
        self.assertEqual(set(r["source_effects"]),set(AXES))
        for axis in AXES:
            s=r["source_effects"][axis]
            self.assertEqual(s["status"],"PAIRED_TAG_ESTIMABLE")
            self.assertEqual(s["common_original_tagged_animals"],5)
            self.assertAlmostEqual(s["rank2_minus_rank1_tag_mean"],.125)
        self.assertNotIn("secret",json.dumps(r))

    def test_unmatched_animals_do_not_contribute_ecological_zero(self):
        first={"a":{k:1.0 for k in AXES}}
        second={"b":{k:-1.0 for k in AXES}}
        r=compare_private_tags(first,second)
        self.assertEqual(r["same_original_tagged_animals_all_axes"],0)
        self.assertEqual(r["source_effects"]["joint_circular"]["status"],"PAIRED_TAG_UNAVAILABLE")
        self.assertIsNone(r["source_effects"]["joint_circular"]["rank2_minus_rank1_tag_mean"])
        self.assertIsNone(r["within_tag_rho_delta_q_vs_delta_raw_sun_cos"])

    def test_directional_components_are_a_sum_and_do_not_confuse_cos_and_sin(self):
        f={str(i):{k:0.1*i for k in AXES} for i in range(6)}
        g={str(i):{**f[str(i)],"sun_cos":0.1*i+0.05,
                    "geographic_marginal":0.1*i+0.02,
                    "joint_circular":0.1*i+0.03,"sun_sin":0.1*i-0.04}
           for i in range(6)}
        r=compare_private_tags(f,g)
        a=r["source_effects"]
        self.assertAlmostEqual(a["sun_cos"]["rank2_minus_rank1_tag_mean"],
                               a["geographic_marginal"]["rank2_minus_rank1_tag_mean"]
                               +a["joint_circular"]["rank2_minus_rank1_tag_mean"])
        self.assertAlmostEqual(a["sun_sin"]["rank2_minus_rank1_tag_mean"],-.04)


if __name__=="__main__":
    unittest.main()
