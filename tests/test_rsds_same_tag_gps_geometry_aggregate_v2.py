"""Contract/denominator tests for RSDS paired original GPS-only result."""
import unittest
from copy import deepcopy
from scripts.aggregate_rsds_same_tag_gps_geometry_v2 import aggregate


def fake(region):
    taxa = ["Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis"] if region=="Utah" else ["Odocoileus hemionus","Ovis canadensis"]
    return {
        "version":"RSDS_SAME_TAG_GPS_GEOMETRY_V2","region":region,
        "source_doi":"10.5441/001/1.711" if region=="Utah" else "10.5441/001/1.703",
        "status":"ORIGINAL_SOURCE_GPS_GEOMETRY_COMPLETE",
        "original_provider_and_q_daylight_replay_before_heading":True,
        "source_taxa":[{
            "taxon":t,"original_source_replay_success":True,
            "geometry_source_status":"PAIRED_TAG_UNAVAILABLE",
            "selection_is_coverage_only":True,
            "all_original_q_gridyears":8,"eligible_q_coverage_gridyears":5,
            "eligible_context_pairs":10,"original_q_tags_in_both_selected_contexts":4,
            "paired_context_effect":None
        } for t in taxa]
    }


class OriginalGPSGeometryAtlasTest(unittest.TestCase):
    def test_all_original_studies_and_unsupported_not_zero(self):
        r=aggregate({"Utah":fake("Utah"),"Wyoming":fake("Wyoming")})
        self.assertEqual(r["source_taxon_rows"],6)
        self.assertEqual(r["original_source_count"],2)
        self.assertEqual(r["source_taxa_with_paired_tag_estimates"],0)
        self.assertTrue(r["previous_raw_ERA5_SHA_HOLD_unchanged"])
        self.assertTrue(all(x["source_geometry_paired"] is None
                            for x in r["all_original_source_taxon_results"]))
        self.assertTrue(all(v["axes"]["q"]["status"]=="NOT_ESTIMABLE_TWO_SOURCES"
                            for v in r["two_region_same_native_taxon_repeatability"]))

    def test_no_invention_of_new_species_or_overclaim(self):
        a,b=fake("Utah"),fake("Wyoming")
        b["source_taxa"].pop()
        with self.assertRaises(ValueError):
            aggregate({"Utah":a,"Wyoming":b})
        b=fake("Wyoming"); b["source_doi"]="10.5441/001/1.999"
        with self.assertRaises(ValueError):
            aggregate({"Utah":a,"Wyoming":b})

    def test_original_tag_n_must_support_pair_effect(self):
        a,b=fake("Utah"),fake("Wyoming")
        source={}
        for axis in ("q","sun_cos","sun_sin","joint_circular","geographic_marginal"):
            source[axis]={"status":"PAIRED_TAG_ESTIMABLE",
                          "common_original_tagged_animals":4,
                          "rank2_minus_rank1_tag_mean":.1}
        a["source_taxa"][0]["paired_context_effect"]={
            "same_original_tagged_animals_all_axes":4,
            "source_effects":source
        }
        with self.assertRaises(ValueError):
            aggregate({"Utah":a,"Wyoming":b})


if __name__=="__main__":
    unittest.main()
