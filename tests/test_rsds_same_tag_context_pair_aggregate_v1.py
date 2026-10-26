"""Source-complete aggregate tests for RSDS paired original context effects."""
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_same_tag_context_pair_v1 import aggregate


def fake_source(region):
    taxa = ["Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis"] if region == "Utah" else ["Odocoileus hemionus","Ovis canadensis"]
    doi = "10.5441/001/1.711" if region == "Utah" else "10.5441/001/1.703"
    return {
        "version":"RSDS_SAME_TAG_CONTEXT_PAIR_V1",
        "region":region,
        "original_source_doi":doi,
        "state":"SOURCE_PAIRED_CONTEXT_REASSESSMENT_COMPLETE",
        "source_taxa":[{
            "taxon":t,
            "original_source_replay_success":True,
            "new_paired_sensitivity_status":"INSUFFICIENT_PAIRED_ANIMALS",
            "frozen_first_contexts_verified":True,
            "paired_tag_effect":None,
            "context_coverage":[{"status":"FIRST_CONTEXT_NUMERICAL_REPLAY_PASS"},
                                {"status":"FIRST_CONTEXT_NUMERICAL_REPLAY_PASS"}]
        } for t in taxa]
    }


class FixedSourceOriginalPairedAggregate(unittest.TestCase):
    def test_full_source_frame_keeps_all_null_not_zero(self):
        d = aggregate(fake_source("Utah"), fake_source("Wyoming"))
        self.assertEqual(d["source_original_studies"],2)
        self.assertEqual(d["native_study_taxon_rows"],6)
        self.assertEqual(d["source_taxa_with_two_comparable_original_contexts"],0)
        self.assertEqual(d["number_of_source_taxa_with_estimate_by_axis"],{
            "sun_cos":0,"sun_sin":0,"q":0
        })
        self.assertTrue(all(row["paired_same_tag_effect"] is None
                            for row in d["all_native_source_taxon_results"]))

    def test_retains_exact_one_valid_tag_source_and_no_unsupported_overclaim(self):
        a,b = fake_source("Utah"),fake_source("Wyoming")
        a["source_taxa"][0]["paired_tag_effect"]={
            "per_axis":{k:{
                "paired_source_animals":6,
                "status":"PAIRED_ESTIMABLE",
                "rank2_minus_rank1_source_tag_mean":1,
                "paired_animal_bootstrap95":[0.5,1.5]
            } for k in ("sun_cos","sun_sin","q")}
        }
        r=aggregate(a,b)
        self.assertEqual(r["number_of_source_taxa_with_estimate_by_axis"],{
            "sun_cos":1,"sun_sin":1,"q":1})
        self.assertEqual(r["source_taxa_with_two_comparable_original_contexts"],1)
        a["source_taxa"][0]["paired_tag_effect"]["per_axis"]["q"]["paired_source_animals"]=4
        with self.assertRaises(ValueError):
            aggregate(a,b)

    def test_cannot_substitute_a_new_study_or_missing_taxon(self):
        a,b=fake_source("Utah"),fake_source("Wyoming")
        b["original_source_doi"]="10.5441/001/1.711"
        with self.assertRaises(ValueError):
            aggregate(a,b)
        b=fake_source("Wyoming")
        b["source_taxa"].pop()
        with self.assertRaises(ValueError):
            aggregate(a,b)

    def test_previous_negative_directional_result_not_reclassified(self):
        a,b=fake_source("Utah"),fake_source("Wyoming")
        d=aggregate(a,b)
        self.assertEqual(d["heldout_negative_radiation_nonreplication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")


if __name__=="__main__":
    unittest.main()
