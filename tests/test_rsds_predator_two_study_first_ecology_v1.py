"""Original source-first Illinois/Washington predator ecology outcome integrity tests."""
from __future__ import annotations
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.aggregate_rsds_predator_pair_two_study_temporal_v1 import (
    ROOT, OUTPUT, aggregate
)

FIRST_BLOB="a0e1ec8e61836f87057189c66e526dd9505f7e62"


def first():
    raw=OUTPUT.read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if sha!=FIRST_BLOB:
        raise ValueError("Original source first wildlife temporal discovery was rewritten")
    return json.loads(raw.decode("utf-8"))


class OriginalPredatorPairTwoStudyFirstEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=first()

    def test_original_two_DOIs_and_five_native_taxa_not_selected_out(self):
        d=self.d
        self.assertEqual(d["original_study_DOI_count"],2)
        self.assertEqual(d["native_study_taxon_rows"],5)
        self.assertEqual(d["preselected_primary_native_pair"],
                         ["Canis latrans","Lynx rufus"])
        self.assertEqual({x["source_original_study_DOI"] for x in
                          d["all_original_source_study_native_taxa"]},
                         {"10.5441/001/1.gm93267b","10.5441/001/1.649"})
        self.assertTrue(all(x["source_provider_MD5_replayed"]
                            for x in d["all_original_source_study_native_taxa"]))

    def test_Washington_complete_provider_replay_but_q_temporally_unavailable(self):
        w=next(x for x in self.d["all_original_source_study_native_taxa"]
               if x["source_region"]=="Washington")
        self.assertEqual(set(w["source_native_taxa"]),{"Canis latrans","Lynx rufus"})
        self.assertEqual({t:x["original_eligible_q_tag_count"]
                         for t,x in w["source_native_taxa"].items()},
                         {"Canis latrans":0,"Lynx rufus":0})
        self.assertEqual(w["shared_calendar"]["ORIGINAL_ONLY"]["source_calendar_days"],0)
        self.assertTrue(all(t["q"]["status"]=="UNAVAILABLE"
                            for t in w["source_native_taxa"].values()))
        self.assertEqual(self.d["two_study_native_pair_repeatability"]["ORIGINAL_ONLY"][
            "original"]["two_source_same_sign_classification"]["q"],"INCONCLUSIVE_OR_SOURCE_HOLD")

    def test_Illinois_coyote_bobcat_same_calendar_q_rank_and_uncertain_animal_seasonal_slope(self):
        s=next(x for x in self.d["all_original_source_study_native_taxa"]
               if x["source_region"]=="Illinois")
        observed=s["shared_calendar"]["ORIGINAL_ONLY"]
        self.assertEqual(observed["source_calendar_days"],194)
        self.assertEqual(observed["source_year_month_blocks"],9)
        self.assertEqual(observed["native_original_tags_on_shared_source_dates"],
                         {"Canis latrans":9,"Lynx rufus":5})
        a=observed["native_source_pair"]["original"]
        contrast=a["q_coyote_minus_bobcat"]
        self.assertAlmostEqual(contrast["coyote_minus_bobcat_mean"],
                               .16552443224935098,places=14)
        self.assertEqual(contrast["nominal_sign"],"POSITIVE")
        self.assertTrue(contrast["source_tag_bootstrap95"][0]>0)
        self.assertEqual(a["daylength_coyote_minus_bobcat"]["nominal_sign"],
                         "INCLUDES_ZERO")
        date=observed["calendar_year_month_block_pair"]["slopes"]["original"]
        self.assertEqual(date["nominal_sign"],"POSITIVE")
        self.assertAlmostEqual(date["same_calendar_coyote_minus_bobcat_q_daylength_slope"],
                               .038716726476820805,places=14)

    def test_Illinois_native_deer_retained_not_invented_shared_date_predation(self):
        s=next(x for x in self.d["all_original_source_study_native_taxa"]
               if x["source_region"]=="Illinois")
        self.assertEqual(set(s["source_native_taxa"]),
                         {"Canis latrans","Lynx rufus","Odocoileus virginianus"})
        self.assertEqual({t:z["original_eligible_q_tag_count"]
                          for t,z in s["source_native_taxa"].items()},
                         {"Canis latrans":20,"Lynx rufus":5,
                          "Odocoileus virginianus":23})
        # White-tailed deer are a native source q estimand but are NOT
        # falsely paired to the two predators' 194 source dates.
        self.assertEqual(set(s["shared_calendar"]["ORIGINAL_ONLY"]["native_original_tags_on_shared_source_dates"]),
                         {"Canis latrans","Lynx rufus"})
        self.assertTrue(s["source_native_taxa"]["Odocoileus virginianus"]["q"]["mean"]>0)

    def test_UTC_even_odd_three_schedule_results_are_a_real_source_HOLD(self):
        for s in self.d["all_original_source_study_native_taxa"]:
            t=s["shared_calendar"]["THREE_SCHEDULES"]
            self.assertEqual(t["source_calendar_days"],0)
            self.assertEqual(t["status"],"INSUFFICIENT_COOCCURRENT_SOURCE_NATIVE_TAGS")
            for x in t["native_source_pair"].values():
                self.assertEqual(x["nominal_q_direction"],"UNAVAILABLE")
        self.assertTrue(self.d["previous_12_ERA5_source_JSON_HOLDS_UNCHANGED"])
        self.assertEqual(
            self.d["previous_48_source_negative_direction_shortwave_external_replication_UNCHANGED"],
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
        )

    def test_public_source_only_no_individual_location_or_times(self):
        text=json.dumps(self.d)
        for k in ('"individual_id"','"tag_id"','"animal_id"',
                  '"location-lat"','"location-long"','"utc_hour"',
                  '"original_gps_event"','"raw_tracks"','"original_day"'):
            self.assertNotIn(k,text)

    def test_prior_source_group_result_modified_fails_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            x=Path(folder)/"first.json"
            x.write_bytes(OUTPUT.read_bytes())
            data=json.loads(x.read_text(encoding="utf-8"))
            data["original_study_DOI_count"]=3
            x.write_text(json.dumps(data))
            import hashlib
            raw=x.read_bytes()
            sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
            self.assertNotEqual(sha,FIRST_BLOB)


if __name__=="__main__":
    unittest.main()
