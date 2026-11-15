"""Known-truth, source-locked source-frame unit tests; no GPS events."""
from __future__ import annotations

import copy
import unittest
from unittest.mock import patch

from scripts import rsds_437_systematic_source_frame_v1 as m

class SourceFrameIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.packages,cls.previous=m.inputs()
        cls.result=m.census(cls.contract,cls.packages,cls.previous)

    def test_every_package_and_sensor_metadata_denominator(self):
        summary=self.result["summary"]
        self.assertEqual(summary["all_original_packages"],437)
        self.assertEqual(summary["gps_gnss_packages"],259)
        self.assertEqual(summary["previous_outcome_examined_DOI"],13)
        self.assertEqual(summary["previous_outcome_examined_study_ids"],13)
        self.assertEqual(summary["duplicate_study_ids_among_GPS_packages"],3)
        self.assertEqual(len(self.result["all_437_metadata_status_rows"]),437)
        self.assertEqual(sum(summary["status_counts"].values()),437)
        self.assertEqual(summary["new_animal_GPS_events_or_effects_read"],0)

    def test_metadata_lottery_unique_strata_and_ids(self):
        s=self.result["wave1_selected_in_rank_order"]
        self.assertEqual(len(s),16)
        self.assertEqual(len({x["package_uuid"] for x in s}),16)
        self.assertEqual(len({x["repository_study_ids"][0] for x in s}),16)
        self.assertEqual(len({x["doi_labels"][0] for x in s}),16)
        self.assertEqual(sum(x["selection_stratum"]=="multi_native" for x in s),8)
        self.assertEqual(sum(x["selection_stratum"]=="single_native" for x in s),8)
        self.assertTrue(all(x["stage1_status"]=="FIRST_WAVE_METADATA_CANDIDATE"
                            and x["wave1_selected"] for x in s))

    def test_originally_examined_previous_DOI_and_study_alias_excluded(self):
        old=self.previous
        for row in self.result["all_437_metadata_status_rows"]:
            if any(m.normalize_doi(d) in old for d in row["doi_labels"]):
                self.assertEqual(row["stage1_status"],
                                 "PREVIOUSLY_OUTCOME_EXAMINED_SOURCE_OR_ALIAS")
                self.assertFalse(row["wave1_selected"])
            if row["repository_study_ids"]==["2919708"]:
                self.assertEqual(row["stage1_status"],
                                 "PREVIOUSLY_OUTCOME_EXAMINED_SOURCE_OR_ALIAS")
                self.assertFalse(row["wave1_selected"])

    def test_stable_against_catalog_order_and_titles(self):
        reversed_result=m.census(self.contract,list(reversed(self.packages)),self.previous)
        self.assertEqual(self.result,reversed_result)
        mutated=copy.deepcopy(self.packages)
        for p in mutated:
            p["title"]="An intentionally unrelated changed metadata title"
        self.assertEqual(self.result,m.census(self.contract,mutated,self.previous))

    def test_all_unselected_are_accounted_for_not_zero_effect(self):
        x=self.result["all_437_metadata_status_rows"]
        self.assertEqual(len([a for a in x if a["wave1_selected"]]),16)
        self.assertTrue(all(a["candidate_file_header_or_GPS_event_read"] is False
                            for a in x))
        self.assertNotIn("q_effect",str(self.result))
        self.assertTrue(self.result["old_48_original_source_sunlight_direction_external_NOT_REPLICATED"])
        self.assertTrue(self.result["old_12_original_source_ERA5_raw_JSON_HOLD"])

    def test_full_Git_blob_pins_and_method_consistency(self):
        self.assertEqual(m.run(),self.result)
        bad=copy.deepcopy(self.contract)
        bad["wave1_rank_salt"]="changed-after-outcomes"
        self.assertNotEqual(self.result,m.census(bad,self.packages,self.previous))
        with patch.object(m,"git_blob",return_value="bad"):
            with self.assertRaisesRegex(ValueError,"Previously fixed"):
                m.inputs()

if __name__=="__main__":
    unittest.main()
