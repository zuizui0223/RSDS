"""Archived study native heat q signs: all 6 originals and no fabricated meta."""
import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import unittest

from scripts import rsds_ungulate_native_Tmax_archived_six_audit_v1 as m


class TestOriginalUngulateSixTaxa(unittest.TestCase):
    def test_exact_archived_six_source_native_ungulates(self):
        d=m.audit()
        self.assertEqual(d["distinct_original_field_study_DOIs"],2)
        self.assertEqual(d["original_source_native_taxon_study_rows"],6)
        self.assertEqual(d["source_native_q_Tmax_nominal95_CI_sign_counts"],
                         {"INCLUDES_ZERO":2,"NEGATIVE":1,"POSITIVE":3})
        self.assertEqual(d["within_Utah_negative_species"],["Ovis canadensis"])
        self.assertEqual(d["within_Utah_positive_species"],
                         ["Cervus elaphus","Odocoileus hemionus"])
        self.assertTrue(d["same_source_does_not_establish_same_grid_month_animal_cohort"])
        self.assertTrue(d["source_bootstrap_not_interaction_CI_or_world_meta_test"])
        self.assertEqual(d["original_GPS_and_weather_event_rows_reopened"],0)

    def test_original_source_signed_artifact_altered_fails_closed(self):
        a=json.loads(m.PATH.read_text(encoding="utf-8"))
        a["source_provenance"]["utah"]["original_artifact_digest"]="sha256:unverified"
        with TemporaryDirectory() as td:
            p=Path(td)/"altered.json";p.write_text(json.dumps(a))
            with patch.object(m,"PATH",p):
                with self.assertRaisesRegex(ValueError,"artifact provenance"):
                    m.audit()

    def test_cherry_picking_or_reclassifying_source_sign_fails(self):
        a=json.loads(m.PATH.read_text(encoding="utf-8"))
        a["original_tmax_q_native_rows"]=a["original_tmax_q_native_rows"][:5]
        with TemporaryDirectory() as td:
            p=Path(td)/"altered.json";p.write_text(json.dumps(a))
            with patch.object(m,"PATH",p):
                with self.assertRaisesRegex(ValueError,"missing"):
                    m.audit()

if __name__=="__main__":unittest.main()
