"""Original file metadata-only gate for all six previously frozen source candidates."""
import unittest
from unittest.mock import patch

from scripts import rsds_six_followup_original_file_metadata_v1 as m


def file(size=15000000,checksum="a"*32):
    return {"name":"original_gps.csv","bytes":size,"bitstream_uuid":"original-uuid",
            "checksum":{"checkSumAlgorithm":"MD5","value":checksum}}


class SixFileMetadataBeforeGPS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources,cls.prior=m.freeze()

    def test_original_six_DOI_and_one_old_source_HOLD(self):
        self.assertEqual(len(self.sources),6)
        self.assertEqual(len(self.prior),6)
        self.assertEqual(len(set(x["original_DOI"] for x in self.sources)),6)
        self.assertEqual(sum(p["prior_data_exposure_status"]==
                             "KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE"
                             for p in self.prior),1)

    def test_previous_mallard_metadata_HOLD_never_calls_provider(self):
        i=next(i for i,x in enumerate(self.prior) if x["prior_data_exposure_status"]==
               "KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE")
        with patch.object(m,"original_files",side_effect=RuntimeError("BAD_PROVIDER")) as p:
            r=m.inspect(self.sources[i],self.prior[i])
            p.assert_not_called()
        self.assertEqual(r["source_metadata_status"],"PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED")
        self.assertFalse(r["provider_original_file_list_newly_queried"])
        self.assertEqual(r["original_GPS_events_opened"],0)

    def test_verified_original_file_and_unchanged_known_fox_exposure(self):
        i=next(i for i,x in enumerate(self.prior) if x["prior_data_exposure_status"]==
               "KNOWN_PREVIOUS_ORIGINAL_GPS_DATA_ANALYZED_FOR_SOLAR_DIRECTION")
        with patch.object(m,"original_files",return_value=[file()]):
            r=m.inspect(self.sources[i],self.prior[i])
        self.assertEqual(r["source_metadata_status"],"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE")
        self.assertEqual(r["original_file"]["bytes"],15000000)
        self.assertEqual(r["prior_original_source_exposure_status"],
                         "KNOWN_PREVIOUS_ORIGINAL_GPS_DATA_ANALYZED_FOR_SOLAR_DIRECTION")

    def test_provider_overbudget_remains_source_HOLD(self):
        i=next(i for i,x in enumerate(self.prior) if x["prior_data_exposure_status"]==
               "NO_MATCH_IN_EIGHT_PINNED_PRIOR_RESULT_SURFACES_NOT_EXHAUSTIVE")
        with patch.object(m,"original_files",return_value=[file(size=200000001)]):
            r=m.inspect(self.sources[i],self.prior[i])
        self.assertEqual(r["source_metadata_status"],"SOURCE_METADATA_HOLD")
        self.assertEqual(r["source_metadata_hold_reason"],
                         "FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES")

    def test_full_fixed_denominator_with_source_only_file_reads(self):
        with patch.object(m,"original_files",return_value=[file()]) as fn:
            d=m.evaluate()
        self.assertEqual(d["original_followup_source_study_slots"],6)
        self.assertEqual(len(d["original_source_metadata_rows"]),6)
        self.assertEqual(fn.call_count,5)
        self.assertEqual(d["source_original_metadata_status_counts"],{
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":5,
            "PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED":1})
        self.assertTrue(d["prior_mallard_source_metadata_HOLD_preserved"])
        self.assertTrue(d["prior_fox_original_GPS_read_known"])
        self.assertEqual(d["all_new_original_animal_GPS_events_opened"],0)
        self.assertTrue(d["target_selected_after_63_archive_q_outcomes"])

if __name__=="__main__":
    unittest.main()
