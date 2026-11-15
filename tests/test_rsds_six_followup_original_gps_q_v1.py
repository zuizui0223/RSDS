"""Frozen source IDs, previous GPS exposure and original q-support tests."""
import unittest
from unittest.mock import patch
from scripts import run_rsds_six_same_species_followup_original_gps_q_v1 as runner
from scripts import aggregate_rsds_six_followup_original_gps_q_v1 as agg

class SixTargetFirstOriginalGPS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.sources=runner.freeze_original()

    def test_all_six_with_only_four_first_source_download_candidates(self):
        self.assertEqual(len(self.sources),6)
        self.assertEqual(sum(x["source_metadata_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE" for x in self.sources),4)
        self.assertEqual(sum(x["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE" for x in self.sources),2)
        self.assertEqual(len(set(x["original_source_DOI"] for x in self.sources)),6)
        self.assertEqual(self.contract["source_q_primary_prediction_sign_by_DOI"],{
            "10.5441/001/1.46b0mq21":"NEGATIVE",
            "10.5441/001/1.708":"NEGATIVE",
            "10.5441/001/1.72hh609t":"POSITIVE",
            "10.5441/001/1.ps244r11":"NEGATIVE"})

    def test_original_mallard_and_big_file_HOLD_do_not_open_provider_GPS(self):
        hold=[i for i,r in enumerate(self.sources) if r["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"]
        self.assertEqual(hold,[3,4])
        with patch.object(runner,"original_files",side_effect=RuntimeError("ORIGINAL_GPS_FORBIDDEN")) as p:
            for index in hold:
                x=runner.process(index)
                self.assertEqual(x["source_q_science_status"],"SOURCE_FILE_METADATA_HOLD")
                self.assertEqual(x["original_animal_GPS_event_rows_opened"],0)
                self.assertFalse(x["original_file_MD5_bytes_verified"])
                self.assertIsNone(x["native_group_results"][0]["q"])
            p.assert_not_called()

    def test_md5_stale_provider_fails_closed_not_reselect(self):
        i=0
        with patch.object(runner,"original_files",side_effect=OSError("REMOTE_URL")):
            x=runner.process(i)
        self.assertEqual(x["source_q_science_status"],"SOURCE_GPS_HOLD")
        self.assertEqual(x["source_hold_stage"],"ORIGINAL_BITSTREAM_UUID_BYTES_MD5")
        self.assertEqual(x["source_hold_error_type"],"OSError")
        self.assertIsNone(x["native_group_results"][0]["q"])
        self.assertNotIn("REMOTE_URL",str(x))

    def test_invalid_original_source_slot_fails(self):
        for n in [-1,6,7]:
            with self.assertRaisesRegex(ValueError,"Only six"):
                runner.process(n)

class SixTargetAggregatedSourceTruth(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.sources=runner.freeze_original()

    def make(self,index,mean=None):
        source=self.sources[index]
        h=source["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
        q=None if mean is None else {
            "mean":mean,"n_individuals":8,
            "animal_bootstrap_95_percentile":[mean-0.01,mean+0.01]
        }
        return {
            "version":runner.RESULT_VERSION,
            "source_index":index,
            "source_original_DOI":source["original_source_DOI"],
            "source_original_repository_study_id":source["original_repository_study_id"],
            "source_original_public_package_uuid":source["original_public_package_uuid"],
            "source_original_registered_native_taxa":source["registered_original_native_labels"],
            "source_original_file_metadata_status":source["source_metadata_status"],
            "prior_original_source_exposure_status":source["prior_original_source_exposure_status"],
            "prior_original_archive_selected_prediction_sign":self.contract["source_q_primary_prediction_sign_by_DOI"].get(source["original_source_DOI"]),
            "original_file_MD5_bytes_verified":not h,
            "source_q_science_status":"SOURCE_FILE_METADATA_HOLD" if h else ("SOURCE_Q_ESTIMABLE" if q else "SOURCE_Q_UNAVAILABLE_NOT_ZERO"),
            "source_hold_stage":"ORIGINAL_PROVIDER_PRE_EVENT_METADATA" if h else None,
            "source_hold_error_type":source["source_metadata_hold_reason"] if h else None,
            "native_group_results":[{"taxon":source["registered_original_native_labels"][0],
                                     "q":q,"daylength":None,"status":"SOURCE_UNEVALUATED"}],
            "original_animal_GPS_event_rows_opened":0 if h else 200,
            "no_raw_animal_GPS_IDs_coordinates_exact_dates_or_tag_effects_exported":True,
            "earlier_original_48_solar_direction_signed_external_not_replicated":True,
            "earlier_original_12_ERA5_source_SHA_HOLD":True
        }

    def fixture(self):
        # One outcome per original source, with both favorable and unfavorable examples.
        return [
            self.make(0,mean=-.10),
            self.make(1,mean=+.08), # predeclared negative, supported contradiction
            self.make(2,mean=+.11), # predeclared positive, but already GPS-exposed in solar project
            self.make(3),
            self.make(4),
            self.make(5,mean=-.03)
        ]

    def test_known_truth_sign_decision_and_holds(self):
        z=agg.aggregate(self.fixture())
        self.assertEqual(z["original_followup_source_DOIs_all"],6)
        self.assertEqual(z["new_original_source_q_supported"],4)
        self.assertEqual(z["new_source_effect_repeatability_nominal_verdict_counts"],{
            "NOMINAL_SUPPORTED_SAME_SIGN":3,
            "NOMINAL_SUPPORTED_OPPOSITE_SIGN":1,
            "UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD":2
        })
        self.assertEqual(len(z["all_six_original_source_taxon_results_including_HOLD"]),6)
        self.assertTrue(z["new_source_study_IDs_selected_after_previous_same_species_results"])
        self.assertTrue(z["source_taxon_same_original_tag_population_independence_unverified"])

    def test_preserve_all_six_and_frozen_source_labels(self):
        with self.assertRaisesRegex(ValueError,"six"):
            agg.aggregate(self.fixture()[:-1])
        z=self.fixture()
        z[0]["source_original_DOI"]=z[1]["source_original_DOI"]
        with self.assertRaisesRegex(ValueError,"Original source"):
            agg.aggregate(z)

    def test_reclassify_previous_SOURCE_HOLD_fails_closed(self):
        z=self.fixture()
        z[3]["source_q_science_status"]="SOURCE_Q_ESTIMABLE"
        with self.assertRaisesRegex(ValueError,"HOLD"):
            agg.aggregate(z)

    def test_nominal_ci_point_support_distinct(self):
        e={"available":True,"point_sign":"NEGATIVE","nominal_source_CI_sign":"INCLUDES_ZERO"}
        self.assertEqual(agg.compare(e,"NEGATIVE"),"SAME_POINT_SIGN_CI_INCLUDES_ZERO")
        self.assertEqual(agg.compare(e,"POSITIVE"),"OPPOSITE_POINT_SIGN_CI_INCLUDES_ZERO")

if __name__=="__main__":
    unittest.main()
