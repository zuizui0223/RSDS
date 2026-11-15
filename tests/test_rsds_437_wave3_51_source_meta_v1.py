"""Source-study aggregation of the frozen Wave3 original GPS case without pseudoreplication."""
import unittest
from copy import deepcopy

from scripts import run_rsds_437_wave3_original_gps_first_effects_v1 as gps
from scripts import aggregate_rsds_437_wave3_51_source_meta_v1 as agg


def mock_aggregate_source(source, index, value=None):
    if source["source_status"] == "SOURCE_METADATA_HOLD":
        status = "SOURCE_FILE_METADATA_HOLD"
        q = None
        verified = False
        tag_status = "SOURCE_UNEVALUATED"
    elif value is None:
        status = "SOURCE_NO_NATIVE_Q_ESTIMABLE"
        q = None
        verified = True
        tag_status = "NATIVE_Q_UNAVAILABLE"
    else:
        status = "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE"
        q = {
            "n_individuals": 6,
            "mean": value,
            "animal_bootstrap_95_percentile": [value-0.01, value+0.01],
        }
        verified = True
        tag_status = "NATIVE_Q_ESTIMABLE"
    return {
        "version": gps.FIRST_RESULT_VERSION,
        "wave3_index": index,
        "wave3_original_source_rank_1_based": index+25,
        "source_package_uuid": source["package_uuid"],
        "source_original_repository_study_id": source["original_repository_study_id"],
        "source_original_DOI_labels": source["original_DOI_labels"],
        "source_metadata_stratum": "single_native",
        "original_registered_native_taxa": source["registered_native_taxa"],
        "original_file_source_status": source["source_status"],
        "original_file_MD5_bytes_verified": verified,
        "source_science_status": status,
        "source_hold_stage": (
            "SOURCE_ONLY_ORIGINAL_FILE_METADATA"
            if source["source_status"]=="SOURCE_METADATA_HOLD" else None
        ),
        "source_hold_reason": (
            source.get("source_hold_reason")
            if source["source_status"]=="SOURCE_METADATA_HOLD" else None
        ),
        "native_group_results": [{
            "taxon": source["registered_native_taxa"][0],
            "q": q,
            "daylength": None,
            "status": tag_status,
        }],
        "native_joint_calendar_covariance": {"status":"NOT_APPLICABLE_SINGLE_NATIVE_SOURCE"},
        "original_source_GPS_event_rows_opened": 0 if not verified else 150,
        "human_mobility_GPS_data_opened": False,
        "old_48_signed_solar_external_nonreplication_UNCHANGED": gps.NEG,
        "old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED": True,
        "original_raw_GPS_dates_tag_IDs_coordinates_exported": False,
    }


class SourceStudyWave3Aggregation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, cls.sources = gps.freeze_original()

    def sample(self):
        rows = [mock_aggregate_source(s,i) for i,s in enumerate(self.sources)]
        rows[0] = mock_aggregate_source(self.sources[0],0,value=0.15)
        rows[1] = mock_aggregate_source(self.sources[1],1,value=-0.03)
        return rows

    def test_fixed_old_35_plus_new_16_without_reassigning_failed_sources(self):
        x = agg.build(self.sample())
        self.assertEqual(x["original_source_study_units"],51)
        self.assertEqual(x["original_previous_source_study_units"],35)
        self.assertEqual(x["new_wave3_source_study_units"],16)
        self.assertEqual(x["original_source_study_unit_all_registered_native_rows"],66)
        self.assertEqual(x["source_complete_all_native_q"],12)
        self.assertEqual(x["source_incomplete_q_count"],39)
        self.assertEqual(x["source_complete_all_native_daylength"],8)
        self.assertEqual(x["source_respondent_only_q_descriptive"]["new_wave3_source_sign_counts"],
                         {"positive":1,"negative":1,"exact_zero":0})
        self.assertEqual(len(x["all_51_source_study_rows"]),51)
        self.assertEqual(len({z["source_original_repository_study_id"]
                              for z in x["all_51_source_study_rows"]}),51)
        self.assertTrue(x["previous_48_source_sunlight_signed_external_not_replicated"])

    def test_all_16_sources_mandatory(self):
        with self.assertRaisesRegex(ValueError,"sixteen"):
            agg.build(self.sample()[:-1])

    def test_no_dropping_or_reclassifying_source_metadata_HOLD(self):
        rows=self.sample()
        rows[2]["source_science_status"]="SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE"
        with self.assertRaisesRegex(ValueError,"HOLD"):
            agg.build(rows)

    def test_no_source_relabeling_or_duplicate(self):
        rows=self.sample()
        rows[0]["source_original_repository_study_id"]=rows[1]["source_original_repository_study_id"]
        with self.assertRaisesRegex(ValueError,"metadata/native"):
            agg.build(rows)

    def test_no_duplicate_prior_source_or_sign_only_acceptance(self):
        rows=self.sample()
        rows[0]["native_group_results"][0]["q"]["n_individuals"]=3
        with self.assertRaisesRegex(ValueError,"Invalid original"):
            agg.build(rows)

    def test_historical_solar_negative_cannot_be_modified(self):
        rows=self.sample()
        rows[0]["old_48_signed_solar_external_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaisesRegex(ValueError,"metadata/native"):
            agg.build(rows)


if __name__=="__main__":
    unittest.main()
