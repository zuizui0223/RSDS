"""All four + two original native source temperature comparison contract tests."""
from __future__ import annotations
from datetime import date,timedelta
import json
from unittest.mock import patch
import unittest

from scripts import rsds_two_source_all_native_ungulate_month_heat_v1 as m
from scripts import aggregate_rsds_two_source_six_native_month_heat_v1 as a
from scripts.movebank_native_ecology_20261019 import daylength


def source_hold(s,i):
    return {
        "version":m.SOURCE_VERSION,
        "source_index":i,
        "original_study_DOI":s["doi"],
        "original_repository_study_id":s["id"],
        "registered_native_taxa":s["taxa"],
        "source_status":"SOURCE_HOLD",
        "source_native_group_results":[],
        "source_ERA5_grid_count":0,"source_ERA5_grid_HOLDs":0,
        "source_ERA5_daily_physical_SHA256":None,
        "source_new_weather_is_not_previous_20261023_or_20261116_physical_replay":True,
        "previous_48_source_irradiance_heading_direction_NOT_REPLICATED":True,
        "previous_12_source_ERA5_raw_SHA_HOLD":True,
        "no_original_animal_ID_tracked_coordinates_exact_date_GPS_weather_grid_or_individual_effect_exported":True
    }


class TwoOriginalSourceSixNative(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c,cls.sources,cls.baseline=m.frozen()

    def test_source_IDs_six_natives_every_species_registered(self):
        self.assertEqual(len(self.sources),2)
        self.assertEqual([len(s["taxa"]) for s in self.sources],[4,2])
        self.assertEqual(sum(len(s["taxa"]) for s in self.sources),6)
        self.assertEqual([s["doi"] for s in self.sources],
                         ["10.5441/001/1.711","10.5441/001/1.703"])
        self.assertEqual(self.sources[0]["reader"],"UTAH_ZIP")
        self.assertEqual(self.sources[1]["reader"],"SOURCE_CSV")
        self.assertEqual(self.sources[0]["taxa"],
                         ["Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis"])
        self.assertEqual(self.sources[1]["taxa"],["Odocoileus hemionus","Ovis canadensis"])
        self.assertTrue(self.c["source_inference"].startswith("Two source DOI studies"))

    def test_frozen_source_q_and_daylength_baselines_all_six_nonzero_support(self):
        self.assertEqual(len(self.baseline),2)
        for s in self.sources:
            self.assertEqual(set(s["taxa"]),set(self.baseline[s["doi"]]))
            for taxon in s["taxa"]:
                p=self.baseline[s["doi"]][taxon]
                self.assertGreaterEqual(p["q"]["source_tags"],5)
                self.assertGreaterEqual(p["daylength"]["source_tags"],5)
                self.assertIsNotNone(p["q"]["mean"])
                self.assertIsNotNone(p["daylength"]["mean"])

    def test_no_weather_missingness_filled_from_source_original(self):
        rows={"Odocoileus hemionus":[
            ("OriginalTag",2020,[(date(2020,3,3),11,-.1,(40,-100),(2020,3))])]}
        with patch.object(m,"climate_at_gps_day",return_value=None):
            with self.assertRaisesRegex(ValueError,"complete past-30d"):
                m.source_eligible_climate_or_HOLD(rows,{(40,-100):{}})

    def test_original_q_component_temperature_identity_on_synthetic_tag_year(self):
        from scripts import rsds_three_mule_within_month_joint_weather_v1 as w
        from scripts import rsds_three_mule_day_night_heat_components_v1 as original
        rows=[];weather={}
        for mo in (2,4,10):
            for i in range(16):
                d=date(2020,mo,3)+timedelta(days=i)
                L=daylength(42,d)
                Tmax=12+i*.61+mo/10
                W=.6+((i*7)%9)*.08+mo/100
                DAY=.3-.004*Tmax+.015*L-.03*W+mo*.03
                NIGHT=.6+.007*Tmax-.021*L+.02*W-mo*.03
                rows.append((d,L,NIGHT-DAY,(40,-100),(2020,mo),DAY,NIGHT))
                weather[d]=(Tmax-1.2,Tmax,W)
        q=[("tag"+str(i),2020,[r[:5] for r in rows]) for i in range(6)]
        parts=[("tag"+str(i),2020,rows) for i in range(6)]
        with patch.object(w,"climate_at_gps_day",side_effect=lambda _,day:weather[day]):
            first=w.per_source_same_sample_models(q,{(40,-100):{}},20261117)
        fake={"new_month_same_tag_model":first}
        with patch.object(original,"climate_at_gps_day",side_effect=lambda _,day:weather[day]):
            r=original.summarize_original_tag_components(parts,{(40,-100):{}},fake,20261117)
        self.assertEqual(r["source_status"],"SOURCE_COMPONENTS_ESTIMABLE")
        self.assertTrue(r["q_night_minus_day_algebra_exact_in_all_tag_years_and_bootstrap"])
        e=r["components"]
        self.assertAlmostEqual(e["month_q_tmax"]["estimate"],
                               e["month_night_tmax"]["estimate"]-e["month_day_tmax"]["estimate"],places=10)
        self.assertAlmostEqual(e["month_day_tmax"]["estimate"],-.004,places=6)
        self.assertAlmostEqual(e["month_night_tmax"]["estimate"],.007,places=6)

    def test_original_provider_error_source_HOLD_and_no_partial_taxa(self):
        with patch.object(m,"original_files",side_effect=OSError("PRIVATE ORIGINAL SOURCE URL")) as f:
            r=m.source(0)
            f.assert_called_once()
        self.assertEqual(r["source_status"],"SOURCE_HOLD")
        self.assertEqual(r["hold_stage"],"ORIGINAL_PROVIDER_SOURCE_UUID_MD5_BYTES")
        self.assertEqual(r["source_native_group_results"],[])
        self.assertNotIn("PRIVATE ORIGINAL SOURCE URL",str(r))

    def test_two_source_HOLD_retains_six_original_taxa(self):
        out=a.aggregate([source_hold(s,i) for i,s in enumerate(self.sources)])
        self.assertEqual(out["original_source_studies"],2)
        self.assertEqual(out["original_native_taxon_study_rows"],6)
        self.assertEqual(out["source_native_q_Tmax_sign_counts"],{"UNAVAILABLE":6})
        self.assertEqual(out["source_native_day_Tmax_sign_counts"],{"UNAVAILABLE":6})
        with self.assertRaisesRegex(ValueError,"Both original"):
            a.aggregate([source_hold(self.sources[0],0)])

if __name__=="__main__":
    unittest.main()
