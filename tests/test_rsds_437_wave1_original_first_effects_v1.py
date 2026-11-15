"""Source frame, strict HOLDS and known-truth joint covariance tests."""
from __future__ import annotations

from datetime import date
from itertools import combinations
import math
import unittest
from unittest.mock import patch

from scripts import run_rsds_437_wave1_original_gps_first_effects_v1 as run
from scripts import rsds_437_wave1_joint_source_covariance_v1 as joint

class FrozenWaveFirstOutcome(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract,cls.sources=run.freeze_original()

    def test_same_original_pre_event_source_roster(self):
        self.assertEqual(len(self.sources),16)
        self.assertEqual(sum(z["source_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
                             for z in self.sources),13)
        self.assertEqual(len({z["repository_study_id"] for z in self.sources}),16)
        self.assertEqual(sum(z["source_status"]=="SOURCE_METADATA_HOLD"
                             for z in self.sources),3)

    def test_heldout_metadata_slots_never_access_GPS(self):
        for ix in (0,1,13):
            with patch.object(run,"original_files") as original:
                r=run.source_effects(ix)
                original.assert_not_called()
            self.assertEqual(r["source_status"],"PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD")
            self.assertEqual(r["observed_original_source_event_rows_if_opened"],0)
            self.assertTrue(all(v["q"] is None for v in r["native_effects"]))
            self.assertEqual(r["source_hold_stage"],"SOURCE_ONLY_ORIGINAL_FILE_METADATA")

    def test_new_source_file_provider_failure_preserves_all_taxa(self):
        with patch.object(run,"original_files",side_effect=RuntimeError("provider unavailable")):
            r=run.source_effects(8)  # two registered marine turtles
        self.assertEqual(r["source_status"],"ORIGINAL_SOURCE_GPS_HOLD")
        self.assertEqual(len(r["native_effects"]),2)
        self.assertEqual(r["source_hold_stage"],"PROVIDER_ORIGINAL_FILE_UUID_BYTES_MD5")
        self.assertEqual(r["source_hold_reason"],"RuntimeError")
        self.assertNotIn("provider unavailable",str(r))

    def test_old_radiation_HOLD_and_source_replacement_never_qualified(self):
        with self.assertRaisesRegex(ValueError,"sixteen"):
            run.source_effects(16)
        with self.assertRaisesRegex(ValueError,"sixteen"):
            run.source_effects(-1)
        with patch.object(run,"NEGATIVE","REPLICATED"):
            r=run.source_effects(0)
        self.assertEqual(r["old_48_signed_solar_shortwave_external_NONREPLICATED_UNCHANGED"],
                         "REPLICATED")  # Detect why constant must remain pinned via contract
        self.assertEqual(self.contract["old_48_radiation_external_nonreplication"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_joint_shared_month_coupling_with_synthetic_source_not_GPS(self):
        # Five 10-day months and six private tagged animals of each native.
        # Synthetic original daily q changes with daylight, common month
        # weather shock and tag-specific noise.
        taxa=["Test alpha","Test beta","Test gamma"]
        days=[(2025,date(2025,month,day)) for month in (1,2,3,4,5)
              for day in range(2,12)]
        groups={}
        for i,name in enumerate(taxa):
            original={}
            for tag in range(6):
                values={}
                for y,d in days:
                    L=joint.daylength(43,d)
                    q=(i*.1 + (i+1)*.02*L + .08*d.month
                       + .01*math.sin(d.day*.3) + .005*(tag-2))
                    values.setdefault(y,{})[d]=(q,L)
                original[f"private_tag_{i}_{tag}"]=values
            groups[name]=original
        obj=joint.original_joint(groups,set(days),taxa,43,replicates=105,seed=7)
        self.assertEqual(obj["status"],"SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED")
        self.assertEqual(obj["eligible_common_month_blocks"],5)
        self.assertEqual(len(obj["source_date_equal_effect_vector"]),6)
        self.assertEqual(len(obj["joint_covariance_modes"]),2)
        for mode in obj["joint_covariance_modes"].values():
            V=mode["original_source_native_covariance_matrix"]
            self.assertEqual(len(V),6)
            self.assertEqual(len(mode["species_pair_contrasts"]),3)
            self.assertTrue(all(abs(V[i][j]-V[j][i])<1e-10
                                for i in range(6) for j in range(6)))
            for i,j in combinations(range(6),2):
                c=[0.0]*6;c[i]=1.0;c[j]=-1.0
                self.assertGreaterEqual(joint.quadratic(V,c),-1e-10)
        self.assertNotIn("private_tag",str(obj))
        self.assertNotIn("2025-01",str(obj))

    def test_joint_requires_every_registered_species_and_calendar_day(self):
        taxa=["SpA","SpB"]
        days={(2025,date(2025,m,d)) for m in (1,2,3,4,5)
              for d in range(2,12)}
        cohort={sp:{f"animal{a}":{2025:{date(2025,m,d):(0.5,12.0)
                for m in (1,2,3,4,5) for d in range(2,12)}}
                for a in range(6)} for sp in taxa}
        x=cohort["SpB"]["animal0"][2025]
        del x[date(2025,1,2)]
        for a in range(1,6):
            del cohort["SpB"][f"animal{a}"][2025][date(2025,1,2)]
        with self.assertRaisesRegex(ValueError,"common dates"):
            joint.private_daily(cohort,days,taxa)
        with self.assertRaisesRegex(ValueError,"registered multi-native"):
            joint.private_daily({"SpA":cohort["SpA"]},days,taxa)

if __name__=="__main__":
    unittest.main()
