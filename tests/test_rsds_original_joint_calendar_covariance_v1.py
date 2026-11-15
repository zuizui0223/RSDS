"""Known-truth joint covariance checks without original telemetry or IDs."""
from __future__ import annotations

import copy
import datetime as dt
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import rsds_original_joint_calendar_covariance_v1 as m
from scripts import run_rsds_original_joint_covariance_v1 as runner

TAXA=("Sp.alpha","Sp.beta","Sp.gamma")


def synthetic():
    dates=[(2025,dt.date(2025,month,day)) for month in (1,2,3,4,5)
           for day in range(2,12)]
    groups={}
    for s,sp in enumerate(TAXA):
        tags={}
        for t in range(6):
            years={2025:{}}
            for yr,d in dates:
                x=m.daylength(43.0,d)
                # Shared month shock and common within-month fluctuations.
                monthshock=(d.month-3)*.09
                common=.02*math.sin(d.day/2.3)
                q=(s*.4 + (s+1)*.065*x + monthshock+common+
                   .007*(t-2.5)+.005*(s+1)*math.cos(d.day/3))
                years[yr][d]={"original":q}
            tags[f"private_{s}_{t}"]=years
        groups[sp]=tags
    return groups,set(dates)


class PairedMonthOriginalCovariance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool,cls.dates=synthetic()
        dates,daily,tags=m.private_daily(cls.pool,cls.dates,TAXA)
        cls.blocks,cls.qa=m.build_calendar_blocks(dates,43)
        cls.point=m.estimate(cls.blocks,daily,TAXA,
                             {sp:[1.0]*tags[sp] for sp in TAXA},
                             list(range(len(cls.blocks))))
        cls.dates_stored=dates
        cls.daily=daily
        cls.tags=tags

    def test_source_parent_pins_exist_and_histories_fail_closed(self):
        for source,doi in runner.EXPECTED.items():
            prev,comparisons=m.frozen_parent(source)
            self.assertEqual(len(comparisons),3 if source in ("Illinois","Geese") else 1)
            self.assertTrue(doi["days"]>=120)
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(m,"HERE",Path(tmp)):
                with self.assertRaises(FileNotFoundError):
                    m.frozen_parent("Illinois")

    def test_common_month_thresholds(self):
        self.assertEqual(len(self.dates_stored),50)
        self.assertEqual(self.qa["eligible_source_month_blocks"],5)
        self.assertEqual(self.qa["eligible_source_calendar_days"],50)
        self.assertGreater(self.qa["within_month_daylength_centered_sum_squares_hours2"],.25)

    def test_shared_date_coverage_fail_closed(self):
        no=copy.deepcopy(self.pool)
        last=sorted(self.dates)[0]
        for tag in no[TAXA[0]]:
            del no[TAXA[0]][tag][last[0]][last[1]]
        with self.assertRaisesRegex(ValueError,"same calendar days"):
            m.private_daily(no,self.dates,TAXA)
        no=copy.deepcopy(self.pool)
        no[TAXA[0]].pop(next(iter(no[TAXA[0]])))
        no[TAXA[0]].pop(next(iter(no[TAXA[0]])))
        with self.assertRaisesRegex(ValueError,"Insufficient native"):
            m.private_daily(no,self.dates,TAXA)

    def test_shared_q_and_slope_species_contrasts_algebra(self):
        # Construct sealed-like expected coefficients from synthetic joint point.
        expected=[]
        for a,b in ((TAXA[0],TAXA[1]),(TAXA[1],TAXA[2]),
                    (TAXA[0],TAXA[2])):
            i,j=TAXA.index(a),TAXA.index(b)
            expected.append({"original_native_pair":[a,b],
                             "month_slope":self.point[3+i]-self.point[3+j]})
        output=m.summarize_joint("Synthetic","source-known-truth",self.pool,
                                 self.dates,TAXA,43.0,{},expected,nboot=140,seed=3)
        self.assertEqual(output["eligible_common_month_days"],50)
        self.assertEqual(output["eligible_common_year_month_blocks"],5)
        for typ in ("paired_month_only","paired_month_plus_animal_weights"):
            model=output["covariance_sensitivity"][typ]
            self.assertEqual(len(model["drawn_parameter_order"]),6)
            matrix=model["joint_sampling_covariance_matrix"]
            self.assertEqual(len(matrix),6)
            for i in range(6):
                self.assertGreaterEqual(matrix[i][i],0)
                for j in range(6):
                    self.assertAlmostEqual(matrix[i][j],matrix[j][i],places=12)
            # Common calendar month innovations induce cross-taxon dependence.
            self.assertGreater(abs(matrix[0][1]),1e-8)
            for coeff in ([1,1,1,1,1,1],[2,-1,0,0,4,-3],
                          [0,0,0,1,-1,0]):
                self.assertGreaterEqual(m.quadratic(matrix,coeff),-1e-10)
            for item in model["species_contrasts"]:
                self.assertTrue(item["source_month_slope_prior_replayed"])
                self.assertGreaterEqual(item["joint_bootstrap_q_difference_se"],0)
                self.assertGreaterEqual(item["joint_bootstrap_month_slope_difference_se"],0)

    def test_calendar_month_resampling_and_positive_animal_weights_are_distinct(self):
        expected=[{"original_native_pair":[TAXA[0],TAXA[1]],
                  "month_slope":self.point[3]-self.point[4]}]
        a=m.summarize_joint("Synthetic","source-known-truth",self.pool,
                            self.dates,TAXA,43,{},expected,nboot=110,seed=17)
        b=m.summarize_joint("Synthetic","source-known-truth",self.pool,
                            self.dates,TAXA,43,{},expected,nboot=110,seed=17)
        self.assertEqual(a,b)
        block=a["covariance_sensitivity"]["paired_month_only"][
            "joint_sampling_covariance_matrix"]
        two=a["covariance_sensitivity"]["paired_month_plus_animal_weights"][
            "joint_sampling_covariance_matrix"]
        self.assertNotEqual(block,two)
        for key in ("original_DOI","date_equal_estimates","source_species_in_fixed_order"):
            self.assertIn(key,a)
        private=str(a)
        self.assertNotIn("private_0_0",private)
        self.assertNotIn("2025-01-02",private)

    def test_joint_false_independence_rejected_and_old_biological_gates_invariant(self):
        no=copy.deepcopy(self.pool)
        no[TAXA[1]][next(iter(no[TAXA[1]]))][2025][dt.date(2025,1,2)]["original"]=float("nan")
        with self.assertRaisesRegex(ValueError,"Unreliable original"):
            m.private_daily(no,self.dates,TAXA)
        with self.assertRaisesRegex(ValueError,"Only four"):
            runner.compute("new_selected_favorable_study")


if __name__=="__main__":
    unittest.main()
