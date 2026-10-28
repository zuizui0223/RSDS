"""Three-way exact original tag+calendar year+date UTC GPS cadence known-truth."""
from __future__ import annotations

import json
import unittest
from datetime import date,timedelta

from scripts.run_rsds_gps_three_phase_exact_dates_v2 import (
    AXES, NEGATIVE, frozen, triple_from_private_original_source,
    aggregate_private_triple
)


def three(n=6,days=24):
    full,even,odd={},{},{}
    start=date(2020,4,1)
    for j in range(n):
        orig={};ev={};od={}
        for i in range(days):
            d=start+timedelta(days=i)
            q=0.02*j+0.005*i
            x=10.0+.06*i
            orig[d]=(q,x)
            ev[d]=(q+.03+j*.001,x+.0002)
            od[d]=(q-.01-j*.0004,x-.0001)
        tag="PRIVATE-tag-"+str(j)
        full[tag]={2020:orig}
        even[tag]={2020:ev}
        odd[tag]={2020:od}
    return full,even,odd


class TripleOriginalUTCCadenceKnownTruth(unittest.TestCase):
    def test_source_first_original_eight_studies_phase_alias_SHAs_frozen(self):
        c,first=frozen()
        self.assertEqual(c["version"],"RSDS_EXACT_TRIPLE_TAG_DATE_UTC_PHASE_GPS_V2")
        self.assertEqual(first["original_field_studies"],8)
        self.assertEqual(first["original_source_taxon_rows"],12)
        self.assertTrue(first["historic_ERA5_old_12_source_raw_SHA_HOLDS_UNCHANGED"])
        self.assertEqual(first["historical_48_source_negative_radiation_replication_UNCHANGED"],
                         NEGATIVE)

    def test_actual_three_schedule_identical_animal_year_day_set_q_identity(self):
        full,even,odd=three()
        result,coverage=triple_from_private_original_source(full,even,odd)
        self.assertEqual(coverage["source_original_tags_with_at_least_ten_triple_days"],6)
        self.assertEqual(coverage["triply_qualified_tag_years_with_at_least_ten_identical_dates"],6)
        self.assertEqual(coverage["triply_qualified_original_animal_days"],144)
        a=aggregate_private_triple(result)
        self.assertEqual(set(a),set(AXES))
        self.assertAlmostEqual(a["q_even_minus_original"]["group_source_mean_delta"],.0325,places=12)
        self.assertAlmostEqual(a["q_odd_minus_original"]["group_source_mean_delta"],-.011,places=12)
        self.assertAlmostEqual(a["q_even_minus_odd"]["group_source_mean_delta"],.0435,places=12)
        self.assertEqual(a["q_even_minus_odd"]["status"],"NOMINAL_POSITIVE")
        self.assertEqual(a["q_odd_minus_original"]["status"],"NOMINAL_NEGATIVE")
        self.assertAlmostEqual(a["daylength_even_minus_odd"]["group_source_mean_delta"],0,places=12)
        self.assertNotIn("PRIVATE",json.dumps(a))
        self.assertEqual(a,aggregate_private_triple(result))

    def test_one_parity_not_on_same_dates_means_no_artificial_triple_match(self):
        full,even,odd=three(n=6)
        for tag in list(odd):
            odd[tag][2020]=dict(list(odd[tag][2020].items())[:8])
        result,coverage=triple_from_private_original_source(full,even,odd)
        self.assertEqual(coverage["source_original_tags_with_at_least_ten_triple_days"],0)
        a=aggregate_private_triple(result)
        self.assertEqual(a["q_even_minus_odd"]["status"],"UNAVAILABLE")
        self.assertIsNone(a["q_even_minus_odd"]["group_source_mean_delta"])
        self.assertIsNone(a["q_even_minus_odd"]["descriptive_source_tag_bootstrap95"])

    def test_one_original_tag_not_available_in_odd_or_even_not_pseudoreplicated(self):
        full,even,odd=three(n=6)
        del odd["PRIVATE-tag-0"]
        del even["PRIVATE-tag-1"]
        result,coverage=triple_from_private_original_source(full,even,odd)
        self.assertEqual(coverage["source_original_tags_with_at_least_ten_triple_days"],4)
        a=aggregate_private_triple(result)
        self.assertEqual(a["q_even_minus_odd"]["source_original_same_tags_and_triple_dates"],4)
        self.assertEqual(a["q_even_minus_odd"]["status"],"UNAVAILABLE")
        self.assertIsNone(a["q_even_minus_odd"]["group_source_mean_delta"])

    def test_common_original_daylength_x_reused_across_three_schedules(self):
        full,even,odd=three(n=6)
        for tag in even:
            even[tag][2020]={d:(q,L+10) for d,(q,L) in even[tag][2020].items()}
        for tag in odd:
            odd[tag][2020]={d:(q,L-10) for d,(q,L) in odd[tag][2020].items()}
        out,qa=triple_from_private_original_source(full,even,odd)
        summaries=aggregate_private_triple(out)
        self.assertAlmostEqual(summaries["daylength_even_minus_original"]["group_source_mean_delta"],0,places=12)
        self.assertAlmostEqual(summaries["daylength_odd_minus_original"]["group_source_mean_delta"],0,places=12)


if __name__=="__main__":
    unittest.main()
