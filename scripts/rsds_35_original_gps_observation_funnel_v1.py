"""Decompose source-native GPS q *measurability* into sequential gates.

Calculates only privacy-safe aggregate counts from PRIVATE, original selected
UTC-hour GPS fixes. Each source/tag original observation series is examined
without emitting animal identifiers, exact GPS dates, timestamps, coordinates,
or movement effect values. All stage cutoffs exactly reuse the previously
published RSDS native source q acceptance (0.5-2h, identical solar period,
speed <=60 km/h, >=3 DAY and 3 NIGHT steps/day, >=10 days/year).

Technical cadence step potential is separated from astronomy/daypart and
movement-speed plausibility. These are descriptive conditional denominators:
solar state and movement gates still depend on geography/biology and source
measurement. This is post-outcome method analysis, not a causal mechanism.
"""
from __future__ import annotations

from collections import defaultdict
from statistics import mean
import math

from scripts.movebank_native_ecology_20261019 import daylength, solar_period
from scripts.movebank_solar_movement_20261019 import great_circle_km, valid_interval
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon


def profile_native(native_taxon, animals):
    """Aggregate one original native source, no tag identity or GPS is returned."""
    counts={
        "source_original_tagged_animals_in_hourly_parsed_track":len(animals),
        "original_selected_UTC_hour_fix_count":0,
        "UTC_even_hour_selected_fix_count":0,
        "UTC_odd_hour_selected_fix_count":0,
        "original_adjacent_source_fix_intervals":0,
        "accepted_original_0p5_to_2h_time_intervals":0,
        "accepted_same_solar_date_and_daypart_intervals":0,
        "accepted_finite_speed_max60kmh_intervals":0,
        "accepted_native_tag_solar_dates_with_ge3_each_daypart":0,
        "accepted_native_tag_years_with_ge10_daypart_dates":0,
        "accepted_native_tags_with_ge10_daypart_dates":0,
        "accepted_native_tag_years_with_ge10_dates_and_daylength_span_ge0p5h":0,
        "accepted_native_tags_with_within_year_daylength_slope":0
    }
    for tag,fixes in animals.items():
        # Private original UTC-hour selection is inherited; source timestamp
        # within the selected hourly bin is the actual original event time.
        obs=sorted(fixes.values(),key=lambda z:z[0])
        counts["original_selected_UTC_hour_fix_count"]+=len(obs)
        counts["UTC_even_hour_selected_fix_count"]+=sum(z[0].hour%2==0 for z in obs)
        counts["UTC_odd_hour_selected_fix_count"]+=sum(z[0].hour%2==1 for z in obs)
        daypart=defaultdict(lambda:{"day":0,"night":0,"lat":[]})
        for a,b in zip(obs,obs[1:]):
            counts["original_adjacent_source_fix_intervals"]+=1
            t1,lat1,lon1=a
            t2,lat2,lon2=b
            hours=(t2-t1).total_seconds()/3600
            if not math.isfinite(hours) or not 0.5<=hours<=2.0:
                continue
            counts["accepted_original_0p5_to_2h_time_intervals"]+=1
            mid=t1+(t2-t1)/2
            lat=(lat1+lat2)/2
            lon=(lon1+lon2)/2
            start=solar_period(t1,lat1,lon1)
            end=solar_period(t2,lat2,lon2)
            center=solar_period(mid,lat,lon)
            if start!=end or start!=center:
                continue
            counts["accepted_same_solar_date_and_daypart_intervals"]+=1
            speed=great_circle_km(lat1,lon1,lat2,lon2)/hours
            if not math.isfinite(speed) or speed>60:
                continue
            original=valid_interval(t1,lat1,lon1,t2,lat2,lon2)
            if original is None or original[0]!=center[0] or original[1]!=center[1]:
                raise ValueError("Cadence funnel diverged from ORIGINAL valid_interval")
            counts["accepted_finite_speed_max60kmh_intervals"]+=1
            part=center[1]
            daypart[center[0]][part]+=1
            daypart[center[0]]["lat"].append(lat)
        year=defaultdict(list)
        for d,report in daypart.items():
            if report["day"]<3 or report["night"]<3:
                continue
            counts["accepted_native_tag_solar_dates_with_ge3_each_daypart"]+=1
            avg_lat=mean(report["lat"])
            year[d.year].append(daylength(avg_lat,d))
        supported=0
        supported_L=0
        for year_number,L in year.items():
            if len(L)<10:
                continue
            supported+=1
            counts["accepted_native_tag_years_with_ge10_daypart_dates"]+=1
            if max(L)-min(L)>=0.5:
                supported_L+=1
                counts["accepted_native_tag_years_with_ge10_dates_and_daylength_span_ge0p5h"]+=1
        counts["accepted_native_tags_with_ge10_daypart_dates"]+=int(supported>0)
        counts["accepted_native_tags_with_within_year_daylength_slope"]+=int(supported_L>0)
    if (counts["UTC_even_hour_selected_fix_count"]
        +counts["UTC_odd_hour_selected_fix_count"]
        !=counts["original_selected_UTC_hour_fix_count"]):
        raise ValueError("Original UTC-hour parity counts fail to partition source")
    keys=[
        "original_adjacent_source_fix_intervals",
        "accepted_original_0p5_to_2h_time_intervals",
        "accepted_same_solar_date_and_daypart_intervals",
        "accepted_finite_speed_max60kmh_intervals"
    ]
    funnel=[counts[k] for k in keys]
    if funnel!=sorted(funnel,reverse=True):
        raise ValueError("Original GPS stage-wise denominator is nonmonotone")
    frozen,_=one_taxon(native_taxon,animals)
    if (frozen["total_original_animals"]!=len(animals)
        or frozen["original_adjacent_hour_steps"]!=counts[keys[0]]
        or frozen["qualified_daypart_step_observations"]!=counts[keys[3]]
        or frozen["animal_local_dates_with_day_and_night"]!=
           counts["accepted_native_tag_solar_dates_with_ge3_each_daypart"]
        or frozen["admitted_animals_diel"]!=
           counts["accepted_native_tags_with_ge10_daypart_dates"]
        or frozen["admitted_animals_photoperiod_slope"]!=
           counts["accepted_native_tags_with_within_year_daylength_slope"]):
        raise ValueError("Cadence funnel does not replay original native q group support")
    c=counts
    n=c["original_adjacent_source_fix_intervals"]
    time=c["accepted_original_0p5_to_2h_time_intervals"]
    solar=c["accepted_same_solar_date_and_daypart_intervals"]
    accepted=c["accepted_finite_speed_max60kmh_intervals"]
    time_rate=time/n if n else None
    solar_given_time=solar/time if time else None
    geometry_given_solar=accepted/solar if solar else None
    end_to_end=accepted/n if n else None
    return {
        "native_taxon":native_taxon,
        "source_original_support_counts":c,
        "original_observation_gate_ratios":{
            "original_time_0p5_to_2h_over_adjacent":time_rate,
            "same_solar_date_daypart_over_time_eligible":solar_given_time,
            "plausible_speed_over_solar_eligible":geometry_given_solar,
            "total_accepted_over_adjacent":end_to_end,
            "UTC_even_hour_fraction":(
                c["UTC_even_hour_selected_fix_count"]/
                c["original_selected_UTC_hour_fix_count"]
                if c["original_selected_UTC_hour_fix_count"] else None)
        },
        "replayed_archived_native_q_support":{
            "n_tag_q":frozen["admitted_animals_diel"],
            "n_tag_daylength":frozen["admitted_animals_photoperiod_slope"],
            "q_estimable":frozen["day_night_speed_contrast"] is not None,
            "daylength_estimable":frozen["daylength_association"] is not None,
            "source_q_mean_for_private_replay_only":(
                frozen["day_night_speed_contrast"]["mean"]
                if frozen["day_night_speed_contrast"] is not None else None),
            "source_daylength_mean_for_private_replay_only":(
                frozen["daylength_association"]["mean"]
                if frozen["daylength_association"] is not None else None)
        },
        "no_animal_ID_date_location_GPS_row_or_movement_values_exported":True
    }


def redact_native_replay(report):
    """First-result q means are checked in private runner but not new outputs."""
    x={k:v for k,v in report.items() if k!="replayed_archived_native_q_support"}
    x["replayed_archived_native_q_support"]={
        k:v for k,v in report["replayed_archived_native_q_support"].items()
        if k not in ("source_q_mean_for_private_replay_only",
                     "source_daylength_mean_for_private_replay_only")
    }
    return x
