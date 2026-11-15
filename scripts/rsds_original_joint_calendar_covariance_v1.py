"""Joint study-native covariance from PRIVATE, coobserved original GPS calendar days.

This is a post-outcome, descriptive, covariance-aware extension. Two
sensitivity measures are reported:
- joint month-block bootstrap, conditional on the observed animals;
- joint month-block bootstrap with strictly positive independent Bayesian
  bootstrap animal weights by species (preserves common date coverage).
No resampled individual IDs, date/UTC, tracks, latitude or replicate vectors
are exported. The matrices capture the same *daily date-level* source
estimands, not the original equal-tag full-source q coefficients.

Important: within-calendar-month daylength is inseparable from calendar date.
The covariance estimates do NOT by themselves have nominal frequentist
coverage, and a DOI does not guarantee biological source independence.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import HERE, daylength
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import sufficient_months

SOURCE_SHA = {
    "Illinois": ("research/RSDS_ILLINOIS_THREE_NATIVE_SHARED_DATE_V1_FIRST_RESULT_20261029.json",
                 "924156e8f98d0bb9fda2f9d23aabc03060af2a31"),
    "Wyoming": ("research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json",
                "9b04b18ee73227b4e7a1715ebfb97f610af66952"),
    "Utah": ("research/RSDS_UTAH_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json",
             "87bfaf257d5c68676342118146db5141c76ae0a0"),
    "Geese": ("research/RSDS_FOUR_NEW_MULTINATIVE_GPS_SHARED_CALENDAR_V1_FIRST_RESULT_20261029.json",
              "4034d5318c914dc789a7fcd6535eedbf3a229d53"),
}
PARENT_MONTH = ("research/RSDS_FOUR_ORIGINAL_SHARED_CALENDAR_WITHIN_MONTH_DAYLENGTH_V1_FIRST_RESULT_20261108.json",
                "1071267425d23e609f8831a982f5cb334c66a2db")
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
B = 2000
SEED = 20261115

def pinned(path, sha):
    raw = (HERE/path).read_bytes()
    actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if actual != sha:
        raise ValueError("Original source result SHA changed: " + path)
    return json.loads(raw)

def frozen_parent(study):
    if study not in SOURCE_SHA:
        raise ValueError("Only four historically declared source studies")
    first = pinned(*SOURCE_SHA[study])
    combined = pinned(*PARENT_MONTH)
    if (combined["original_field_studies"] != 4
        or combined["original_study_native_pair_rows"] != 8
        or combined["previous_48_source_external_signed_radiation_NONREPLICATED_UNCHANGED"] != NEGATIVE
        or not combined["original_12_ERA5_raw_source_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Historical original shared-date studies or nulls changed")
    expected = [r for r in combined["all_four_studies_native_pair_month_results"]
                if r["original_source_field_study"] == study]
    if len(expected) != (3 if study in ("Illinois","Geese") else 1):
        raise ValueError("Historical native pair support changed")
    return first,expected

def private_daily(pool,dates,taxa):
    """Input: taxon -> tag -> year -> date -> {'original': daily-q}.

    Keys are private. Only source DOI, taxon labels, counts and full matrices
    leave the source-specific CI runner. Exclude individual identities.
    """
    if set(pool) != set(taxa) or len(set(taxa)) != len(taxa) or len(taxa)<2:
        raise ValueError("Native taxa missing, duplicated or reselected")
    dates = sorted(dates)
    if not dates or any(yr != d.year for yr,d in dates):
        raise ValueError("Missing or mismatched original source calendar")
    by_taxon = {}
    tags = {}
    for sp in taxa:
        private = pool[sp]
        if len(private)<5:
            raise ValueError("Insufficient native original tagged animals")
        label = {tag:i for i,tag in enumerate(sorted(private))}
        daily = defaultdict(list)
        for tag,years in private.items():
            for yr,ds in years.items():
                for d, rec in ds.items():
                    if (yr,d) not in dates:
                        raise ValueError("Post-selection out-of-frame daily GPS q")
                    val = rec["original"]
                    if not isinstance(val,(int,float)) or not math.isfinite(val):
                        raise ValueError("Unreliable original source q")
                    daily[(yr,d)].append((label[tag],float(val)))
        if set(daily)!=set(dates):
            raise ValueError("Not all source native taxa cover the same calendar days")
        by_taxon[sp] = daily
        tags[sp] = len(private)
    return dates,by_taxon,tags

def cov_matrix(draws):
    """Unbiased sample covariance; PSD by construction for complete joint draws."""
    if len(draws)<3:
        raise ValueError("Insufficient joint replicate draws")
    n=len(draws)
    p=len(draws[0])
    if any(len(v)!=p or not all(math.isfinite(a) for a in v) for v in draws):
        raise ValueError("Nonfinite or ragged joint vector")
    means=[mean(v[j] for v in draws) for j in range(p)]
    cov=[[sum((v[i]-means[i])*(v[j]-means[j]) for v in draws)/(n-1)
          for j in range(p)] for i in range(p)]
    return cov

def quadratic(matrix,coeff):
    if len(matrix)!=len(coeff) or any(len(row)!=len(coeff) for row in matrix):
        raise ValueError("Invalid joint covariance dimension")
    return sum(coeff[i]*matrix[i][j]*coeff[j]
               for i in range(len(coeff)) for j in range(len(coeff)))

def build_calendar_blocks(datekeys,latitude):
    """Common source/date support; thresholds identical to prior month method."""
    null_y={date:0.0 for date in datekeys}
    _, qa = sufficient_months(null_y,latitude,daylength)
    if qa["status"] != "MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE":
        raise ValueError("Frozen original source calendar has insufficient within-month support")
    raw=defaultdict(list)
    for yr,d in datekeys:
        x=daylength(latitude,d)
        if not math.isfinite(x):
            raise ValueError("Nonfinite original astronomy")
        raw[(yr,d.month)].append(((yr,d),x))
    out=[]
    for _, values in sorted(raw.items()):
        xs=[x for _,x in values]
        if len(values)<5 or max(xs)-min(xs)<.05:
            continue
        xx=mean(xs)
        sxx=sum((v-xx)**2 for v in xs)
        if sxx<=0:
            continue
        out.append({"dates":[key for key,_ in values],
                    "centered_x":{key:x-xx for key,x in values},
                    "sxx":sxx,"n":len(values)})
    if (len(out)!=qa["eligible_source_month_blocks"]
        or sum(b["n"] for b in out)!=qa["eligible_source_calendar_days"]
        or not math.isclose(sum(b["sxx"] for b in out),
                            qa["within_month_daylength_centered_sum_squares_hours2"],
                            abs_tol=1e-9,rel_tol=0)):
        raise ValueError("Changed original month-centered q/daylength estimator")
    return out, qa

def estimate(blocks,by_taxon,taxa,tag_weights,selected):
    """All taxa and both q/slope endpoints use EXACT same sampled month blocks."""
    p=len(taxa)
    qtot=[0.0]*p
    stot=[0.0]*p
    n=0
    xx=0.0
    for ib in selected:
        block=blocks[ib]
        xx+=block["sxx"]
        n+=block["n"]
        for key in block["dates"]:
            centered=block["centered_x"][key]
            for i,taxon in enumerate(taxa):
                samples=by_taxon[taxon][key]
                weights=tag_weights[taxon]
                den=sum(weights[tag] for tag,_ in samples)
                if den<=0:
                    raise ValueError("Uncovered original common date in joint reweighting")
                value=sum(weights[tag]*v for tag,v in samples)/den
                qtot[i]+=value
                stot[i]+=centered*value
    if n<=0 or xx<=0:
        raise ValueError("Insufficient source monthly daylength slope")
    return [a/n for a in qtot] + [b/xx for b in stot]

def comparison(point, matrix, taxa, expected):
    p=len(taxa)
    differences=[]
    for pair in expected:
        a,b=pair["original_native_pair"]
        i,j=taxa.index(a),taxa.index(b)
        coef=[0.0]*(2*p)
        coef[i]=1.0;coef[j]=-1.0
        dc=point[i]-point[j]
        coef_slope=[0.0]*(2*p)
        coef_slope[i+p]=1.0;coef_slope[j+p]=-1.0
        ds=point[i+p]-point[j+p]
        if not math.isclose(ds,pair["month_slope"],rel_tol=0,abs_tol=2e-10):
            raise ValueError("Recomputed native pair month slope disagrees with immutable source result: "+a+" vs "+b)
        v1=quadratic(matrix,coef)
        v2=quadratic(matrix,coef_slope)
        if min(v1,v2)<-1e-12:
            raise ValueError("Negative joint contrast variance")
        differences.append({
            "first_taxon":a,"second_taxon":b,
            "source_date_equal_q_difference":dc,
            "source_month_slope_difference_per_hour":ds,
            "source_month_slope_prior_replayed":True,
            "joint_bootstrap_q_difference_se":math.sqrt(max(0,v1)),
            "joint_bootstrap_month_slope_difference_se":math.sqrt(max(0,v2))
        })
    return differences

def summarize_joint(study,doi,pool,dates,taxa,latitude,first,expected,
                    nboot=B,seed=SEED):
    """Return public aggregate-only source-native covariance and exact checks."""
    if nboot<100:
        raise ValueError("Insufficient covariance replicates")
    dates,within,tags=private_daily(pool,dates,taxa)
    blocks,qa=build_calendar_blocks(dates,latitude)
    indices=list(range(len(blocks)))
    same={sp:[1.0]*tags[sp] for sp in taxa}
    point=estimate(blocks,within,taxa,same,indices)
    dim=2*len(taxa)
    modes={}
    for offset,mode in enumerate(("paired_month_only","paired_month_plus_animal_weights")):
        rng=random.Random(seed+offset*100003)
        draws=[]
        for _ in range(nboot):
            selected=[rng.randrange(len(blocks)) for _ in blocks]
            if mode=="paired_month_only":
                weights=same
            else:
                weights={sp:[rng.expovariate(1.0) for j in range(tags[sp])]
                         for sp in taxa}
            draws.append(estimate(blocks,within,taxa,weights,selected))
        v=cov_matrix(draws)
        for i in range(dim):
            if v[i][i]<0 or not math.isfinite(v[i][i]):
                raise ValueError("Nonpositive empirical variance")
            for j in range(dim):
                if not math.isclose(v[i][j],v[j][i],abs_tol=1e-12,rel_tol=0):
                    raise ValueError("Non-symmetric covariance")
        for coeff in ([1.0]*dim,
                      [1.0 if z%2==0 else -1.0 for z in range(dim)]):
            if quadratic(v,coeff)<-1e-10:
                raise ValueError("Not a valid PSD joint covariance")
        means=[mean(x[i] for x in draws) for i in range(dim)]
        modes[mode]={
            "replicates":nboot,
            "drawn_parameter_order":[("q_date_equal:" if i<len(taxa)
                                     else "q_within_month_slope_per_hour:")+
                                      taxa[i%len(taxa)] for i in range(dim)],
            "joint_sampling_covariance_matrix":v,
            "replicate_means":means,
            "species_contrasts":comparison(point,v,taxa,expected),
        }
    return {
        "version":"RSDS_ORIGINAL_GPS_JOINT_MONTH_ANIMAL_COVARIANCE_V1",
        "research_stage":"POST_OUTCOME_EXPLORATORY_COVARIANCE_NOT_NEW_FIELD_REPLICATION",
        "study":study,"original_DOI":doi,
        "status":"ORIGINAL_SOURCE_JOINT_COVARIANCE_ESTIMATED",
        "original_source_result_SHA":SOURCE_SHA[study][1],
        "original_four_month_result_SHA":PARENT_MONTH[1],
        "source_species_in_fixed_order":list(taxa),
        "source_original_tag_counts_by_taxon":tags,
        "original_source_shared_calendar_days":len(dates),
        "eligible_common_month_days":qa["eligible_source_calendar_days"],
        "eligible_common_year_month_blocks":qa["eligible_source_month_blocks"],
        "source_month_centered_sxx_hours2":qa["within_month_daylength_centered_sum_squares_hours2"],
        "date_equal_estimates":{
            sp:{"q_date_equal":point[i],
                "q_within_month_slope_per_hour":point[i+len(taxa)]}
            for i,sp in enumerate(taxa)
        },
        "covariance_sensitivity":modes,
        "historical_external_shortwave_nonreplication_UNCHANGED":NEGATIVE,
        "historic_12_ERA5_raw_JSON_HOLD_UNCHANGED":True,
        "no_individual_date_location_or_weight_draws_exported":True,
        "limitations":[
            "Post-outcome study/species design and no independent field replication",
            "Raw GPS q is observed straight-line step displacement, not true activity time",
            "Joint month-block covariance conditional on observed animal set",
            "Animal-weight sensitivity is Bayesian bootstrap-style, not calibrated frequentist CI",
            "Only common-date and eligible-month source support, not original full-year equal-tag q",
            "Same source calendar date does not prove sympatry or equal solar exposure",
            "No causal attribution to daylength or between-species interaction",
            "Four source DOI are purposely selected and no pooled global meta effect authorized"
        ]
    }
