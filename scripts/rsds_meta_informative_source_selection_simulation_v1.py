"""Known-truth source-measurability and meta-sign selection simulation.

A METHOD-ONLY, source-level DGP; no original animal GPS is opened.
Finite source frame is exactly 109 single-native + 11 multi-native, and
the counterfactual selected sample uses fixed 8+8 then 24+11 study quotas.
Wave2 source quota in REAL data was adaptive after Wave1, so this
simulation is intentionally *not* a faithful prospective replication.

The true source q distribution is +0.04, and we vary ONLY whether
source-complete observation R also depends on the (unobserved) q.
A true-prensity Horvitz–Thompson estimator is an oracle benchmark,
NOT executable on unmeasurable source q in real ecology. The cadence-
only benchmark knows the synthetic DGP marginal propensity exactly
but cannot identify non-ignorable selection on the missing q.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from statistics import NormalDist, mean, stdev

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

SPEC=HERE/"research/RSDS_META_OBSERVATION_SELECTION_KNOWN_TRUTH_SIM_V1_CONTRACT.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
GROUPS=("single_native","multi_native")
Z_QUANTILES=tuple(NormalDist().inv_cdf((i+.5)/41) for i in range(41))


def freeze():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    if (c["version"]!="RSDS_META_OBSERVATION_SELECTION_KNOWN_TRUTH_SIMULATION_V1"
        or c["original_source_frame"]!={"single_native":109,"multi_native":11}
        or c["hypothetical_wave1_sample"]!={"single_native":8,"multi_native":8}
        or c["hypothetical_combined_sample"]!={"single_native":24,"multi_native":11}
        or c["repetitions_per_scenario"]!=1000
        or c["random_seed"]!=20261115
        or c["true_mean_log1p_night_minus_day_speed_q"]!=.04
        or c["true_between_source_standard_deviation"]!=.12
        or c["logit_source_observation_intercept"]!={"single_native":-.7,"multi_native":-1.55}
        or [x["beta_q"] for x in c["source_q_response_dependence_scaled_betas"]]!=[0,-1.5,1.5]
        or c["historical_48_shortwave_external"]!=NEGATIVE
        or c["historical_12_ERA5_raw_source_HOLD"] is not True):
        raise ValueError("Known-truth simulation DGP not frozen")
    d=pinned(c["parent_source_35_meta_path"],c["parent_source_35_git_blob"])
    f=pinned(c["parent_observation_result_path"],c["parent_observation_result_git_blob"])
    if (d["original_source_study_units"]!=35
        or d["original_metadata_candidate_frame"]!=120
        or d["finite_120_metadata_frame_source_measurability"]["q"]["source_complete_within_sample_count"]!=10
        or f["source_original_registered_native_rows"]!=50
        or f["original_source_funnel_replayed"]!=25
        or f["original_human_GPS_source_opened"] is not False
        or f["original_48_solar_negative_external_NONREPLICATION_UNCHANGED"]!=NEGATIVE
        or not f["old_12_ERA5_raw_source_versions_HOLD_UNCHANGED"]):
        raise ValueError("Original observed source sample or historical HOLD changed")
    return c


def logistic(x):
    if x>=0:
        return 1/(1+math.exp(-x))
    z=math.exp(x)
    return z/(1+z)


def probability(c,group,cadence,q,beta):
    """Synthetic Bernoulli propensity; source level, not GPS event level."""
    alpha=c["logit_source_observation_intercept"][group]
    eta=(alpha+c["logit_cadence_beta"]*
         (cadence-c["cadence_center"])/c["cadence_divisor"]+
         beta*(q-c["true_mean_log1p_night_minus_day_speed_q"])/
         c["true_between_source_standard_deviation"])
    return logistic(eta)


def marginal_probability(c,group,cadence,beta):
    """Oracle known-DGP P(R | cadence, stratum), integrating out q."""
    alpha=c["logit_source_observation_intercept"][group]
    eta=(alpha+c["logit_cadence_beta"]*
         (cadence-c["cadence_center"])/c["cadence_divisor"])
    return mean(logistic(eta+beta*z) for z in Z_QUANTILES)


def generate_world(c,world):
    """Truth, cadence, response uniforms and source order independent of scenario."""
    rng=random.Random(c["random_seed"]+100003*world)
    frame={}
    for group in GROUPS:
        n=c["original_source_frame"][group]
        low,high=c["cadence_metadata_uniform_range"]
        values=[]
        for source_index in range(n):
            q=(c["true_mean_log1p_night_minus_day_speed_q"]+
               c["true_between_source_standard_deviation"]*rng.gauss(0,1))
            cadence=low+(high-low)*rng.random()
            values.append({"q":q,"cadence":cadence,"R_uniform":rng.random()})
        rng.shuffle(values)
        frame[group]=values
    return frame


def target(c,frame):
    return sum(sum(v["q"] for v in frame[g]) for g in GROUPS)/sum(
        c["original_source_frame"].values())


def sample(c,frame,wave):
    key="hypothetical_wave1_sample" if wave==1 else "hypothetical_combined_sample"
    return {g:frame[g][:c[key][g]] for g in GROUPS}


def observed(c,sampled,beta):
    result={}
    for g,units in sampled.items():
        result[g]=[]
        for v in units:
            p=probability(c,g,v["cadence"],v["q"],beta)
            r=v["R_uniform"]<p
            marginal=marginal_probability(c,g,v["cadence"],beta)
            result[g].append((v["q"],r,p,marginal))
    return result


def poststratified(c,groups,indicator):
    """Respondent-only source means, undefined if any stratum has no response."""
    weighted=0.0
    for g in GROUPS:
        q=[x[0] for x in groups[g] if indicator(x)]
        if not q:
            return None
        weighted+=c["original_source_frame"][g]*mean(q)
    return weighted/sum(c["original_source_frame"].values())


def HT(c,groups,weight_field):
    """Unbiased source-frame Horvitz-Thompson under known selection probability."""
    if weight_field not in (2,3):
        raise ValueError("True conditional or cadence-only oracle propensity")
    total=0.0
    for g in GROUPS:
        vals=groups[g]
        total+=c["original_source_frame"][g]/len(vals)*sum(
            q/p for q,r,*rest in vals for p in [rest[weight_field-2]]
            if r
        )
    return total/sum(c["original_source_frame"].values())


def world_scenario(c,frame,beta):
    complete=sample(c,frame,2)
    early=sample(c,frame,1)
    obs=observed(c,complete,beta)
    w1_obs={g:obs[g][:c["hypothetical_wave1_sample"][g]] for g in GROUPS}
    w1=poststratified(c,w1_obs,lambda x:x[1])
    w2=poststratified(c,obs,lambda x:x[1])
    full=poststratified(c,{g:[(x["q"],True,1.0,1.0) for x in complete[g]]
                            for g in GROUPS},lambda x:True)
    return {
      "truth":target(c,frame),
      "all_sampled_q_observed_oracle":full,
      "respondent_only_wave1":w1,
      "respondent_only_35":w2,
      "true_source_q_response_probability_oracle_HT":HT(c,obs,2),
      "cadence_only_marginal_response_probability_oracle_HT":HT(c,obs,3),
      "respondent_counts":{g:sum(x[1] for x in obs[g]) for g in GROUPS}
    }


def metric(name,worlds):
    obs=[x for x in worlds if x[name] is not None]
    err=[x[name]-x["truth"] for x in obs]
    n=len(err)
    if n<2:raise ValueError("Insufficient simulated evaluations")
    bias=mean(err)
    se=stdev(err)/math.sqrt(n)
    mismatch=sum((x[name]<0) != (x["truth"]<0) for x in obs)
    return {
       "worlds_with_defined_estimator":n,
       "worlds_missing_estimator":len(worlds)-n,
       "bias_vs_actual_finite_120_source_truth":bias,
       "monte_carlo_SE_of_bias_over_independent_worlds":se,
       "monte_carlo_95_interval_for_mean_bias_only":[bias-1.96*se,bias+1.96*se],
       "RMSE_vs_finite_frame_truth":math.sqrt(mean(z*z for z in err)),
       "wrong_signed_source_mean_fraction_among_estimator_defined":mismatch/n,
       "not_a_biological_or_model_coverage_CI":True
    }


def analyze(c):
    n=c["repetitions_per_scenario"]
    worlds=[generate_world(c,i) for i in range(n)]
    out={}
    for sc in c["source_q_response_dependence_scaled_betas"]:
        beta=sc["beta_q"]
        simulated=[world_scenario(c,f,beta) for f in worlds]
        both=[r for r in simulated if (r["respondent_only_wave1"] is not None
                                      and r["respondent_only_35"] is not None)]
        flip=sum(r["respondent_only_wave1"]<0<r["respondent_only_35"]
                 for r in both)
        total_selected={g:c["hypothetical_combined_sample"][g]*n for g in GROUPS}
        response_rate={g:sum(r["respondent_counts"][g] for r in simulated)/
                       total_selected[g] for g in GROUPS}
        estimators={}
        for label in (
            "all_sampled_q_observed_oracle",
            "respondent_only_35",
            "true_source_q_response_probability_oracle_HT",
            "cadence_only_marginal_response_probability_oracle_HT"
        ):
            estimators[label]=metric(label,simulated)
        out[sc["label"]]={
           "beta_true_source_q_in_measurability_logit":beta,
           "synthetic_source_q_frame_truth_mean_across_worlds":
               mean(x["truth"] for x in simulated),
           "synthetic_true_source_q_positive_truth_fraction":
               sum(x["truth"]>0 for x in simulated)/n,
           "actual_source_response_fraction_by_stratum":response_rate,
           "two_wave_respondent_only_estimates_both_available_worlds":len(both),
           "wave1_negative_to_wave2_positive_sign_flip_fraction_given_both_available":
               flip/len(both) if both else None,
           "estimators":estimators,
           "scientific_scope":"KNOWN_TRUTH_SYNTHETIC_SOURCE_SELECTION_ONLY_NO_ORIGINAL_GPS"
        }
    return {
       "version":"RSDS_META_INFORMATIVE_SOURCE_MEASURABILITY_KNOWN_TRUTH_V1",
       "stage":"SYNTHETIC_SEPARATE_FROM_ARCHIVED_35_SOURCE_GPS_ECOLOGY",
       "source_metadata_frame":{"single_native":109,"multi_native":11},
       "counterfactual_predeclared_sample":{"single_native":24,"multi_native":11},
       "worlds_per_scenario":n,
       "source_truth_and_cadence_shared_across_all_three_scenarios":True,
       "source_response_condition_depends_on_missing_q_in_MNAR":True,
       "synthetic_scenarios":out,
       "source_file_MD5_ORIGINAL_35_GPS_NOT_ACCESSED_BY_THIS_SIMULATION":True,
       "old_48_source_shortwave_direction_external_NOT_REPLICATED":True,
       "old_12_raw_ERA5_source_HOLD_UNCHANGED":True,
       "claims_not_supported":[
           "Simulation parameters are NOT estimated from the 35 original study sources",
           "Simulated source observation R is a single summary proxy, not mechanistic GPS temporal cadence",
           "Positive/negative MNAR dependence is a possible failure mechanism, not evidence of real causal MNAR",
           "Even oracle propensity weights may have high variance with a small multi-native source stratum",
           "Neither a synthetic pooled sign nor modeled source admission generalizes to global wildlife"
       ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=analyze(freeze())
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    report={name:{
        "response":v["actual_source_response_fraction_by_stratum"],
        "bias":{
            label:round(d["bias_vs_actual_finite_120_source_truth"],6)
            for label,d in v["estimators"].items()},
        "wave1_to_wave2_negative_to_positive_flip_fraction":
            v["wave1_negative_to_wave2_positive_sign_flip_fraction_given_both_available"]
    } for name,v in result["synthetic_scenarios"].items()}
    print("RSDS_MEASURABILITY_KNOWN_TRUTH",json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
