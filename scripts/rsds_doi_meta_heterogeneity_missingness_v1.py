"""Post-outcome study-level heterogeneity and missingness audit for RSDS.

Reads ONLY the checksum-pinned original twelve-study summary. It does not open
GPS, reselect animals/taxa, alter original endpoints, or claim a global survey.
First/second-stage sampling and study-level ecological selection are not random.
All CI/PI outputs are exploratory *model sensitivity*, not nominal coverage
guarantees: tagged-animal percentile bootstrap CIs are only SE proxies, and
within-field-study covariance is unavailable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, stdev

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json"
SOURCE_GIT_BLOB_SHA = "af0574c6bb6383402a97b5bebadb8cf0408e9740"
EXPECTED = {"studies": 12, "rows": 22, "estimable_studies": 9,
            "estimable_rows": 15}
T975_DF8 = 2.306004135204166  # Frozen nine-DOI analysis, two-sided 95%.


def read_pinned():
    raw = SOURCE.read_bytes()
    sha = hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()
    if sha != SOURCE_GIT_BLOB_SHA:
        raise ValueError("Frozen original 12-study Git blob changed; stop")
    return json.loads(raw)


def observed_rows(document):
    rows = document["all_original_study_native_taxon_rows_including_HOLD"]
    if len(rows) != EXPECTED["rows"]:
        raise ValueError("Original 22 native study-taxon rows changed")
    identities = [(r["original_study_doi"], r["native_taxon"]) for r in rows]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate native study-taxon rows")
    if len(set(x[0] for x in identities)) != EXPECTED["studies"]:
        raise ValueError("Original 12 DOI frame changed")
    if document["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"] != (
        "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
    ) or document["original_12_ERA5_raw_JSON_source_HOLDS_UNCHANGED"] is not True:
        raise ValueError("Historical negative replication or weather HOLD was changed")
    for field in ("q", "daylength"):
        available = 0
        for row in rows:
            e = row[field]
            if e["mean"] is None:
                if e["status"] != "UNAVAILABLE" or e["tag_bootstrap95"] is not None:
                    raise ValueError("HOLD must remain unavailable, never zero")
            else:
                lo, hi = e["tag_bootstrap95"]
                if not (math.isfinite(e["mean"]) and lo <= e["mean"] <= hi
                        and lo < hi and e["source_tags"] >= 5):
                    raise ValueError("Invalid original individual bootstrap summary")
                available += 1
        if available != EXPECTED["estimable_rows"]:
            raise ValueError("15 supported native source rows changed")
    return rows


def grouped(rows):
    d = {}
    for r in rows:
        d.setdefault(r["original_study_doi"], []).append(r)
    return {doi: d[doi] for doi in sorted(d)}


def study_effects(rows, field, correlation):
    """Source-balanced means; unknown within-source correlations: rho = 0 or 1.

    Each tagged-animal percentile bootstrap CI width / 3.92 is merely a
    normal-approximation SE proxy, NOT an SE estimated from joint tag bootstrap
    draws. Equicorrelation rho=1 is the upper-variance extreme.
    """
    if correlation not in (0.0, 1.0):
        raise ValueError("Only two prespecified covariance endpoints allowed")
    result = []
    for doi, cohort in grouped(rows).items():
        available = [r[field] for r in cohort if r[field]["mean"] is not None]
        datum = {"doi": doi, "registered_native_taxa": len(cohort),
                 "estimable_native_taxa": len(available)}
        if available:
            y = mean(e["mean"] for e in available)
            se = [(e["tag_bootstrap95"][1] - e["tag_bootstrap95"][0]) / 3.92
                  for e in available]
            ss = sum(z*z for z in se)
            variance = (ss + correlation*(sum(se)**2 - ss)) / len(se)**2
            datum.update({"estimate": y, "se_proxy": math.sqrt(variance),
                          "status": "SUPPORTED"})
        else:
            datum.update({"estimate": None, "se_proxy": None,
                          "status": "HOLD"})
        result.append(datum)
    return result


def random_effects(studies):
    """Paule-Mandel heterogeneity and modified Hartung-Knapp t8 intervals.

    This is a *diagnostic*, not inferential calibration for the purposively
    selected studies, correlated source taxa, or unknown nonresponse mechanism.
    """
    values = [s for s in studies if s["status"] == "SUPPORTED"]
    n = len(values)
    if n != EXPECTED["estimable_studies"]:
        raise ValueError("Random effects exploratory route requires nine DOI")
    y = [s["estimate"] for s in values]
    v = [s["se_proxy"]**2 for s in values]

    def weighted(tau2):
        w = [1.0/(z+tau2) for z in v]
        total = sum(w)
        pooled = sum(wi*yi for wi, yi in zip(w,y))/total
        q = sum(wi*(yi-pooled)**2 for wi, yi in zip(w,y))
        return pooled, q, total

    q0 = weighted(0.0)[1]
    if q0 <= n-1:
        tau2 = 0.0
    else:
        lo, hi = 0.0, 1e-8
        while weighted(hi)[1] > n-1:
            hi *= 2.0
        for _ in range(90):
            mid = (lo+hi)/2
            if weighted(mid)[1] > n-1:
                lo = mid
            else:
                hi = mid
        tau2 = (lo+hi)/2
    pooled, q_tau, sw = weighted(tau2)
    hk = max(1.0, q_tau/(n-1))  # avoid HK variance below RE model variance
    semu = math.sqrt(hk/sw)
    ci_half = T975_DF8*semu
    pi_half = T975_DF8*math.sqrt(tau2+semu**2)
    return {
        "method": "Paule-Mandel tau2; modified Hartung-Knapp t(df=8)",
        "studies": n, "study_equal_mean": mean(y),
        "observed_study_sd": stdev(y), "random_effects_mean": pooled,
        "approx_se_of_mean": semu, "tau2": tau2, "tau": math.sqrt(tau2),
        "Q0_se_proxy": q0,
        "I2_proxy": max(0.0, (q0-(n-1))/q0) if q0 > 0 else 0.0,
        "model_95_interval": [pooled-ci_half, pooled+ci_half],
        "model_95_prediction_interval": [pooled-pi_half, pooled+pi_half],
        "assumption": "NOT world-representative, causal, or calibrated 95%-coverage"
    }


def missingness(rows, field, studies):
    """Two explicit, conditional sensitivity questions; no zero-imputation."""
    group = grouped(rows)
    available = [s["estimate"] for s in studies if s["estimate"] is not None]
    nmissing = len(group)-len(available)
    if nmissing != 3:
        raise ValueError("Missing three whole-source DOI is part of frozen audit")
    observed_native = [r[field]["mean"] for r in rows
                       if r[field]["mean"] is not None]
    native_low, native_high = min(observed_native), max(observed_native)

    def hypothetical_full_frame(fill):
        return mean(mean(fill if x[field]["mean"] is None else x[field]["mean"]
                         for x in group[doi]) for doi in sorted(group))

    return {
        "unsupported_native_rows": sum(r[field]["mean"] is None for r in rows),
        "fully_unsupported_source_DOI": nmissing,
        "missing_source_equal_common_effect_at_zero_pooled_mean":
            -sum(available)/nmissing,
        "all_registered_12_studies_hypothetical_native_range": {
            "observed_native_min": native_low,
            "observed_native_max": native_high,
            "all_missing_native_at_observed_min": hypothetical_full_frame(native_low),
            "all_missing_native_at_observed_max": hypothetical_full_frame(native_high),
        },
        "interpretation": "Illustrative pattern-mixture/tipping diagnostic, NOT "
                          "identified missing-effect bounds; unobserved effects "
                          "need not be within observed native range."
    }


def analyze(document):
    rows = observed_rows(document)
    out = {
        "version": "RSDS_DOI_META_HETEROGENEITY_MISSINGNESS_V1",
        "research_stage": "POST_OUTCOME_METHOD_SENSITIVITY_NO_NEW_GPS",
        "source_git_blob_sha": SOURCE_GIT_BLOB_SHA,
        "original_studies": EXPECTED["studies"],
        "original_native_rows": EXPECTED["rows"],
        "limitations": [
            "Purposefully selected original GPS field studies, not probability samples",
            "Reported percentile tag intervals are converted to approximate SE",
            "Within-DOI covariance not recoverable; only rho=0 and rho=1 endpoints",
            "Study DOI may not certify biological source independence",
            "HOLD not missing at random; bound depends on explicit hypothetical range",
            "Q is observed GPS straight-line displacement, not latent activity",
            "No confirmatory p-values, uncorrected meta-regressions, or causal solar claims",
            "Old 48-original-source signed radiation nonreplication unchanged",
            "No original GPS events, tag keys, or exact location/time disclosed"
        ],
        "endpoints": {}
    }
    for field in ("q", "daylength"):
        s0 = study_effects(rows, field, 0.0)
        s1 = study_effects(rows, field, 1.0)
        if [s["doi"] for s in s0] != [s["doi"] for s in s1]:
            raise ValueError("Covariance scenarios changed the DOI universe")
        available = [s for s in s0 if s["estimate"] is not None]
        if len(available) != 9:
            raise ValueError("Unexpected supported study count")
        out["endpoints"][field] = {
            "unit": "log1p median GPS speed NIGHT-DAY" if field == "q"
                    else "within-tag-year q association per hour of astronomical daylength",
            "study_details_rho0": s0,
            "study_details_rho1": s1,
            "models": {
                "assumed_within_DOI_rho0": random_effects(s0),
                "assumed_within_DOI_rho1": random_effects(s1)
            },
            "missingness": missingness(rows, field, s0),
            "source_sign_counts": {
                "positive": sum(s["estimate"] > 0 for s in available),
                "negative": sum(s["estimate"] < 0 for s in available)
            }
        }
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, help="write JSON (default: stdout)")
    args = p.parse_args()
    result = analyze(read_pinned())
    rendered = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
