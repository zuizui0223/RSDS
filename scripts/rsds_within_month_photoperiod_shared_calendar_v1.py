"""Calendar-month conditional q~photoperiod slope on source-shared animal dates.

Only public source study group aggregates leave runtime. Astronomy daylength is
deterministically computed at a predeclared 43N reference; within-month
centering distinguishes associations based primarily on BETWEEN-month
seasonal shifts from associations also detectable WITHIN a source calendar
month. It does NOT identify biological causality from changing daylight.
"""
from __future__ import annotations

from collections import defaultdict
from statistics import mean
import math
import random

NBOOT=2000
SEED=20261029

def sign(ci):
    if ci is None: return "UNAVAILABLE"
    return "POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"


def sufficient_months(date_values, latitude, daylength_fn,
                      min_block_days=5, min_block_L_span=.05,
                      min_blocks=5, min_total_days=40, min_sxx=.25):
    """Private map (year, date) -> response q contrast for fixed source pair.

    No q value enters support-selection, which depends only on original
    common calendar dates and astronomical daylength.
    """
    if not date_values:
        return [], {
            "source_calendar_dates":0, "source_year_month_blocks":0,
            "eligible_source_month_blocks":0,
            "eligible_source_calendar_days":0,
            "within_month_daylength_centered_sum_squares_hours2":0.,
            "status":"INSUFFICIENT_COMMON_SOURCE_MONTH_SUPPORT"
        }
    groups=defaultdict(list)
    for (year,d),q in date_values.items():
        if not isinstance(q,(int,float)) or not math.isfinite(q):
            raise ValueError("Nonfinite original source taxon daily q difference")
        if d.year != year:
            raise ValueError("Year/date source identity changed")
        L=daylength_fn(latitude,d)
        if not math.isfinite(L):
            raise ValueError("Nonfinite astronomical daylength")
        groups[(year,d.month)].append((L,float(q)))
    usable=[]
    for key,series in sorted(groups.items()):
        x=[k[0] for k in series]
        if len(series)<min_block_days or max(x)-min(x)<min_block_L_span:
            continue
        mx=mean(x)
        my=mean(y for _,y in series)
        sxx=sum((a-mx)**2 for a,_ in series)
        sxy=sum((a-mx)*(b-my) for a,b in series)
        if sxx<=0: continue
        usable.append((len(series),sxx,sxy))
    n=sum(x[0] for x in usable)
    sxx_total=sum(x[1] for x in usable)
    ok=len(usable)>=min_blocks and n>=min_total_days and sxx_total>=min_sxx
    return usable,{
        "source_calendar_dates":len(date_values),
        "source_year_month_blocks":len(groups),
        "eligible_source_month_blocks":len(usable),
        "eligible_source_calendar_days":n,
        "within_month_daylength_centered_sum_squares_hours2":sxx_total,
        "status":"MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE" if ok else
                 "INSUFFICIENT_COMMON_SOURCE_MONTH_SUPPORT",
    }


def month_conditional_slope(date_values, latitude, daylength_fn, seed=SEED):
    blocks,qa=sufficient_months(date_values,latitude,daylength_fn)
    if qa["status"]!="MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE":
        return {**qa,"within_month_daylength_slope":None,
                "month_cluster_bootstrap95":None,
                "leave_one_month_block_out_range":None,
                "bootstrap_finite_replicates":0,
                "nominal_interval_sign":"UNAVAILABLE",
                "not_causal_daylength":True}
    sx=sum(x[1] for x in blocks)
    sy=sum(x[2] for x in blocks)
    beta=sy/sx
    if not math.isfinite(beta): raise ValueError("Nonfinite source fixed-month slope")
    loo=[]
    for n,sxx,sxy in blocks:
        if sx-sxx<=0:
            raise ValueError("One month uniquely determines source fixed-month slope")
        loo.append((sy-sxy)/(sx-sxx))
    rng=random.Random(seed)
    draws=[]
    size=len(blocks)
    for _ in range(NBOOT):
        chosen=[blocks[rng.randrange(size)] for j in range(size)]
        denom=sum(z[1] for z in chosen)
        if denom>0:
            v=sum(z[2] for z in chosen)/denom
            if math.isfinite(v):
                draws.append(v)
    if len(draws)<1900:
        ci=None
    else:
        draws.sort()
        ci=[draws[int(.025*len(draws))],
            draws[max(0,int(.975*len(draws))-1)]]
    return {
        **qa,
        "within_month_daylength_slope":beta,
        "month_cluster_bootstrap95":ci,
        "leave_one_month_block_out_range":[min(loo),max(loo)],
        "bootstrap_finite_replicates":len(draws),
        "nominal_interval_sign":sign(ci),
        "not_causal_daylength":True,
    }


def compare_frozen_datewise_prior(actual,before,tol=1e-12):
    """Assert previously fixed unadjusted season effect to machine precision."""
    if actual is None or before is None:
        if actual!=before: raise ValueError("Frozen original unadjusted daylength slope lost")
        return "BOTH_UNESTIMABLE"
    if not math.isfinite(actual) or not math.isclose(
        actual,before,rel_tol=0,abs_tol=tol
    ):
        raise ValueError("Archived original shared-calendar daylength contrast was modified")
    return "ARCHIVED_UNADJUSTED_SLOPE_REPLAYED"


def interpret_month_adjustment(previous_sign,current_sign):
    if current_sign=="UNAVAILABLE":
        return "WITHIN_MONTH_SUPPORT_UNAVAILABLE"
    if current_sign=="INCLUDES_ZERO":
        return "WITHIN_MONTH_INTERVAL_INCLUDES_ZERO"
    if previous_sign in ("POSITIVE","NOMINAL_POSITIVE") and current_sign=="POSITIVE":
        return "PREVIOUS_POSITIVE_ASSOCIATION_PERSISTS_WITHIN_MONTH_EXPLORATORY"
    if previous_sign in ("NEGATIVE","NOMINAL_NEGATIVE") and current_sign=="NEGATIVE":
        return "PREVIOUS_NEGATIVE_ASSOCIATION_PERSISTS_WITHIN_MONTH_EXPLORATORY"
    if previous_sign in ("POSITIVE","NEGATIVE","NOMINAL_POSITIVE","NOMINAL_NEGATIVE"):
        return "PREVIOUS_SIGN_NOT_RETAINED_WITHIN_MONTH"
    return "PREVIOUS_SLOPE_UNRESOLVED"
