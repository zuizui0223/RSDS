"""Joint month/block covariance for NEW source-selected native GPS effects.

New 16-slot source frame is pre-event fixed. Conditional on multi-native
common original GPS calendar support, every species and both axes share
the SAME source month resamples. The positive animal weights prevent loss of
source calendar support; their uncertainty is Bayesian-bootstrap-style and
does not have calibrated frequentist coverage. Never serialize dates, GPS,
tag IDs or individual draws.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
import math
import random
from statistics import mean

from scripts.movebank_native_ecology_20261019 import daylength
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import sufficient_months

B=1000
SEED=20261115

def private_daily(common, dates, taxa):
    if len(taxa)<2 or len(set(taxa))!=len(taxa) or set(common)!=set(taxa):
        raise ValueError("Incomplete frozen registered multi-native source")
    dates=set(dates)
    if len(dates)<40 or any(y!=d.year for y,d in dates):
        raise ValueError("Insufficient or invalid source dates")
    vals={};counts={}
    for sp in taxa:
        if len(common[sp])<5:
            raise ValueError("Source registered native taxon has under five original tags")
        keyed={}
        for i,(tag,years) in enumerate(sorted(common[sp].items())):
            for y,day in years.items():
                for d,tuple_value in day.items():
                    if (y,d) not in dates:
                        raise ValueError("Different native source calendar selection")
                    q=float(tuple_value[0])
                    if not math.isfinite(q):
                        raise ValueError("Invalid original native daily q")
                    keyed.setdefault((y,d),[]).append((i,q))
        if set(keyed)!=dates:
            raise ValueError("Source native does not span all common dates")
        vals[sp]=keyed
        counts[sp]=len(common[sp])
    return vals,counts

def month_blocks(dates,latitude):
    report_input={key:0.0 for key in dates}
    _,r=sufficient_months(report_input,latitude,daylength)
    if r["status"]!="MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE":
        raise ValueError("Joint covariance month-support HOLD")
    grouped=defaultdict(list)
    for y,d in dates:
        grouped[(y,d.month)].append((y,d,daylength(latitude,d)))
    blocks=[]
    for k,rows in sorted(grouped.items()):
        xs=[z[2] for z in rows]
        if len(xs)<5 or max(xs)-min(xs)<.05:
            continue
        mu=mean(xs);sxx=sum((x-mu)**2 for x in xs)
        if sxx>0:
            blocks.append(([(y,d,x-mu) for y,d,x in rows],sxx,len(xs)))
    if (len(blocks)!=r["eligible_source_month_blocks"]
        or sum(b[2] for b in blocks)!=r["eligible_source_calendar_days"]
        or not math.isclose(sum(b[1] for b in blocks),
                            r["within_month_daylength_centered_sum_squares_hours2"],
                            rel_tol=0,abs_tol=1e-9)):
        raise ValueError("Source original month estimand definition shifted")
    return blocks,r

def vector(blocks,vals,taxa,weights,selected):
    n=0;sxx=0;K=len(taxa)
    y=[0.0]*K;xy=[0.0]*K
    for ix in selected:
        rows,ss,n_rows=blocks[ix]
        n+=n_rows;sxx+=ss
        for yr,date,x in rows:
            key=(yr,date)
            for j,sp in enumerate(taxa):
                samples=vals[sp][key]
                den=sum(weights[sp][tag] for tag,q in samples)
                if den<=0:raise ValueError("No source month tag weight support")
                v=sum(weights[sp][tag]*q for tag,q in samples)/den
                y[j]+=v;xy[j]+=x*v
    if n<=0 or sxx<=0:raise ValueError("Degenerate source month resampling")
    return [v/n for v in y]+[a/sxx for a in xy]

def covariance(draws):
    if len(draws)<3:
        raise ValueError("Too few source joint draws")
    d=len(draws[0]);mu=[mean(row[j] for row in draws) for j in range(d)]
    out=[[sum((x[i]-mu[i])*(x[j]-mu[j]) for x in draws)/(len(draws)-1)
          for j in range(d)] for i in range(d)]
    if any(z<0 or not math.isfinite(z) for i,row in enumerate(out) for j,z in enumerate(row)
           if i==j):
        raise ValueError("Nonfinite/negative joint source variance")
    return out

def quadratic(V,c):
    return sum(c[i]*V[i][j]*c[j] for i in range(len(c)) for j in range(len(c)))

def original_joint(common,dates,taxa,latitude,seed=SEED,replicates=B):
    if not 100<=replicates<=4000:
        raise ValueError("Out-of-contract source covariance replication count")
    vals,counts=private_daily(common,dates,taxa)
    blocks,qa=month_blocks(dates,latitude)
    idx=list(range(len(blocks)))
    base={sp:[1.0]*counts[sp] for sp in taxa}
    q=vector(blocks,vals,taxa,base,idx)
    K=len(taxa)
    modes={}
    for mode_index,mode in enumerate(("month_only","month_plus_positive_tag_weights")):
        rng=random.Random(seed+mode_index*100003)
        draws=[]
        for _ in range(replicates):
            sample=[rng.randrange(len(blocks)) for z in blocks]
            w=(base if mode=="month_only" else
               {sp:[rng.expovariate(1.0) for i in range(counts[sp])]
                for sp in taxa})
            draws.append(vector(blocks,vals,taxa,w,sample))
        V=covariance(draws)
        # covariance from centered complete joint vectors is PSD, validate
        # independent signed contrasts instead of source-native independence.
        for a in range(2*K):
            for b in range(a+1,2*K):
                c=[0.0]*(2*K);c[a]=1;c[b]=-1
                if quadratic(V,c)<-1e-9:
                    raise ValueError("Negative joint source difference variance")
        pairs=[]
        for i,j in combinations(range(K),2):
            coef=[0.0]*(2*K);coef[i]=1;coef[j]=-1
            slope=[0.0]*(2*K);slope[i+K]=1;slope[j+K]=-1
            pairs.append({
                "taxa":[taxa[i],taxa[j]],
                "source_common_month_q_difference":q[i]-q[j],
                "q_difference_joint_standard_deviation":math.sqrt(max(0,quadratic(V,coef))),
                "month_daylength_slope_difference_per_hour":q[i+K]-q[j+K],
                "month_slope_difference_joint_standard_deviation":math.sqrt(max(0,quadratic(V,slope))),
                "NOT_independent_source_studies":True
            })
        modes[mode]={
            "original_source_native_covariance_matrix":V,
            "replicates":replicates,
            "species_pair_contrasts":pairs,
        }
    return {
        "status":"SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED",
        "original_native_taxa":list(taxa),
        "source_tag_counts":counts,
        "original_common_calendar_days":len(dates),
        "eligible_common_month_days":qa["eligible_source_calendar_days"],
        "eligible_common_month_blocks":qa["eligible_source_month_blocks"],
        "vector_order":["q_date_equal:"+s for s in taxa]+[
            "q_month_daylength_slope_per_hour:"+s for s in taxa],
        "source_date_equal_effect_vector":q,
        "joint_covariance_modes":modes,
        "neither_global_meta_analysis_nor_activity_duration":True,
        "source_dates_tags_positions_and_draws_never_exported":True,
    }
