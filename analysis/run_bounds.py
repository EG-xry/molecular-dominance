#!/usr/bin/env python3
"""Conditional profile bounds and an exploratory SCN5A application.

Writes only bounds.csv, bounds.json and Figure3.{png,pdf}.
This is a specialization of shape-restricted partial identification
(Mogstad, Santos & Torgovitsky 2018, doi:10.3982/ECTA15463).
For fixed p the coordinate extrema solve a linear program analytically.
Both extrema decrease with p on its feasible interval, by stochastic ordering
of the binomial expectations of the extremal monotone profiles. This gives
continuous interval-p extrema without a finite grid. Marginal extrema need
not coexist in one profile.

SCN5A assumes d=(W-H)/(W-s*A). The displayed control-based pair is
W=218.4, s=2.184, p=1-100/218.4. The additive calibration pair W=200,
s=2, p=.5 and historical W=218.4,s=2 are retained as named sensitivities.
The loss-of-function reading requires W-s*A>0 (m=s*A/W<1).
p ranges are analyst scenarios, not measured or probabilistic intervals.
At fixed q, d>q is assessed using C(q)=(1-q)*W+q*s*A-H>0.
W cancels at q=1, but dose uncertainty, shared normalization, batches and model
discrepancy do not disappear. All tests use independent original input means.
Ratio confidence sets invert Gaussian linear contrasts (Fieller 1954,
doi:10.1111/j.2517-6161.1954.tb00159.x). Published SEs are estimated, so coverage
is approximate. Covariance induced by W in numerator and denominator is included.
No activity bounds are calculated from uncalibrated Kv7 or TP53 readouts.
"""
import json
from datetime import date
from math import comb, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
P_LO, P_HI = 0.40, 0.60
INERT_CONTROL_P = 1 - 134.2 / 218.4
P_WIDE = (0.35, 0.65)
TOL = 1e-12


def _order(n):
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("assembly order must be an integer >= 2")


def _range(p_lo, p_hi):
    if not (np.isfinite(p_lo) and np.isfinite(p_hi) and 0 < p_lo <= p_hi < 1):
        raise ValueError("fractions must satisfy 0 < p_lo <= p_hi < 1")


def weights(n, p):
    _order(n)
    if not np.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("p must lie in [0,1]")
    return np.array([comb(n,k)*p**k*(1-p)**(n-k) for k in range(n+1)])


def bounds_at_p(n, p, d):
    """Sharp marginal bounds under 0=g0<=...<=gn=1 at interior p.
    Infeasible bounds encode conflicting inequalities; do not clip to a fit.
    """
    _range(p,p)
    if not np.isfinite(d):
        raise ValueError("d must be finite")
    w=weights(n,p)
    lo=np.array([max(0.,(d-w[k+1:].sum())/w[1:k+1].sum()) for k in range(1,n)])
    hi=np.array([min(1.,(d-w[n])/w[k:n].sum()) for k in range(1,n)])
    return lo,hi,bool(w[n]-TOL <= d <= 1-w[0]+TOL)


def feasible_band(n,p):
    w=weights(n,p)
    return float(w[n]),float(1-w[0])


def feasible_p_interval(n,d):
    """Compatibility interval in the closed domain; pure endpoints may occur."""
    _order(n)
    if not np.isfinite(d):
        raise ValueError("d must be finite")
    if not 0 <= d <= 1:
        return None
    lo=-np.expm1(np.log1p(-d)/n) if d < 1 else 1.
    return float(lo),float(d**(1/n))


def restricted_p_interval(n,d,p_lo=P_LO,p_hi=P_HI):
    _range(p_lo,p_hi)
    f=feasible_p_interval(n,d)
    if f is None:
        return None
    lo,hi=max(p_lo,f[0]),min(p_hi,f[1])
    if lo > hi+TOL:
        return None
    if lo > hi:
        lo=hi=(lo+hi)/2
    return float(lo),float(hi)


def bounds_over_p(n,d,p_lo=P_LO,p_hi=P_HI):
    """Exact marginal hull over feasible p; third value is its length fraction.
    This geometric fraction is not a probability. A feasible singleton in a
    nonzero-width scenario has fraction zero but retains nonempty bounds.
    """
    f=restricted_p_interval(n,d,p_lo,p_hi)
    if f is None:
        return np.full(n-1,np.nan),np.full(n-1,np.nan),0.
    a,b=f
    lo,hi=bounds_at_p(n,b,d)[0],bounds_at_p(n,a,d)[1]
    share=(b-a)/(p_hi-p_lo) if p_hi > p_lo else 1.
    return lo,hi,float(share)


def lp_bound(n,p,d,k,sense):
    if sense not in ("min","max") or not 1 <= k < n:
        raise ValueError("invalid objective")
    w=weights(n,p)
    c=np.zeros(n-1)
    c[k-1]=1 if sense=="min" else -1
    rows=np.eye(n-1)[:-1]-np.eye(n-1)[1:]
    r=linprog(c,A_ub=rows if len(rows) else None,
              b_ub=np.zeros(len(rows)) if len(rows) else None,
              A_eq=w[None,1:n],b_eq=[d-w[n]],bounds=[(0,1)]*(n-1),method="highs")
    return float(r.x[k-1]) if r.success else None


def lp_selfcheck():
    worst=0.
    for n in range(2,7):
        for p in (.1,.25,.4,.5,.6,.75,.9):
            for d in (.05,.2,.3645,.4227,.5342,.619,.8,.95):
                lo,hi,feasible=bounds_at_p(n,p,d)
                if feasible:
                    for k in range(1,n):
                        worst=max(worst,abs(lo[k-1]-lp_bound(n,p,d,k,"min")),
                                  abs(hi[k-1]-lp_bound(n,p,d,k,"max")))
    return float(worst)


def dimer_g1(p,d):
    p=np.asarray(p,float)
    if np.any(~np.isfinite(p)) or np.any((p<=0)|(p>=1)):
        raise ValueError("dimer inversion requires interior fractions")
    return (np.asarray(d)-p*p)/(2*p*(1-p))


def dimer_turning_point(d):
    if 0 <= d <= 1:
        return None
    r=sqrt(d*d-d)
    p=d/(d+r) if d > 1 else d/(d-r)
    return float(p),float(dimer_g1(p,d))


def dimer_g1_hull(d,p_lo=P_LO,p_hi=P_HI):
    """Unrestricted profile sensitivity range, possibly outside [0,1]."""
    _range(p_lo,p_hi)
    vals=[float(dimer_g1(p_lo,d)),float(dimer_g1(p_hi,d))]
    turn=dimer_turning_point(d)
    if turn is not None and p_lo < turn[0] < p_hi:
        vals.append(turn[1])
    return min(vals),max(vals)


def allele_ratio_tolerance(d):
    """Lower compatibility boundary, only defined for 0<=d<=1.
    At equality g1=1; strict upper violation occurs below the boundary.
    g1>1 means a1<m only when m<1.
    """
    a=np.asarray(d,float)
    return np.where((a>=0)&(a<=1),1-np.sqrt(np.clip(1-a,0,None)),np.nan)


def supra_wildtype_tolerance(d):
    a=np.asarray(d,float)
    return np.where((a>=0)&(a<=1),np.sqrt(np.clip(a,0,None)),np.nan)


def ratio_confidence_set(x,y,vx,vy,cov,z):
    """Invert (y-q*x)^2 <= z^2 Var(y-q*x), retaining disjoint/unbounded sets.
    Estimated input variances and normal approximation give approximate coverage.
    """
    if vx < 0 or vy < 0 or cov*cov > vx*vy+TOL or z <= 0:
        raise ValueError("invalid covariance or critical value")
    a,b,c=x*x-z*z*vx,-2*x*y+2*z*z*cov,y*y-z*z*vy
    eps=np.finfo(float).eps*max(1,abs(a),abs(b),abs(c))*8
    if abs(a)<=eps:
        if abs(b)<=eps:
            return [(-np.inf,np.inf)] if c<=eps else []
        return [(-np.inf,-c/b)] if b>0 else [(-c/b,np.inf)]
    disc=b*b-4*a*c
    if disc<0:
        return [(-np.inf,np.inf)] if a<0 else []
    r1,r2=sorted(((-b-sqrt(disc))/(2*a),(-b+sqrt(disc))/(2*a)))
    return [(r1,r2)] if a>0 else [(-np.inf,r1),(r2,np.inf)]


def threshold_test(frame,W,Wse,Wn,scale,q,discrepancy_sd=0.):
    """One-sided C(q)>0 normal/Welch sensitivity tests from independent means."""
    coeff=np.array([1-q,q*scale,-1.])
    means=np.vstack([np.full(len(frame),W),frame.homo,frame.het])
    variances=np.vstack([np.full(len(frame),Wse**2),frame.homo_se**2,frame.het_se**2])
    ns=np.vstack([np.full(len(frame),Wn),frame.homo_n,frame.het_n])
    parts=coeff[:,None]**2*variances
    contrast=coeff@means
    se=np.sqrt(parts.sum(axis=0)+discrepancy_sd**2)
    dof=parts.sum(axis=0)**2/(parts**2/(ns-1)).sum(axis=0)
    z=contrast/se
    return {"contrast":contrast,"se":se,"z":z,"df":dof,
            "normal_p":stats.norm.sf(z),"welch_p":stats.t.sf(z,dof)}


def _classify(d):
    if d>1:
        return "above unit range at point estimate"
    if d<0:
        return "below unit range at point estimate"
    f=restricted_p_interval(2,d)
    if f is None:
        return "no monotone profile in scenario at point estimate"
    if f==(P_LO,P_HI):
        return "all scenario fractions compatible at point estimate"
    return "some scenario fractions compatible at point estimate"


def scn5a_table(scale=2.184, reference=None, reference_se=None, p_reference=None):
    if not np.isfinite(scale) or scale<=0:
        raise ValueError("dose scale must be fixed and positive")
    t=pd.read_csv(HERE/"scn5a_peak.csv")
    ctl=pd.read_csv(HERE/"scn5a_controls.csv").set_index("variant").loc["WT+WT"]
    W,Wse,Wn=float(ctl.het),float(ctl.het_se),int(ctl.het_n)
    if reference is not None:
        W = float(reference)
        Wse = float(reference_se if reference_se is not None else 0.)
    if p_reference is None:
        p_reference = 1 - 100 / W
    t["x"],t["y"]=W-scale*t.homo,W-t.het
    if (t.x<=0).any():
        raise ValueError("loss-of-function interpretation requires W-scale*alone>0")
    t["d"]=t.y/t.x
    t["m_assumed"]=scale*t.homo/W
    t["p_star"]=allele_ratio_tolerance(t.d)
    t["p_zero"]=supra_wildtype_tolerance(t.d)
    t["scenario_class"]=t.d.map(_classify)
    vx,vy,cov=Wse**2+(scale*t.homo_se)**2,Wse**2+t.het_se**2,Wse**2
    zs={"95":stats.norm.ppf(.975),"sim95":stats.norm.ppf(1-.05/(2*len(t)))}
    for name,z in zs.items():
        sets=[ratio_confidence_set(x,y,a,b,cov,z) for x,y,a,b in zip(t.x,t.y,vx,vy)]
        bounded=[len(v)==1 and np.isfinite(v[0]).all() for v in sets]
        t[f"d_{name}_bounded"]=bounded
        t[f"d_{name}_lo"]=[v[0][0] if ok else np.nan for v,ok in zip(sets,bounded)]
        t[f"d_{name}_hi"]=[v[0][1] if ok else np.nan for v,ok in zip(sets,bounded)]
    point=[dimer_g1_hull(d) for d in t.d]
    t["g1_hull_lo"],t["g1_hull_hi"]=np.array(point).T
    # g1 increases with d at fixed p. Union over d confidence set and the p
    # scenario gives a conditional envelope, distinct from the pointwise hull.
    t["g1_95_envelope_lo"]=[dimer_g1_hull(d)[0] if np.isfinite(d) else -np.inf for d in t.d_95_lo]
    t["g1_95_envelope_hi"]=[dimer_g1_hull(d)[1] if np.isfinite(d) else np.inf for d in t.d_95_hi]
    for p in (P_LO,.5,P_HI):
        t[f"g1_at_p{round(100*p)}"]=dimer_g1(p,t.d)
    t["g1_at_reference_p"] = dimer_g1(p_reference, t.d)
    t["g1_at_inert_control_p"] = dimer_g1(INERT_CONTROL_P, t.d)
    allowed=[restricted_p_interval(2,d) for d in t.d]
    t["monotone_p_lo"]=[v[0] if v is not None else np.nan for v in allowed]
    t["monotone_p_hi"]=[v[1] if v is not None else np.nan for v in allowed]
    for name,q in (("unit",1.),("scenario_upper",1-(1-P_HI)**2)):
        result=threshold_test(t,W,Wse,Wn,scale,q)
        for key,values in result.items():
            t[f"{name}_{key}"]=values
        for method in ("normal","welch"):
            t[f"{name}_{method}_bonf"]=np.minimum(1.,len(t)*result[f"{method}_p"])
    return t,{"WW":W,"WW_se":Wse,"WW_n":Wn,"scale":scale,"p_reference":p_reference,
              "critical_values":{k:float(v) for k,v in zs.items()},
              "min_simultaneous_denominator_lower":float((t.x-zs["sim95"]*np.sqrt(vx)).min())}


def _test_summary(frame,result):
    return {"point_count":int((result["contrast"]>0).sum()),
            "normal_nominal_set":sorted(frame.loc[result["normal_p"]<.05,"variant"]),
            "normal_bonferroni_set":sorted(frame.loc[result["normal_p"]<.05/len(frame),"variant"]),
            "welch_bonferroni_set":sorted(frame.loc[result["welch_p"]<.05/len(frame),"variant"])}


def dose_threshold(frame,variant,family,alpha=.05):
    """Smallest dose factor s at which the Bonferroni adjusted one sided test of
    s*A-H>0 rejects, under the normal and the Welch approximations.

    At q=1 the reference W cancels from the contrast, so only the alone and
    mixed means enter. The family is the same 51 variants as every other test.
    """
    from scipy.optimize import brentq
    r=frame.set_index("variant").loc[variant]
    def se(s_):
        return sqrt((s_*r.homo_se)**2+r.het_se**2)
    def dof(s_):
        parts=np.array([(s_*r.homo_se)**2,r.het_se**2])
        return parts.sum()**2/(parts**2/(np.array([r.homo_n,r.het_n])-1)).sum()
    z=stats.norm.ppf(1-alpha/family)
    lo,hi=1.,2.184
    normal=brentq(lambda s_:s_*r.homo-r.het-z*se(s_),lo,hi)
    welch=brentq(lambda s_:s_*r.homo-r.het-stats.t.ppf(1-alpha/family,dof(s_))*se(s_),lo,hi)
    return {"variant":variant,"normal":float(normal),"welch":float(welch),
            "critical_normal":float(z),"family_size":int(family),
            "definition":("smallest s at which the one sided lower confidence bound for s*A-H, "
                          "Bonferroni adjusted across the family, reaches zero")}


def figure3(sc,scale_rows,threshold,ctl):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from palette import (AXIS, BLUE, BLUE_DARK, BLUE_LIGHT, FIG_WIDTH, GRID, INK, PILL, TEAL,
                         TEAL_LIGHT, marker_handle, panel_label, rounded_bar, soft_grid, tint,
                         use_figure_style)
    from hgvs import scn5a_hgvs
    use_figure_style()
    fig=plt.figure(figsize=(FIG_WIDTH,3.5))
    # Fixed positions leave room below each panel for its legend.
    a=fig.add_axes([0.135,0.335,0.365,0.565])
    b=fig.add_axes([0.62,0.335,0.365,0.565])
    ring=dict(mec="white",mew=0.6)

    # A: paired W=218.4, s=2.184 reference for the four variants discussed.
    names=["D785N","D356N","W822X","E746K"]
    rows=sc.set_index("variant").loc[names]
    a.set_xlim(-1.35,3.7)
    a.set_ylim(len(names)-.5,-.5)
    a.axvspan(0,1,color=tint(TEAL_LIGHT,0.24),lw=0,zorder=0)
    a.hlines(range(len(names)),*a.get_xlim(),color=GRID,lw=.6,zorder=1)
    for edge in (0,1):
        a.axvline(edge,color=TEAL,lw=.6,alpha=.55,zorder=1)
    for i,row in enumerate(rows.itertuples()):
        rounded_bar(a,row.g1_95_envelope_lo,row.g1_95_envelope_hi,i,0.26,radius=2.2,both=True,
                    color=BLUE_LIGHT,zorder=2)
        a.plot([row.g1_hull_lo,row.g1_hull_hi],[i,i],color=BLUE_DARK,lw=2.6,
               solid_capstyle="butt",zorder=3)
        a.plot(row.g1_at_reference_p,i,"o",color=BLUE_DARK,ms=5,zorder=4,**ring)
        a.plot(row.g1_at_inert_control_p,i,"D",color=TEAL,ms=4.5,zorder=5,**ring)
    a.set_yticks(range(len(names)),[scn5a_hgvs(n) for n in names])
    a.tick_params(axis="y",length=0)
    a.spines["left"].set_visible(False)
    a.set_xlabel("conditional dimer profile entry $g_1$")
    panel_label(a,"A",[("paired scale: ",False),("s",True),(" = 2.184, ",False),("W",True),(" = 218.4",False)], size=8)
    a.legend(handles=[Line2D([],[],color=BLUE_DARK,lw=2.6),marker_handle("o",BLUE_DARK,5,**ring),
                      Patch(facecolor=BLUE_LIGHT,lw=0),marker_handle("D",TEAL,4.5,**ring),
                      Patch(facecolor=tint(TEAL_LIGHT,0.24),lw=0)],
             labels=["$p=0.40$ to $0.60$",f"WT+WT: $p={ctl['p_reference']:.3f}$",
                     "conditional 95% envelope",f"inert control: $p={INERT_CONTROL_P:.3f}$",
                     "$0\\leq g_1\\leq 1$"],
             loc="upper center",bbox_to_anchor=(0.45,-0.23),ncol=2,
             handlelength=1.6,columnspacing=1.0,labelspacing=0.4)

    # B: the s*A-H contrast for p.Asp785Asn as a continuous function of the
    # assumed dose factor, with the Bonferroni adjusted one sided lower bound
    # that the reported test uses. Dots mark the four tabulated scenarios.
    r=sc.set_index("variant").loc["D785N"]
    grid=np.linspace(1.,2.25,251)
    est=grid*r.homo-r.het
    se=np.sqrt((grid*r.homo_se)**2+r.het_se**2)
    lower=est-threshold["critical_normal"]*se
    scales=np.array([q["scale"] for q in scale_rows])
    contrasts=np.array([q["D785N"]["contrast"] for q in scale_rows])
    if not np.allclose(contrasts,scales*r.homo-r.het):
        raise ValueError("tabulated contrasts disagree with the continuous curve")
    soft_grid(b,"y")
    b.fill_between(grid,lower,est,color=tint(BLUE_LIGHT,0.45),lw=0,zorder=1)
    b.axhline(0,color=AXIS,lw=.7,zorder=2)
    b.plot(grid,est,color=BLUE,lw=1.5,zorder=3)
    b.plot(grid,lower,color=TEAL,lw=1.2,ls=(0,(4,2)),zorder=3)
    b.plot(scales,contrasts,"o",color=BLUE,ms=5.4,zorder=5,**ring)
    s_star=threshold["normal"]
    b.plot([s_star],[0],"o",mfc="white",mec=TEAL,mew=1.3,ms=5.6,zorder=6)
    b.annotate(f"$s$ = {s_star:.2f}",(s_star,0),xytext=(7,-12),textcoords="offset points",
               fontsize=8,color=INK,ha="left",va="center",bbox=PILL)
    b.set_xticks(scales,[f"{q:g}" for q in scales])
    b.set_xlim(.95,2.25)
    b.set_xlabel("assumed combined scale factor $s$")
    b.set_ylabel("$s\\,I_{\\mathrm{alone}}-I_{\\mathrm{mixed}}$ (percentage points)")
    panel_label(b,"B","p.Asp785Asn: the scale assumption")
    b.legend(handles=[Line2D([],[],color=BLUE,lw=1.5,marker="o",ms=5.4,**ring),
                      Line2D([],[],color=TEAL,lw=1.2,ls=(0,(4,2)))],
             labels=["estimate (dots: fixed scenarios)",
                     "one sided 95% lower bound, Bonferroni across 51"],
             loc="upper center",bbox_to_anchor=(0.45,-0.23),ncol=1,
             handlelength=2.0,labelspacing=0.4)
    fig.savefig(HERE/"Figure3.png",dpi=300)
    fig.savefig(HERE/"Figure3.pdf")
    plt.close(fig)


def main():
    sc,ctl=scn5a_table()
    W,Wse,Wn,s=(ctl[k] for k in ("WW","WW_se","WW_n","scale"))
    profile=pd.read_csv(HERE/"profile_scn5a.csv").set_index("variant").g1_p50
    legacy, _ = scn5a_table(scale=2., p_reference=.5)
    check_error=float((legacy.set_index("variant").g1_at_p50.reindex(profile.index)-profile).abs().max())
    paired, paired_frames = {}, []
    for name, scale, reference, reference_se, p_reference in (
            ("observed WT+WT reference with matched scale", 2.184, 218.4, 7.7, 1-100/218.4),
            ("additive reference with matched scale", 2., 200., 0., .5),
            ("inert p.Trp822Ter fraction with matched scale", 2.184, 218.4, 7.7, INERT_CONTROL_P),
            ("historical unmatched scale sensitivity", 2., 218.4, 7.7, .5)):
        table, control = scn5a_table(scale, reference, reference_se, p_reference)
        panel = table[table.variant != "W822X"]
        x, y = panel.x.to_numpy(), panel.y.to_numpy()
        slope = float(np.dot(x, y) / np.dot(x, x))
        positions = panel.variant.str.extract(r"(\d+)")[0].astype(int).to_numpy()
        groups = [np.where(positions == position)[0] for position in np.unique(positions)]
        rng = np.random.default_rng(20261120)
        slopes = []
        for _ in range(20000):
            indices = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
            slopes.append(float(np.dot(x[indices], y[indices]) / np.dot(x[indices], x[indices])))
        ci = np.quantile(slopes, [.025, .975])
        contrast = threshold_test(table, reference, reference_se, Wn, scale, 1.)
        k = table.index[table.variant == "D785N"][0]
        paired[name] = {"reference": reference, "reference_se": reference_se, "scale": scale,
                        "p": p_reference, "slope": slope, "slope_ci": ci.tolist(),
                        "implied_g1": float(dimer_g1(p_reference, slope)),
                        "implied_g1_ci": [float(dimer_g1(p_reference, q)) for q in ci],
                        "R104W": {"g1": float(table.set_index("variant").loc["R104W", "g1_at_reference_p"]),
                                  "mixed_mean": float(table.set_index("variant").loc["R104W", "het"]),
                                  "mixed_se": float(table.set_index("variant").loc["R104W", "het_se"]),
                                  "alone_mean": float(table.set_index("variant").loc["R104W", "homo"]),
                                  "alone_se": float(table.set_index("variant").loc["R104W", "homo_se"])},
                        "D785N": {key: float(value[k]) for key, value in contrast.items()},
                        "D785N_normal_bonferroni": float(min(1, len(table) * contrast['normal_p'][k])),
                        "D785N_welch_bonferroni": float(min(1, len(table) * contrast['welch_p'][k])),
                        "table_J": table.set_index("variant").loc[["D785N", "D356N", "W822X", "E746K"],
                            ["d", "g1_at_reference_p", "g1_at_inert_control_p", "g1_hull_lo", "g1_hull_hi", "g1_95_envelope_lo", "g1_95_envelope_hi"]].to_dict("index"),
                        "limits": "Fixed combined reference-conversion/dose factor and assumed fraction; controls do not measure mutant scaling. Same mixed-current observations retained under each reference assumption."}
        paired_frames.append(table.assign(scenario=name, reference=reference, scale=scale, p_reference=p_reference))
    pd.concat(paired_frames, ignore_index=True).to_csv(HERE/"scn5a_reference_dose.csv",index=False)
    scale_rows=[]
    for scale in (1.,1.5,2.,2.184):
        result=threshold_test(sc,W,Wse,Wn,scale,1.)
        k=sc.index[sc.variant=="D785N"][0]
        scale_rows.append({"scale":scale,**_test_summary(sc,result),
                           "D785N":{key:float(value[k]) for key,value in result.items()}})
    threshold=dose_threshold(sc,"D785N",len(sc))
    scenarios=[]
    for n in (2,3,4):
        for low,high in ((P_LO,P_HI),P_WIDE):
            q=1-(1-high)**n
            result=threshold_test(sc,W,Wse,Wn,s,q)
            scenarios.append({"n":n,"p_range":[low,high],"upper_feasible_d":q,
                              **_test_summary(sc,result)})
    discrepancy=[]
    for tau in (0.,5.,10.,15.):
        result=threshold_test(sc,W,Wse,Wn,s,1.,tau)
        discrepancy.append({"additional_sd":tau,"normal_bonferroni_set":
                            sorted(sc.loc[result["normal_p"]<.05/len(sc),"variant"])})
    examples=[]
    for n,d in ((4,.5),(4,.42)):
        lo,hi,share=bounds_over_p(n,d)
        examples.append({"label":"synthetic mathematical example","n":n,"d":d,
                         "p_range":[P_LO,P_HI],"lower":lo.tolist(),"upper":hi.tolist(),
                         "feasible_length_fraction":share})
    report={"generated":date.today().isoformat(),"method":"conditional shape-restricted partial identification",
            "closed_form_vs_lp_max_error":lp_selfcheck(),
            "g1_at_p_half_vs_profile_scn5a_max_error":check_error,
            "profile_crosscheck_scope":"historical scale2/reference218.4 only; displayed baseline uses scale2.184",
            "family_size":len(sc),"p_scenario":[P_LO,P_HI],"scn5a_scale":ctl,
            "point_d_above_one":sorted(sc.loc[sc.d>1,"variant"]),
            "point_d_below_zero":sorted(sc.loc[sc.d<0,"variant"]),
            "scenario_class_counts":sc.scenario_class.value_counts().to_dict(),
            "dose_scaling_sensitivity":scale_rows,"dose_threshold":threshold,
            "paired_reference_dose_scenarios": paired,
            "paired_bootstrap": {"seed": 20261120, "draws": 20000, "unit": "residue", "scope": "fixed published measurements; shared draws across scenarios"},
            "inert_control_profile_fraction": {"p": INERT_CONTROL_P,
                "inside_displayed_fraction_range": False,
                "figure_marker": "teal diamond; estimates at this fixed p are additional to the p=.40-.60 range and envelope"},
            "fraction_and_order_sensitivity":scenarios,
            "model_discrepancy_sensitivity":discrepancy,"synthetic_examples":examples,
            "ratio_sets_bounded":bool(sc.d_95_bounded.all() and sc.d_sim95_bounded.all()),
            "references":{"partial_identification":"https://doi.org/10.3982/ECTA15463",
                          "ratio_intervals":"https://doi.org/10.1111/j.2517-6161.1954.tb00159.x"},
            "interpretation_limits":[
                "Fraction ranges are analyst scenarios, not measured confidence intervals.",
                "The displayed W=218.4,s=2.184 pair is evaluated at both p=.5421 (WT+WT) and p=.3855 (inert control); the additive pair is W=200,s=2,p=.5. These all assume fixed fractions and linear mutant-dose scaling. Historical W=218.4,s=2 remains explicitly a scale sensitivity.",
                "W-s*A>0 (m<1) is required for the stated loss-of-function sign interpretation.",
                "Known-p dimer inversion needs no monotonicity; g1 in [0,1] is the hypothesis assessed.",
                "For n>2 the bounds impose 0=g0<=...<=gn=1; marginal extrema need not coexist.",
                "d>1 implies some interior gk>1 under convex assembly weights, not specifically g1 for n>2.",
                "Normal/Welch tests and ratio sets assume independent input means and adequate published SEs.",
                "The reference W induces numerator-denominator covariance and is included.",
                "Shared normalization, batches, covariance among original inputs and dose-factor uncertainty remain omitted.",
                "Each sensitivity specification corrects across 51 variants; no combined family across scenarios is claimed.",
                "Empty point-estimate sets do not by themselves constitute statistical rejection or establish a mechanism.",
                "No activity bounds are reported for uncalibrated Kv7 or TP53 axes."]}
    sc.to_csv(HERE/"bounds.csv",index=False)
    (HERE/"bounds.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    figure3(sc,scale_rows,threshold,ctl)
    print(json.dumps(report,indent=2,allow_nan=False))
    print("\nWrote bounds.csv, bounds.json and Figure3.{png,pdf}.")


if __name__=="__main__":
    main()
