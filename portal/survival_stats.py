"""Survival statistics helpers for the portal.

Per-gene analysis convention (as in the manuscript figures):
  * DISPLAY: Kaplan-Meier at the optimal cutpoint (maximally selected log-rank
    statistic, cutpoints scanned between the 10th and 90th expression
    percentiles, >=10 samples per arm).
  * REPORT: Cox proportional-hazards regression on CONTINUOUS expression —
    likelihood-ratio test p value and hazard ratio (per 1 SD of expression, so
    HRs are comparable across genes).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test


def cox_continuous(df: pd.DataFrame, time: str, status: str, value: str = "value") -> dict:
    """Cox PH on z-scored continuous expression. Returns HR per SD, CI, LR p."""
    d = df[[time, status, value]].dropna().copy()
    if len(d) < 10 or d[status].sum() < 3 or d[value].std() == 0:
        return {}
    d["z"] = (d[value] - d[value].mean()) / d[value].std()
    cph = CoxPHFitter()
    cph.fit(d[[time, status, "z"]], duration_col=time, event_col=status)
    s = cph.summary.loc["z"]
    return {
        "n": len(d),
        "events": int(d[status].sum()),
        "hr": float(np.exp(s["coef"])),
        "hr_lo": float(np.exp(s["coef lower 95%"])),
        "hr_hi": float(np.exp(s["coef upper 95%"])),
        "lr_p": float(cph.log_likelihood_ratio_test().p_value),
        "wald_p": float(s["p"]),
    }


def optimal_cutpoint(df: pd.DataFrame, time: str, status: str, value: str = "value",
                     min_per_arm: int = 10) -> dict:
    """Maximally selected log-rank cutpoint (surv_cutpoint-style).

    Scans unique expression values between the 10th and 90th percentiles,
    keeps splits with >= min_per_arm samples in each arm, and returns the
    cutpoint maximizing the log-rank statistic.
    """
    d = df[[time, status, value]].dropna()
    if len(d) < 2 * min_per_arm:
        return {}
    lo, hi = d[value].quantile([0.10, 0.90])
    candidates = np.unique(d[value][(d[value] >= lo) & (d[value] <= hi)])
    if len(candidates) > 150:  # grid-thin for speed; resolution is ample
        candidates = np.quantile(candidates, np.linspace(0, 1, 150))
    best = {}
    for c in candidates:
        high = d[value] > c
        if high.sum() < min_per_arm or (~high).sum() < min_per_arm:
            continue
        r = logrank_test(d.loc[high, time], d.loc[~high, time],
                         d.loc[high, status], d.loc[~high, status])
        if not best or r.test_statistic > best["stat"]:
            best = {"cutpoint": float(c), "stat": float(r.test_statistic),
                    "logrank_p": float(r.p_value),
                    "n_high": int(high.sum()), "n_low": int((~high).sum())}
    return best
