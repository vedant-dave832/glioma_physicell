"""Reduce trajectories to one row per replicate, and compare against theory.

Reviewer major issue #4: "Each complete run, rather than each timepoint, should
be treated as an independent replicate. Select one primary outcome at the
replicate level." This module is that reduction. Everything downstream of it
operates on runs, never on timepoints, so the dependence between successive
timepoints cannot leak into a p-value.

It also computes what each arm is EXPECTED to do from its parameters alone. If
the simulation matches the closed-form expectation, the simulation has confirmed
its own bookkeeping and nothing about biology -- which is the honest reading of
this experiment and the thing the manuscript needs to say.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import convert
from .registry import Model

PRIMARY_OUTCOME = "log_growth_rate_per_min"

SECONDARY_OUTCOMES = [
    "n_end",
    "log2_fold_change",
    "observed_doubling_time_days",
    "auc_log_count",
    "final_proportion",
    "final_radius_of_gyration",
]


def _fit_log_growth(times_min: np.ndarray, counts: np.ndarray) -> float:
    """OLS slope of ln(count) on time. NaN if fewer than 3 positive counts.

    Fitted on the whole trajectory, not on a hand-picked window. If crowding
    bends the curve, that shows up as a lower fitted rate in every arm equally.
    """
    ok = counts > 0
    if ok.sum() < 3:
        return float("nan")
    return float(np.polyfit(times_min[ok], np.log(counts[ok]), 1)[0])


def _shannon(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return float("nan")
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def per_replicate(traj: pd.DataFrame) -> pd.DataFrame:
    """One row per (arm, seed, population)."""
    rows = []
    for (arm, seed, pop), g in traj.groupby(["arm", "seed", "population"]):
        g = g.sort_values("time_min")
        t = g["time_min"].to_numpy(float)
        n = g["n_live"].to_numpy(float)
        n0, n_end = n[0], n[-1]

        arm_seed = traj[(traj.arm == arm) & (traj.seed == seed)]
        final = arm_seed[arm_seed.frame == arm_seed.frame.max()]
        total_final = final["n_live"].sum()

        rate = _fit_log_growth(t, n)
        rows.append(
            {
                "arm": arm,
                "seed": seed,
                "population": pop,
                "apoptosis_factor": g["apoptosis_factor"].iloc[0],
                "n_start": n0,
                "n_end": n_end,
                "n_dead_end": g["n_dead"].iloc[-1],
                "fold_change": n_end / n0 if n0 else np.nan,
                "log2_fold_change": np.log2(n_end / n0) if n0 and n_end else np.nan,
                PRIMARY_OUTCOME: rate,
                "observed_doubling_time_days": (
                    convert.min_to_days(np.log(2) / rate)
                    if rate and rate > 0
                    else np.nan
                ),
                # AUC of log-count: sensitive to the whole trajectory, not just
                # the endpoint, so a run that grows then crashes is not scored
                # the same as one that grows steadily.
                "auc_log_count": float(
                    np.trapezoid(np.log(np.maximum(n, 1.0)), t) / (t[-1] - t[0])
                )
                if len(t) > 1
                else np.nan,
                "final_proportion": n_end / total_final if total_final else np.nan,
                "final_radius_of_gyration": g["radius_of_gyration"].iloc[-1],
                "min_oxygen_mmHg": g["min_oxygen_mmHg"].min(),
                "n_timepoints_observed": len(g),
            }
        )
    return pd.DataFrame(rows)


def per_replicate_arm(traj: pd.DataFrame) -> pd.DataFrame:
    """One row per (arm, seed): whole-community outcomes, chiefly diversity."""
    rows = []
    for (arm, seed), g in traj.groupby(["arm", "seed"]):
        first = g[g.frame == g.frame.min()]
        last = g[g.frame == g.frame.max()]
        h0 = _shannon(first["n_live"].to_numpy(float))
        h1 = _shannon(last["n_live"].to_numpy(float))
        rows.append(
            {
                "arm": arm,
                "seed": seed,
                "shannon_start": h0,
                "shannon_end": h1,
                "shannon_change": h1 - h0,
                "n_total_end": last["n_live"].sum(),
                "min_oxygen_mmHg": g["min_oxygen_mmHg"].min(),
            }
        )
    return pd.DataFrame(rows)


def analytic_expectation(model: Model) -> pd.DataFrame:
    """Closed-form prediction for each genotype, from parameters alone.

    ``N(t) = N0 exp((r - d) t)`` in the crowding-free regime. Computed before any
    simulation is read, so the comparison is a genuine check and not a fit.
    """
    n0 = model.seeding["n_cells_per_population"]
    T = model.time["max_time"]
    rows = []
    for key, g in model.genotypes.items():
        net = g.cycle.net_growth_rate_per_min
        rows.append(
            {
                "genotype": key,
                "cycle_rate_per_min": g.cycle.cycle_rate_per_min,
                "apoptosis_rate_per_min": g.apoptosis_rate_per_min,
                "expected_log_growth_rate_per_min": net,
                "expected_doubling_time_days": g.doubling_time_days,
                "expected_n_end": n0 * np.exp(net * T),
            }
        )
    return pd.DataFrame(rows)


def expectation_check(
    reps: pd.DataFrame, model: Model, tolerance: float = 0.15
) -> pd.DataFrame:
    """Compare monoculture arms against the closed form.

    Reviewer recommended issue #2: "you should simulate each parameterized cell
    type by itself to confirm that its growth matches what is mathematically
    expected from its individual settings." A monoculture that misses the closed
    form by more than ``tolerance`` means the config is wrong, or crowding is
    already biting, and either way the mixed-arm results cannot be interpreted
    until it is explained.
    """
    exp = analytic_expectation(model).set_index("genotype")
    mono = reps[reps.arm.str.startswith("mono_")]
    rows = []
    for (arm, pop), g in mono.groupby(["arm", "population"]):
        if pop not in exp.index:
            continue
        e = exp.loc[pop]
        obs_rate = g[PRIMARY_OUTCOME].mean()
        exp_rate = e["expected_log_growth_rate_per_min"]
        rel = (obs_rate - exp_rate) / exp_rate if exp_rate else np.nan
        rows.append(
            {
                "arm": arm,
                "population": pop,
                "n_replicates": len(g),
                "expected_log_growth_rate_per_min": exp_rate,
                "observed_log_growth_rate_per_min": obs_rate,
                "relative_error": rel,
                "expected_n_end": e["expected_n_end"],
                "observed_n_end_mean": g["n_end"].mean(),
                "passes": bool(abs(rel) <= tolerance) if np.isfinite(rel) else False,
            }
        )
    return pd.DataFrame(rows)


def competition_coefficients(reps: pd.DataFrame) -> pd.DataFrame:
    """Growth in co-culture minus growth in monoculture, paired within seed.

    Reviewer major issue #5: a change in relative abundance shows the assigned
    birth and death schedules and "cannot show competition between the cells".
    The only thing that can is a growth rate that differs from the same
    genotype's growth rate alone, at the same seed. If this coefficient is zero,
    there is no competition, and the manuscript must drop the word.
    """
    mono = reps[reps.arm.str.startswith("mono_")][["arm", "population"]].drop_duplicates()
    mono_map = dict(zip(mono["population"], mono["arm"]))
    rows = []
    co = reps[~reps.arm.str.startswith("mono_")]
    for _, r in co.iterrows():
        mono_arm = mono_map.get(r["population"])
        if mono_arm is None:
            continue
        m = reps[
            (reps.arm == mono_arm)
            & (reps.population == r["population"])
            & (reps.seed == r["seed"])
        ]
        if m.empty:
            continue
        rows.append(
            {
                "arm": r["arm"],
                "seed": r["seed"],
                "population": r["population"],
                "rate_coculture": r[PRIMARY_OUTCOME],
                "rate_monoculture": m.iloc[0][PRIMARY_OUTCOME],
                "competition_coefficient": r[PRIMARY_OUTCOME] - m.iloc[0][PRIMARY_OUTCOME],
            }
        )
    return pd.DataFrame(rows)
