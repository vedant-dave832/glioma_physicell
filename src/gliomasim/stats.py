"""Replicate-level statistics.

Rules this module enforces, all of them from reviewer major issue #4:

* The unit of analysis is a run, never a timepoint. Nothing here accepts a
  trajectory frame.
* Arms are compared within seed. Replicate i of every arm shares a seed and a
  starting arrangement, so comparisons are paired and the seed-to-seed variance
  is removed rather than counted as evidence.
* Every comparison is reported as an effect size with a confidence interval.
  A p-value never appears without one.
* The label control sets the false-positive floor. Any genotype effect smaller
  than the effect seen between two identically-parameterised populations is
  noise, and is reported as such.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as sps

RNG = np.random.default_rng(20260907)


@dataclass
class Estimate:
    label: str
    n: int
    mean: float
    ci_low: float
    ci_high: float
    method: str

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "n": self.n,
            "mean": self.mean,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "method": self.method,
        }


def bootstrap_mean(
    values: np.ndarray, n_boot: int = 10000, alpha: float = 0.05, label: str = ""
) -> Estimate:
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return Estimate(label, 0, np.nan, np.nan, np.nan, "bootstrap")
    if len(v) == 1:
        return Estimate(label, 1, float(v[0]), np.nan, np.nan, "single run, no CI")
    idx = RNG.integers(0, len(v), size=(n_boot, len(v)))
    means = v[idx].mean(axis=1)
    return Estimate(
        label,
        len(v),
        float(v.mean()),
        float(np.quantile(means, alpha / 2)),
        float(np.quantile(means, 1 - alpha / 2)),
        f"percentile bootstrap, {n_boot} resamples",
    )


def describe_arms(reps: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """Distribution of one replicate-level outcome, per arm and population."""
    rows = []
    for (arm, pop), g in reps.groupby(["arm", "population"]):
        est = bootstrap_mean(g[outcome].to_numpy(), label=f"{arm}/{pop}")
        rows.append(
            {
                "arm": arm,
                "population": pop,
                "outcome": outcome,
                "n_replicates": est.n,
                "mean": est.mean,
                "sd": float(np.nanstd(g[outcome], ddof=1)) if est.n > 1 else np.nan,
                "ci_low": est.ci_low,
                "ci_high": est.ci_high,
            }
        )
    return pd.DataFrame(rows)


def paired_contrast(
    reps: pd.DataFrame, arm: str, pop_a: str, pop_b: str, outcome: str
) -> dict:
    """Within-seed paired contrast of two populations in the same arm.

    Uses the paired difference per seed. Reports the effect size and CI first;
    the p-value is included but is not the headline, because with a simulator you
    can always buy a smaller p-value by running more replicates, and that does
    not make an effect larger or more meaningful.
    """
    g = reps[reps.arm == arm]
    a = g[g.population == pop_a].set_index("seed")[outcome]
    b = g[g.population == pop_b].set_index("seed")[outcome]
    common = sorted(set(a.index) & set(b.index))
    if len(common) < 2:
        return {
            "arm": arm,
            "comparison": f"{pop_a} - {pop_b}",
            "outcome": outcome,
            "n_pairs": len(common),
            "note": "fewer than 2 paired seeds; no inference attempted",
        }
    d = (a.loc[common] - b.loc[common]).to_numpy(float)
    est = bootstrap_mean(d, label=f"{pop_a}-{pop_b}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_res = sps.ttest_rel(a.loc[common], b.loc[common])
        try:
            w_res = sps.wilcoxon(d)
            w_p = float(w_res.pvalue)
        except ValueError:
            w_p = float("nan")
    sd = float(np.std(d, ddof=1))
    return {
        "arm": arm,
        "comparison": f"{pop_a} - {pop_b}",
        "outcome": outcome,
        "n_pairs": len(common),
        "mean_difference": est.mean,
        "ci_low": est.ci_low,
        "ci_high": est.ci_high,
        "cohens_dz": est.mean / sd if sd else np.nan,
        "p_paired_t": float(t_res.pvalue),
        "p_wilcoxon_signed_rank": w_p,
    }


def mixed_effects(traj_reps: pd.DataFrame, outcome: str) -> dict:
    """Mixed-effects model with seed as a random intercept, across arms.

    Optional: falls back cleanly if statsmodels is not installed, because the
    paired contrasts above are the primary analysis and this is corroboration.
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        return {"available": False, "reason": "statsmodels not installed"}

    df = traj_reps.dropna(subset=[outcome]).copy()
    if df["population"].nunique() < 2 or df["seed"].nunique() < 3:
        return {"available": False, "reason": "not enough populations or seeds"}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.mixedlm(
                f"{outcome} ~ C(population)", df, groups=df["seed"]
            ).fit()
        return {
            "available": True,
            "formula": f"{outcome} ~ C(population) + (1|seed)",
            "n_obs": int(model.nobs),
            "n_groups": int(df["seed"].nunique()),
            "params": {k: float(v) for k, v in model.params.items()},
            "pvalues": {k: float(v) for k, v in model.pvalues.items()},
            "conf_int": {
                k: [float(a), float(b)]
                for k, (a, b) in model.conf_int().iterrows()
            },
        }
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "reason": f"fit failed: {exc}"}


def convergence(reps: pd.DataFrame, arm: str, population: str, outcome: str) -> pd.DataFrame:
    """How the estimate settles as replicates accumulate.

    Reviewer major issue #4: "first run a short convergence check to justify how
    many runs you use". Reports the running mean and the relative half-width of
    the 95% CI as a function of n, so the chosen n can be defended by pointing
    at the n where the half-width drops below the precision you need.
    """
    g = reps[(reps.arm == arm) & (reps.population == population)].sort_values("seed")
    v = g[outcome].to_numpy(float)
    rows = []
    for n in range(2, len(v) + 1):
        sub = v[:n]
        sd = float(np.std(sub, ddof=1))
        sem = sd / np.sqrt(n)
        half = 1.96 * sem
        mean = float(np.mean(sub))
        rows.append(
            {
                "n_replicates": n,
                "running_mean": mean,
                "running_sd": sd,
                "sem": sem,
                "ci_half_width": half,
                "relative_half_width": half / abs(mean) if mean else np.nan,
            }
        )
    return pd.DataFrame(rows)


def determinism_report(reps: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """Is the model deterministic at fixed seed, and stochastic across seeds?

    The manuscript claims one run sufficed because the setup was deterministic.
    Reviewer major issue #4 asks for that to be verified. Non-zero apoptosis is
    sampled probabilistically and the Live cycle transition is exponential, so
    the expectation is: identical within a seed, variable across seeds. This
    reports the across-seed coefficient of variation, which is the number that
    decides whether a single run could ever have been enough.
    """
    rows = []
    for (arm, pop), g in reps.groupby(["arm", "population"]):
        v = g[outcome].to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v) < 2:
            continue
        mean = float(np.mean(v))
        rows.append(
            {
                "arm": arm,
                "population": pop,
                "n_replicates": len(v),
                "mean": mean,
                "sd_across_seeds": float(np.std(v, ddof=1)),
                "cv_across_seeds": float(np.std(v, ddof=1) / abs(mean)) if mean else np.nan,
                "min": float(v.min()),
                "max": float(v.max()),
                "spread_pct_of_mean": float((v.max() - v.min()) / abs(mean) * 100)
                if mean
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def false_positive_floor(reps: pd.DataFrame, outcome: str) -> dict:
    """Effect size between two identically-parameterised populations.

    This is the label control. Whatever difference shows up here is what the
    pipeline produces from nothing. A genotype effect is only interpretable to
    the extent it exceeds this.
    """
    g = reps[reps.arm == "label_control"]
    pops = sorted(g["population"].unique())
    if len(pops) != 2:
        return {"available": False, "reason": "label_control arm not run"}
    res = paired_contrast(reps, "label_control", pops[0], pops[1], outcome)
    res["interpretation"] = (
        "Difference between two populations with identical parameters. Any "
        "genotype contrast of comparable magnitude is indistinguishable from noise."
    )
    return res
