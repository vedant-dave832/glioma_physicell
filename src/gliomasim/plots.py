"""Manuscript figures.

Reviewer major issue #8: the Results section references Figures 1-8 that are not
in the PDF. These are generated from the analysis outputs so the figures cannot
drift from the numbers in the text.

Every figure writes a companion CSV of the exact values plotted. That is the
table view: no value in a figure is reachable only by reading pixels off it, and
a reviewer can check any point without rerunning anything.

Design: light print surface, hairline solid grid, thin marks, legend always
present with selective direct labels at line ends. Series colours come from
model.yaml and are fixed per entity, so a population keeps its colour across
every figure in the manuscript.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .registry import Model  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FALLBACK = ["#2a78d6", "#eb6834", "#1baf7a"]

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK_2,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlecolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _colours(model: Model, populations: list[str]) -> dict[str, str]:
    out = {}
    for i, p in enumerate(sorted(populations)):
        g = model.genotypes.get(p)
        out[p] = g.colour if g else FALLBACK[i % len(FALLBACK)]
    return out


def _label(model: Model, pop: str) -> str:
    g = model.genotypes.get(pop)
    return g.display_name if g else pop


def _finish(ax, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=10)
    ax.grid(True, axis="y", zorder=0)
    ax.set_axisbelow(True)


def _save(fig, out_dir: Path, name: str, table: pd.DataFrame) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{name}.png"
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(out_dir / f"{name}.pdf")
    table.to_csv(out_dir / f"{name}_data.csv", index=False)
    plt.close(fig)
    return png


def trajectories(
    traj: pd.DataFrame, model: Model, out_dir: Path, arm: str, name: str | None = None
) -> Path:
    """Mean live count over time with a bootstrap CI band, one panel per arm.

    Log y-axis: the model is exponential by construction, so a log axis makes a
    constant growth rate a straight line and any departure from it visible.
    """
    g = traj[traj.arm == arm]
    pops = sorted(g["population"].unique())
    colours = _colours(model, pops)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    rows = []
    for pop in pops:
        sub = g[g.population == pop]
        agg = (
            sub.groupby("time_days")["n_live"]
            .agg(["mean", "std", "count", "min", "max"])
            .reset_index()
        )
        sem = agg["std"] / np.sqrt(agg["count"].clip(lower=1))
        lo, hi = agg["mean"] - 1.96 * sem, agg["mean"] + 1.96 * sem
        ax.fill_between(
            agg["time_days"], lo.clip(lower=1), hi, color=colours[pop], alpha=0.16, lw=0
        )
        ax.plot(
            agg["time_days"], agg["mean"], color=colours[pop], lw=2, label=_label(model, pop)
        )
        # selective direct label: the endpoint only
        ax.annotate(
            f"{agg['mean'].iloc[-1]:.0f}",
            (agg["time_days"].iloc[-1], agg["mean"].iloc[-1]),
            textcoords="offset points",
            xytext=(6, 0),
            color=colours[pop],
            fontsize=8,
            fontweight="bold",
            va="center",
        )
        agg.insert(0, "population", pop)
        agg["ci_low"], agg["ci_high"] = lo, hi
        rows.append(agg)

    # Log axis only when the data actually spans enough for it to help. Over a
    # short run the counts barely move, and a log axis then produces ticks like
    # "1.15 x 10^2" that are harder to read than plain numbers.
    all_means = pd.concat(rows)["mean"]
    span = all_means.max() / max(all_means.min(), 1e-9)
    log_scale = span >= 5
    if log_scale:
        ax.set_yscale("log")
    _finish(
        ax,
        "Time (days)",
        "Live agents" + (" (log scale)" if log_scale else ""),
        f"Population growth — {arm}",
    )
    n_reps = g["seed"].nunique()
    ax.legend(loc="upper left", fontsize=8)
    ax.text(
        0.99,
        0.02,
        f"mean of {n_reps} seeded replicates, 95% CI",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color=MUTED,
        fontsize=7.5,
    )
    ax.margins(x=0.12)
    return _save(fig, out_dir, name or f"fig_trajectories_{arm}", pd.concat(rows))


def effect_sizes(contrasts: pd.DataFrame, out_dir: Path, floor: dict | None = None) -> Path:
    """Effect size with CI for each paired contrast, against the noise floor.

    The shaded band is the label control: the effect the pipeline produces
    between two identically-parameterised populations. An interval that overlaps
    it has not shown a genotype effect.
    """
    d = contrasts.dropna(subset=["mean_difference"]).copy()
    d["y"] = np.arange(len(d))[::-1]

    fig, ax = plt.subplots(figsize=(6.4, 0.62 * max(len(d), 3) + 1.5))

    if floor and np.isfinite(floor.get("ci_low", np.nan)):
        ax.axvspan(
            floor["ci_low"],
            floor["ci_high"],
            color=MUTED,
            alpha=0.16,
            lw=0,
            label="label-control noise floor (95% CI)",
        )
    ax.axvline(0, color=AXIS, lw=1)

    ax.errorbar(
        d["mean_difference"],
        d["y"],
        xerr=[
            d["mean_difference"] - d["ci_low"],
            d["ci_high"] - d["mean_difference"],
        ],
        fmt="o",
        color=FALLBACK[0],
        ecolor=FALLBACK[0],
        elinewidth=2,
        capsize=0,
        markersize=8,
        markeredgecolor=SURFACE,
        markeredgewidth=2,
    )
    ax.set_yticks(d["y"])
    ax.set_yticklabels([f"{r.arm}\n{r.comparison}" for r in d.itertuples()], fontsize=8)
    _finish(ax, "Difference in log growth rate (1/min)", "", "Paired effect sizes by arm")
    ax.grid(True, axis="x")
    ax.grid(False, axis="y")
    ax.margins(y=0.22)
    # Legend below the axis: inside the plot it lands on top of the noise-floor
    # band, which is the one thing the reader needs to see unobstructed.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8, ncol=2)
    return _save(fig, out_dir, "fig_effect_sizes", d)


def sensitivity(
    sweep: pd.DataFrame, out_dir: Path, model: Model | None = None
) -> Path:
    """Outcome as a function of the assumed apoptosis suppression factor.

    Reviewer recommended issue #2. The x-axis is the assumption; the point of the
    figure is that the headline result is a function of it, not a finding.

    Series colours come from the model so a population keeps the same colour it
    has in every other figure. Colour follows the entity, never the plotting
    order -- a reader who learned "TP53-mutant is orange" must not meet a blue
    one here.
    """
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    pops = sorted(sweep["population"].unique())
    colours = _colours(model, pops) if model else {
        p: FALLBACK[i % len(FALLBACK)] for i, p in enumerate(pops)
    }

    for pop in pops:
        sub = sweep[sweep.population == pop].sort_values("apoptosis_factor")
        c = colours[pop]
        sem = sub["sd"] / np.sqrt(sub["n_replicates"].clip(lower=1))
        ax.errorbar(
            sub["apoptosis_factor"], sub["mean"], yerr=1.96 * sem,
            fmt="-o", color=c, ecolor=c, lw=2, ms=6, elinewidth=1.4, capsize=3,
            markeredgecolor=SURFACE, markeredgewidth=2,
            label=_label(model, pop) if model else pop,
        )

    # The ceiling: with apoptosis removed entirely, growth is capped by the cycle
    # rate. It is what makes the curve plateau, so drawing it turns an unexplained
    # flattening into a stated bound.
    if model is not None:
        g = model.genotypes.get("IDHmut_TP53mut")
        if g is not None:
            n0 = model.seeding["n_cells_per_population"]
            ceiling = n0 * np.exp(g.cycle.cycle_rate_per_min * model.time["max_time"])
            ax.axhline(ceiling, color=MUTED, lw=1.2, ls=(0, (4, 3)))
            ax.annotate(f"ceiling with apoptosis removed entirely ({ceiling:.0f})",
                        (sweep["apoptosis_factor"].max(), ceiling),
                        textcoords="offset points", xytext=(0, 5), ha="right",
                        color=MUTED, fontsize=7.5)

    ax.axvline(1.0, color=AXIS, lw=1)
    ax.annotate("1x = no effect\n(the null)", (1.0, 0.02), xycoords=("data", "axes fraction"),
                textcoords="offset points", xytext=(5, 0), color=MUTED, fontsize=7.5,
                va="bottom")
    ax.set_xscale("log")
    _finish(
        ax,
        "Assumed apoptosis suppression factor (x-fold reduction)",
        "Mean endpoint count",
        "Sensitivity to the unsourced assumption",
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.margins(x=0.18, y=0.16)
    return _save(fig, out_dir, "fig_sensitivity", sweep)


def convergence(conv: pd.DataFrame, out_dir: Path, target: float = 0.10) -> Path:
    """Relative CI half-width against replicate count, with the target line."""
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(
        conv["n_replicates"],
        conv["relative_half_width"],
        "-o",
        color=FALLBACK[0],
        lw=2,
        ms=6,
        markeredgecolor=SURFACE,
        markeredgewidth=2,
        label="relative 95% CI half-width",
    )
    ax.axhline(target, color=FALLBACK[1], lw=2, label=f"target ({target:.0%})")
    _finish(
        ax,
        "Number of replicates",
        "CI half-width / mean",
        "Replicate-count convergence",
    )
    ax.legend(loc="upper right", fontsize=8)
    return _save(fig, out_dir, "fig_convergence", conv)


def expectation(check: pd.DataFrame, out_dir: Path) -> Path:
    """Observed versus closed-form growth rate for each monoculture arm."""
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    x = np.arange(len(check))
    w = 0.36
    ax.bar(x - w / 2, check["expected_log_growth_rate_per_min"], w,
           color=FALLBACK[0], label="expected (closed form)")
    ax.bar(x + w / 2, check["observed_log_growth_rate_per_min"], w,
           color=FALLBACK[1], label="observed (simulation)")
    ax.set_xticks(x)
    ax.set_xticklabels(check["population"], fontsize=8)
    for xi, r in zip(x, check.itertuples()):
        ax.annotate(f"{r.relative_error:+.1%}", (xi, max(
            r.expected_log_growth_rate_per_min, r.observed_log_growth_rate_per_min)),
            textcoords="offset points", xytext=(0, 4), ha="center",
            color=INK_2, fontsize=8)
    _finish(ax, "", "Log growth rate (1/min)", "Monoculture sanity check")
    ax.legend(loc="upper right", fontsize=8)
    return _save(fig, out_dir, "fig_expectation_check", check)
