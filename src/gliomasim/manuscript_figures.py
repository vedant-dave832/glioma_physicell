"""The eight manuscript figures, rebuilt from the replicate run set.

Reviewer major issue #8: Figures 1-8 are referenced in the Results but absent
from the PDF. This module regenerates that figure set.

It is NOT a literal reproduction of the original legends. Several of them
describe results the new runs contradict, and reproducing them would mean
drawing a picture of something that did not happen:

  Fig 1  "TP53-mutant fit was limited by the Day 3.25 artifact (R2 = 0.132)"
         The artifact was a failed run. It is rerun here, not worked around, so
         there is no artifact to report and no degraded fit.
  Fig 3  "IDH1-mutant competitive exclusion"
         Competition coefficients all span zero. Nothing is excluded by
         competition. The panel now shows population fraction without the claim.
  Fig 4  "Pairwise p-value heatmap from Mann-Whitney U tests" over timepoints
         Timepoints within a run are not independent. Replaced with
         replicate-level distributions and paired effect sizes.
  Fig 6  "the hypoxic threshold of 5 mmHg"
         5 mmHg is PhysiCell's NECROSIS threshold; the hypoxic threshold is
         15 mmHg. Both are drawn, labelled correctly.
  Fig 7  "All voxels remain dark green (34-38 mmHg) with no spatial hypoxic
         core" - measured minimum across the run set is 10.8 mmHg, below the
         15 mmHg hypoxic threshold. The field is plotted on its true range.

Every figure writes a companion CSV of the values plotted.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from .collect import _current_time, _labels  # noqa: E402
from .plots import AXIS, GRID, INK, INK_2, MUTED, SURFACE, _save  # noqa: E402
from .registry import Model  # noqa: E402

# Sequential ramp, one hue light->dark, from the validated palette.
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGE_RAMP = ["#fde3d6", "#f9bd9d", "#f49468", "#eb6834", "#c14f22", "#933a16", "#63250c"]
CMAP_BLUE = LinearSegmentedColormap.from_list("brand_blue", BLUE_RAMP)
CMAP_ORANGE = LinearSegmentedColormap.from_list("brand_orange", ORANGE_RAMP)

# PhysiCell oxygen thresholds (defaults), in mmHg.
HYPOXIC_MMHG = 15.0
NECROSIS_MMHG = 5.0

_OUT_RE = re.compile(r"output(\d+)\.xml$")


# ---------------------------------------------------------------- readers ----
def frame_paths(run_dir: Path) -> list[Path]:
    return sorted(
        p for p in (run_dir / "output").glob("output*.xml") if _OUT_RE.search(p.name)
    )


def read_positions(xml_path: Path, type_names: dict[int, str]) -> pd.DataFrame:
    """Live agent positions for one frame."""
    mat = xml_path.with_name(xml_path.stem + "_cells.mat")
    cells = loadmat(mat)["cells"]
    lab = _labels(xml_path)
    ci, di, pi = lab["cell_type"][0], lab["dead"][0], lab["position"][0]
    live = cells[di, :].astype(int) == 0
    return pd.DataFrame(
        {
            "x": cells[pi, live],
            "y": cells[pi + 1, live],
            "population": [
                type_names.get(int(t), str(int(t))) for t in cells[ci, live]
            ],
            "time_days": _current_time(xml_path) / 1440.0,
        }
    )


def read_oxygen_field(xml_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Oxygen on the mesh for one frame, as (x, y, value) columns."""
    mat = xml_path.with_name(xml_path.stem + "_microenvironment0.mat")
    arr = loadmat(mat)["multiscale_microenvironment"]
    return arr[0, :], arr[1, :], arr[4, :]


def _nearest_frames(paths: list[Path], targets_days: list[float]) -> list[Path]:
    times = np.array([_current_time(p) / 1440.0 for p in paths])
    return [paths[int(np.argmin(np.abs(times - t)))] for t in targets_days]


# ------------------------------------------------------------------ style ----
def _panel(ax, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=9.5, pad=8)
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)


def _colours(model: Model, pops: list[str]) -> dict[str, str]:
    fallback = ["#2a78d6", "#eb6834", "#1baf7a"]
    out = {}
    for i, p in enumerate(sorted(pops)):
        g = model.genotypes.get(p)
        out[p] = g.colour if g else fallback[i % len(fallback)]
    return out


def _label(model: Model, pop: str) -> str:
    g = model.genotypes.get(pop)
    return g.display_name if g else pop


# ---------------------------------------------------------------- figures ----
def fig1_growth_and_fold_change(traj, reps, model, out_dir, arm="co_separated"):
    """Growth curves with replicate CIs, and total fold change."""
    g = traj[traj.arm == arm]
    pops = sorted(g["population"].unique())
    col = _colours(model, pops)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.7))
    rows = []
    for pop in pops:
        sub = g[g.population == pop]
        agg = sub.groupby("time_days")["n_live"].agg(["mean", "std", "count"]).reset_index()
        sem = agg["std"] / np.sqrt(agg["count"].clip(lower=1))
        ax1.fill_between(agg["time_days"], agg["mean"] - 1.96 * sem,
                         agg["mean"] + 1.96 * sem, color=col[pop], alpha=0.16, lw=0)
        ax1.plot(agg["time_days"], agg["mean"], color=col[pop], lw=2,
                 label=_label(model, pop))
        agg.insert(0, "population", pop)
        rows.append(agg)
    _panel(ax1, "Time (days)", "Live agents", "A. Population growth")
    ax1.legend(fontsize=8, loc="upper left")

    r = reps[reps.arm == arm]
    fold = r.groupby("population")["fold_change"].agg(["mean", "std", "count"])
    x = np.arange(len(fold))
    sem = fold["std"] / np.sqrt(fold["count"].clip(lower=1))
    ax2.bar(x, fold["mean"], 0.55, yerr=1.96 * sem, capsize=4,
            color=[col[p] for p in fold.index],
            error_kw={"ecolor": INK_2, "elinewidth": 1.2})
    # Label above the error bar, not the bar top, or the whisker crosses the text.
    for xi, v, e in zip(x, fold["mean"], 1.96 * sem):
        ax2.annotate(f"{v:.2f}x", (xi, v + (e if np.isfinite(e) else 0)),
                     textcoords="offset points", xytext=(0, 6), ha="center",
                     fontsize=8.5, color=INK_2, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels([_label(model, p).replace(" / ", "\n") for p in fold.index],
                        fontsize=8)
    _panel(ax2, "", "Fold change from seeding", "B. Total expansion")
    ax2.margins(y=0.18)

    table = pd.concat(rows).merge(
        fold.reset_index().rename(columns={"mean": "fold_change_mean"}),
        on="population", how="left")
    return _save(fig, out_dir, "fig1_growth_and_fold_change", table)


def fig2_spatial_agents(run_dir, model, out_dir, days=(0.0, 3.3, 5.5)):
    """Agent positions at three timepoints, one representative replicate."""
    import json

    manifest = json.loads((run_dir / "manifest.json").read_text())
    type_names = {p["type_id"]: p["label"] for p in manifest["populations"]}
    frames = _nearest_frames(frame_paths(run_dir), list(days))

    dfs = [read_positions(f, type_names) for f in frames]
    pops = sorted({p for d in dfs for p in d["population"].unique()})
    col = _colours(model, pops)

    fig, axes = plt.subplots(1, len(frames), figsize=(3.15 * len(frames), 3.6),
                             sharex=True, sharey=True)
    for ax, d in zip(np.atleast_1d(axes), dfs):
        for pop in pops:
            s = d[d.population == pop]
            ax.scatter(s["x"], s["y"], s=7, color=col[pop], alpha=0.85,
                       linewidths=0, label=_label(model, pop))
        ax.set_title(f"Day {d['time_days'].iloc[0]:.1f}  (n = {len(d)})",
                     loc="left", fontsize=9, pad=6)
        ax.set_aspect("equal")
        ax.set_xlim(model.domain["x_min"], model.domain["x_max"])
        ax.set_ylim(model.domain["y_min"], model.domain["y_max"])
        ax.set_xlabel("x (µm)")
        ax.grid(True, lw=0.5)
        ax.set_axisbelow(True)
    np.atleast_1d(axes)[0].set_ylabel("y (µm)")
    # Legend inside the first panel: a figure-level legend below the axes
    # collides with the x-axis labels once tight_layout runs.
    np.atleast_1d(axes)[0].legend(fontsize=7.5, loc="upper left", markerscale=1.6)
    fig.suptitle(f"Spatial distribution of agents — {manifest['arm']}, "
                 f"seed {manifest['seed']}", x=0.02, ha="left",
                 fontsize=10, fontweight="bold")
    return _save(fig, out_dir, "fig2_spatial_agents", pd.concat(dfs))


def fig3_fraction_and_doubling(traj, reps, model, out_dir, arm="co_separated"):
    """Population fraction over time, and observed vs expected doubling time."""
    g = traj[traj.arm == arm].copy()
    tot = g.groupby("time_days")["n_live"].transform("sum")
    g["fraction"] = g["n_live"] / tot
    pops = sorted(g["population"].unique())
    col = _colours(model, pops)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.7))
    rows = []
    for pop in pops:
        agg = g[g.population == pop].groupby("time_days")["fraction"].agg(
            ["mean", "std", "count"]).reset_index()
        sem = agg["std"] / np.sqrt(agg["count"].clip(lower=1))
        ax1.fill_between(agg["time_days"], agg["mean"] - 1.96 * sem,
                         agg["mean"] + 1.96 * sem, color=col[pop], alpha=0.16, lw=0)
        ax1.plot(agg["time_days"], agg["mean"], color=col[pop], lw=2,
                 label=_label(model, pop))
        agg.insert(0, "population", pop)
        rows.append(agg)
    ax1.axhline(0.5, color=AXIS, lw=1)
    ax1.set_ylim(0, 1)
    _panel(ax1, "Time (days)", "Fraction of live agents", "A. Population fraction")
    ax1.legend(fontsize=8, loc="center left")

    obs = reps[reps.arm == arm].groupby("population")[
        "observed_doubling_time_days"].agg(["mean", "std", "count"])
    exp = {k: v.doubling_time_days for k, v in model.genotypes.items()}
    x = np.arange(len(obs))
    w = 0.36
    sem = obs["std"] / np.sqrt(obs["count"].clip(lower=1))
    ax2.bar(x - w / 2, [exp.get(p, np.nan) for p in obs.index], w,
            color=MUTED, label="implied by parameters")
    ax2.bar(x + w / 2, obs["mean"], w, yerr=1.96 * sem, capsize=4,
            color=[col[p] for p in obs.index], label="observed in simulation",
            error_kw={"ecolor": INK_2, "elinewidth": 1.2})
    ax2.set_xticks(x)
    ax2.set_xticklabels([_label(model, p).replace(" / ", "\n") for p in obs.index],
                        fontsize=8)
    _panel(ax2, "", "Population doubling time (days)",
           "B. Observed vs implied doubling time")
    ax2.legend(fontsize=8)
    ax2.margins(y=0.18)

    table = pd.concat(rows).merge(
        obs.reset_index().rename(columns={"mean": "observed_doubling_days"}),
        on="population", how="left")
    return _save(fig, out_dir, "fig3_fraction_and_doubling", table)


def fig4_replicate_distributions(reps, contrasts, floor, model, out_dir):
    """Replicate-level distributions, and paired effect sizes vs the noise floor.

    This replaces the original Figure 4. Its notched boxplots over timepoints and
    Mann-Whitney heatmap treated dependent timepoints as independent samples.
    """
    arms = [a for a in sorted(reps["arm"].unique()) if not a.startswith("mono_")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))

    groups, labels, colours = [], [], []
    for arm in arms:
        for pop in sorted(reps[reps.arm == arm]["population"].unique()):
            v = reps[(reps.arm == arm) & (reps.population == pop)]["n_end"]
            groups.append(v.to_numpy())
            labels.append(f"{arm}\n{pop}")
            colours.append(_colours(model, [pop])[pop])

    bp = ax1.boxplot(groups, patch_artist=True, widths=0.55, showfliers=False,
                     medianprops={"color": INK, "lw": 1.4},
                     whiskerprops={"color": AXIS}, capprops={"color": AXIS},
                     boxprops={"edgecolor": SURFACE, "lw": 2})
    for patch, c in zip(bp["boxes"], colours):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    for i, v in enumerate(groups, 1):
        ax1.scatter(np.full(len(v), i) + np.random.uniform(-0.09, 0.09, len(v)),
                    v, s=14, color=INK_2, alpha=0.65, linewidths=0, zorder=3)
    ax1.set_xticks(range(1, len(labels) + 1))
    ax1.set_xticklabels(labels, fontsize=6.5, rotation=30, ha="right")
    _panel(ax1, "", "Endpoint live agents",
           "A. Endpoint count per replicate (each dot is one run)")

    d = contrasts.dropna(subset=["mean_difference"]).copy()
    y = np.arange(len(d))[::-1]
    if floor and np.isfinite(floor.get("ci_low", np.nan)):
        ax2.axvspan(floor["ci_low"], floor["ci_high"], color=MUTED, alpha=0.16,
                    lw=0, label="label-control noise floor")
    ax2.axvline(0, color=AXIS, lw=1)
    ax2.errorbar(d["mean_difference"], y,
                 xerr=[d["mean_difference"] - d["ci_low"],
                       d["ci_high"] - d["mean_difference"]],
                 fmt="o", color="#2a78d6", ecolor="#2a78d6", elinewidth=2,
                 capsize=0, markersize=7, markeredgecolor=SURFACE, markeredgewidth=2)
    ax2.set_yticks(y)
    ax2.set_yticklabels(d["arm"], fontsize=8)
    _panel(ax2, "Difference in log growth rate (1/min)", "",
           "B. Paired effect size, within seed")
    ax2.grid(True, axis="x")
    ax2.grid(False, axis="y")
    ax2.margins(y=0.25)
    ax2.legend(fontsize=8, loc="lower right")

    table = reps[["arm", "seed", "population", "n_end", "log_growth_rate_per_min"]]
    return _save(fig, out_dir, "fig4_replicate_distributions", table)


def fig5_diversity_and_counts(traj, arm_reps, model, out_dir, arm="co_separated"):
    """Shannon diversity over time, and absolute counts."""
    g = traj[traj.arm == arm]
    pops = sorted(g["population"].unique())
    col = _colours(model, pops)

    div = []
    for (t,), sub in g.groupby(["time_days"]):
        per_seed = []
        for _, s in sub.groupby("seed"):
            n = s["n_live"].to_numpy(float)
            tot = n.sum()
            p = n[n > 0] / tot if tot else np.array([])
            per_seed.append(float(-(p * np.log(p)).sum()) if len(p) else np.nan)
        div.append({"time_days": t, "mean": np.nanmean(per_seed),
                    "sem": np.nanstd(per_seed, ddof=1) / np.sqrt(len(per_seed))})
    div = pd.DataFrame(div).sort_values("time_days")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.7))
    ax1.fill_between(div["time_days"], div["mean"] - 1.96 * div["sem"],
                     div["mean"] + 1.96 * div["sem"], color="#2a78d6",
                     alpha=0.16, lw=0)
    ax1.plot(div["time_days"], div["mean"], color="#2a78d6", lw=2)
    ax1.axhline(np.log(len(pops)), color=AXIS, lw=1)
    ax1.annotate(f"maximum for {len(pops)} equal populations (ln {len(pops)})",
                 (div["time_days"].iloc[0], np.log(len(pops))),
                 textcoords="offset points", xytext=(4, -12),
                 fontsize=7.5, color=MUTED)
    _panel(ax1, "Time (days)", "Shannon index H", "A. Diversity over time")
    ax1.annotate(f"{div['mean'].iloc[0]:.3f} → {div['mean'].iloc[-1]:.3f} "
                 f"({div['mean'].iloc[-1] - div['mean'].iloc[0]:+.3f})",
                 (0.98, 0.06), xycoords="axes fraction", ha="right",
                 fontsize=8, color=INK_2)

    rows = []
    for pop in pops:
        agg = g[g.population == pop].groupby("time_days")["n_live"].agg(
            ["mean", "std", "count"]).reset_index()
        sem = agg["std"] / np.sqrt(agg["count"].clip(lower=1))
        ax2.fill_between(agg["time_days"], agg["mean"] - 1.96 * sem,
                         agg["mean"] + 1.96 * sem, color=col[pop], alpha=0.16, lw=0)
        ax2.plot(agg["time_days"], agg["mean"], color=col[pop], lw=2,
                 label=_label(model, pop))
        agg.insert(0, "population", pop)
        rows.append(agg)
    _panel(ax2, "Time (days)", "Live agents", "B. Absolute counts")
    ax2.legend(fontsize=8, loc="upper left")

    table = div.assign(quantity="shannon")
    return _save(fig, out_dir, "fig5_diversity_and_counts",
                 pd.concat([table] + rows, ignore_index=True))


def fig6_oxygen_over_time(run_dirs, out_dir):
    """Mean oxygen and the min-max envelope over time, with correct thresholds."""
    rows = []
    for run_dir in run_dirs:
        for f in frame_paths(run_dir):
            _, _, o2 = read_oxygen_field(f)
            rows.append({"run": run_dir.parent.name, "time_days": _current_time(f) / 1440.0,
                         "mean": float(o2.mean()), "min": float(o2.min()),
                         "max": float(o2.max())})
    df = pd.DataFrame(rows)
    agg = df.groupby("time_days").agg(mean=("mean", "mean"), lo=("min", "min"),
                                      hi=("max", "max")).reset_index()

    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    ax.fill_between(agg["time_days"], agg["lo"], agg["hi"], color="#2a78d6",
                    alpha=0.16, lw=0, label="min–max across domain and runs")
    ax.plot(agg["time_days"], agg["mean"], color="#2a78d6", lw=2, label="domain mean")
    ax.axhline(HYPOXIC_MMHG, color="#eb6834", lw=1.6, label=f"hypoxic threshold ({HYPOXIC_MMHG:.0f} mmHg)")
    ax.axhline(NECROSIS_MMHG, color="#933a16", lw=1.6, label=f"necrosis threshold ({NECROSIS_MMHG:.0f} mmHg)")
    _panel(ax, "Time (days)", "Oxygen (mmHg)", "Oxygen over the simulation")
    ax.set_ylim(0, 40)
    ax.legend(fontsize=7.5, loc="lower left")
    ax.annotate(f"minimum reached: {agg['lo'].min():.1f} mmHg",
                (0.98, 0.55), xycoords="axes fraction", ha="right",
                fontsize=8, color=INK_2)
    return _save(fig, out_dir, "fig6_oxygen_over_time", agg)


def fig7_spatial_oxygen(run_dir, model, out_dir, days=(0.0, 3.3, 5.5)):
    """Oxygen field at three timepoints, on its true measured range."""
    frames = _nearest_frames(frame_paths(run_dir), list(days))
    fields = [read_oxygen_field(f) for f in frames]
    vmin = min(v.min() for _, _, v in fields)
    vmax = max(v.max() for _, _, v in fields)

    fig, axes = plt.subplots(1, len(frames), figsize=(3.15 * len(frames), 3.5),
                             sharey=True)
    rows = []
    for ax, f, (x, y, v) in zip(np.atleast_1d(axes), frames, fields):
        n = int(np.sqrt(len(v)))
        im = ax.imshow(v.reshape(n, n), origin="lower", cmap=CMAP_BLUE,
                       vmin=vmin, vmax=vmax,
                       extent=[model.domain["x_min"], model.domain["x_max"],
                               model.domain["y_min"], model.domain["y_max"]])
        t = _current_time(f) / 1440.0
        ax.set_title(f"Day {t:.1f}   min {v.min():.1f} mmHg", loc="left", fontsize=9)
        ax.set_xlabel("x (µm)")
        rows.append(pd.DataFrame({"time_days": t, "x": x, "y": y, "oxygen_mmHg": v}))
    np.atleast_1d(axes)[0].set_ylabel("y (µm)")
    cb = fig.colorbar(im, ax=np.atleast_1d(axes).tolist(), fraction=0.025, pad=0.05)
    cb.set_label("Oxygen (mmHg)", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("Spatial oxygen distribution", x=0.02, ha="left",
                 fontsize=10, fontweight="bold")
    fig.subplots_adjust(top=0.86)
    png = out_dir / "fig7_spatial_oxygen.png"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "fig7_spatial_oxygen.pdf", bbox_inches="tight")
    pd.concat(rows).to_csv(out_dir / "fig7_spatial_oxygen_data.csv", index=False)
    plt.close(fig)
    return png


def fig8_density_kde(run_dir, model, out_dir, day=5.5):
    """Kernel density of each population at the endpoint."""
    import json

    manifest = json.loads((run_dir / "manifest.json").read_text())
    type_names = {p["type_id"]: p["label"] for p in manifest["populations"]}
    frame = _nearest_frames(frame_paths(run_dir), [day])[0]
    pos = read_positions(frame, type_names)
    pops = sorted(pos["population"].unique())
    cmaps = [CMAP_BLUE, CMAP_ORANGE]

    xmin, xmax = model.domain["x_min"], model.domain["x_max"]
    ymin, ymax = model.domain["y_min"], model.domain["y_max"]
    gx, gy = np.mgrid[xmin:xmax:120j, ymin:ymax:120j]
    grid = np.vstack([gx.ravel(), gy.ravel()])

    fig, axes = plt.subplots(1, len(pops), figsize=(4.1 * len(pops), 3.8), sharey=True)
    rows = []
    for ax, pop, cm in zip(np.atleast_1d(axes), pops, cmaps):
        s = pos[pos.population == pop]
        try:
            z = gaussian_kde(np.vstack([s["x"], s["y"]]))(grid).reshape(gx.shape)
        except Exception:  # noqa: BLE001 - degenerate (too few or collinear points)
            z = np.zeros_like(gx)
        ax.imshow(z.T, origin="lower", cmap=cm,
                  extent=[xmin, xmax, ymin, ymax])
        ax.scatter(s["x"], s["y"], s=3, color=SURFACE, alpha=0.5, linewidths=0)
        ax.set_title(f"{_label(model, pop)}\nn = {len(s)}", loc="left", fontsize=9)
        ax.set_xlabel("x (µm)")
        ax.set_aspect("equal")
        rows.append(pd.DataFrame({"population": pop, "x": s["x"], "y": s["y"]}))
    np.atleast_1d(axes)[0].set_ylabel("y (µm)")
    fig.suptitle(f"Agent density at day {day:.1f} — {manifest['arm']}, "
                 f"seed {manifest['seed']}", x=0.02, ha="left",
                 fontsize=10, fontweight="bold")
    return _save(fig, out_dir, "fig8_density_kde", pd.concat(rows))


# -------------------------------------------------------------- orchestrate --
def build_all(traj, reps, arm_reps, contrasts, floor, model, runs_root: Path,
              out_dir: Path, arm: str = "co_separated") -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    seeds = sorted(reps[reps.arm == arm]["seed"].unique())
    run_dir = runs_root / arm / f"seed_{seeds[0]}"
    mono_dirs = [
        runs_root / a / f"seed_{s}"
        for a in sorted(reps["arm"].unique())
        for s in sorted(reps[reps.arm == a]["seed"].unique())[:3]
    ]

    made = [
        fig1_growth_and_fold_change(traj, reps, model, out_dir, arm),
        fig2_spatial_agents(run_dir, model, out_dir),
        fig3_fraction_and_doubling(traj, reps, model, out_dir, arm),
        fig4_replicate_distributions(reps, contrasts, floor, model, out_dir),
        fig5_diversity_and_counts(traj, arm_reps, model, out_dir, arm),
        fig6_oxygen_over_time(mono_dirs, out_dir),
        fig7_spatial_oxygen(run_dir, model, out_dir),
        fig8_density_kde(run_dir, model, out_dir),
    ]
    return made
