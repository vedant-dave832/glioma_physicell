"""Command line interface.

    glioma check                 validate the model and print what it implies
    glioma table                 write Table 1, Methods, reviewer map
    glioma run                   run seeded replicates
    glioma sweep                 run the apoptosis-factor sensitivity sweep
    glioma analyze               outcomes, statistics, figures, report

``check`` runs no simulation and is the fastest way to see whether a parameter
edit does what you meant.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd

from . import (collect, manuscript_figures, outcomes, plots, registry, report,
               run as runner, stats)

REPO = Path(__file__).resolve().parents[2]
DEFAULT_PARAMS = REPO / "params"
DEFAULT_RUNS = REPO / "runs"
DEFAULT_RESULTS = REPO / "results"


def _physicell_root(args) -> Path:
    root = args.physicell or os.environ.get("PHYSICELL_ROOT")
    if not root:
        sys.exit(
            "PhysiCell location unknown. Pass --physicell /path/to/PhysiCell or set "
            "PHYSICELL_ROOT.\n"
            "  git clone --branch 1.14.2 https://github.com/MathCancer/PhysiCell.git\n"
            "  cd PhysiCell && make template && make -j"
        )
    return Path(root).expanduser().resolve()


def _load(args, factor: float | None = None) -> registry.Model:
    try:
        return registry.load(args.params, apoptosis_factor=factor)
    except registry.RegistryError as exc:
        sys.exit(str(exc))


# ---------------------------------------------------------------- check ------
def cmd_check(args) -> None:
    model = _load(args)
    print(f"Model loaded from {args.params}\n")
    print(f"  simulated time      {model.time['max_time']:.0f} min "
          f"({model.max_time_days:.2f} days)")
    print(f"  save interval       {model.time['save_interval']:.0f} min")
    print(f"  saved timepoints    {model.n_timepoints}  (computed, incl. t=0 and endpoint)")
    print(f"  mesh                {model.n_voxels} voxels at {model.domain['dx']:.0f} um\n")

    print("Derived cell settings (nothing below is typed in by hand):")
    hdr = f"  {'genotype':<22}{'r (1/min)':>13}{'d (1/min)':>13}{'phase (min)':>13}{'T_double (d)':>14}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for k, g in model.genotypes.items():
        print(f"  {k:<22}{g.cycle.cycle_rate_per_min:>13.4e}"
              f"{g.apoptosis_rate_per_min:>13.4e}"
              f"{g.cycle.mean_phase_duration_min:>13.0f}"
              f"{g.doubling_time_days:>14.3f}")

    exp = outcomes.analytic_expectation(model)
    print("\nExpected endpoint counts BY CONSTRUCTION (before any simulation):")
    for r in exp.itertuples():
        print(f"  {r.genotype:<22}{r.expected_n_end:>10.0f} agents")
    prim = [g for g in model.primary_genotypes()]
    if len(prim) == 2:
        a, b = (exp.set_index("genotype").loc[g.key, "expected_n_end"] for g in prim)
        print(f"\n  ratio {prim[1].key} / {prim[0].key} = {b / a:.3f}")
        print("  This is arithmetic from the parameters. The simulation cannot")
        print("  discover it, and the manuscript must not present it as a finding.")

    print("\nProvenance of parameters:")
    for r in report.provenance_summary(model).itertuples():
        print(f"  {r.Provenance:<12}{r.n_parameters:>4}")
    unsourced = [p.name for p in model.params
                 if p.provenance == "assumed" and p.sweep]
    if unsourced:
        print("\nSwept assumptions (results must be reported as a function of these):")
        for u in unsourced:
            print(f"  - {u}")
    print("\nOK: model is internally consistent.")


# ----------------------------------------------------------------- table -----
def cmd_table(args) -> None:
    model = _load(args)
    prov = (
        runner.physicell_provenance(_physicell_root(args))
        if (args.physicell or os.environ.get("PHYSICELL_ROOT"))
        else {"physicell_version": model.meta["physicell_version"],
              "physicell_commit": "unknown (pass --physicell to record it)"}
    )
    paths = report.write_all(model, Path(args.out), prov, args.replicates
                             or model.run["n_replicates"])
    for k, p in paths.items():
        print(f"  {k:<16}{p}")


# ------------------------------------------------------------------- run -----
def cmd_run(args) -> None:
    model = _load(args, factor=args.apoptosis_factor)
    physicell = _physicell_root(args)
    specs = runner.plan(
        model,
        Path(args.runs),
        physicell,
        arms=args.arms,
        n_replicates=args.replicates,
        apoptosis_factor=args.apoptosis_factor,
    )
    print(f"{len(specs)} runs planned "
          f"({len(set(s.arm for s in specs))} arms x "
          f"{len(set(s.seed for s in specs))} seeds)")

    def progress(i, n, spec, result):
        status = result.get("status", "?")
        wall = result.get("wall_seconds")
        extra = f"{wall:.0f}s" if wall else ""
        print(f"  [{i}/{n}] {spec.arm}/seed_{spec.seed:<12} {status:<10} {extra}",
              flush=True)

    results = runner.run_all(specs, model, physicell, force=args.force,
                             on_progress=progress, max_seconds=args.max_seconds)
    bad = [r for r in results if r.get("status") in ("failed", "incomplete")]
    stopped = [r for r in results if r.get("status") == "budget_reached"]
    print(f"\ndone: {len(results) - len(bad) - len(stopped)} ok/skipped, "
          f"{len(bad)} needing a rerun")
    if stopped:
        print(f"  stopped on time budget with {stopped[0]['remaining']} runs left; "
              f"rerun the same command to continue")
    for r in bad:
        print(f"  ! {r['arm']}/seed_{r['seed']}: {r['status']}")


# ----------------------------------------------------------------- sweep -----
def cmd_sweep(args) -> None:
    base = _load(args)
    factors = args.factors or base.raw["genotypes"]["IDHmut_TP53mut"][
        "apoptosis_rate"]["factor"]["sweep"]
    physicell = _physicell_root(args)
    print(f"sweeping apoptosis factor over {factors}")

    # The time budget spans the WHOLE sweep, not each factor. Handing the same
    # max_seconds to every factor in turn would let a "400 second" chunk run for
    # 400 x len(factors) seconds, which is how a chunked run gets killed
    # mid-simulation by an outer timeout.
    started = time.time()
    stopped_early = False
    for f in factors:
        remaining = (
            None if args.max_seconds is None
            else args.max_seconds - (time.time() - started)
        )
        if remaining is not None and remaining <= 0:
            stopped_early = True
            print(f"\nstopped on time budget before factor {f}")
            break
        model = _load(args, factor=float(f))
        root = Path(args.runs) / f"sweep_factor_{f}"
        specs = runner.plan(model, root, physicell, arms=args.arms,
                            n_replicates=args.replicates, apoptosis_factor=float(f))
        print(f"\nfactor {f}: {len(specs)} runs")
        results = runner.run_all(
            specs, model, physicell, force=args.force,
            on_progress=lambda i, n, s, r: print(
                f"  [{i}/{n}] {s.arm}/seed_{s.seed} {r.get('status')}", flush=True),
            max_seconds=remaining,
        )
        if any(r.get("status") == "budget_reached" for r in results):
            stopped_early = True
            break
    if stopped_early:
        print("\nSweep incomplete. Rerun the same command to continue "
              "(completed runs are skipped).")
    else:
        print("\nSweep complete.")


# --------------------------------------------------------------- analyze -----
def cmd_analyze(args) -> None:
    model = _load(args)
    out = Path(args.out)
    fig_dir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)

    traj = collect.read_all(Path(args.runs))
    reps = outcomes.per_replicate(traj)
    arm_reps = outcomes.per_replicate_arm(traj)

    traj.to_csv(out / "trajectories.csv", index=False)
    reps.to_csv(out / "replicate_outcomes.csv", index=False)
    arm_reps.to_csv(out / "replicate_outcomes_by_arm.csv", index=False)

    summary: dict = {
        "n_runs": int(reps.groupby(["arm", "seed"]).ngroups),
        "arms": sorted(reps["arm"].unique()),
        "n_timepoints_expected": model.n_timepoints,
        "n_timepoints_observed": sorted(reps["n_timepoints_observed"].unique().tolist()),
        "primary_outcome": outcomes.PRIMARY_OUTCOME,
    }

    # --- did the simulation do what the parameters say it must? -----------
    check = outcomes.expectation_check(reps, model)
    if not check.empty:
        check.to_csv(out / "expectation_check.csv", index=False)
        plots.expectation(check, fig_dir)
        summary["expectation_check_passes"] = bool(check["passes"].all())

    # --- descriptive, per arm ---------------------------------------------
    desc = pd.concat(
        [stats.describe_arms(reps, o) for o in
         [outcomes.PRIMARY_OUTCOME, *outcomes.SECONDARY_OUTCOMES]
         if o in reps.columns],
        ignore_index=True,
    )
    desc.to_csv(out / "arm_summaries.csv", index=False)

    # --- paired contrasts, within seed ------------------------------------
    contrasts = []
    for arm in sorted(reps["arm"].unique()):
        pops = sorted(reps[reps.arm == arm]["population"].unique())
        if len(pops) == 2:
            contrasts.append(
                stats.paired_contrast(reps, arm, pops[0], pops[1],
                                      outcomes.PRIMARY_OUTCOME))
    contrasts_df = pd.DataFrame(contrasts)
    if not contrasts_df.empty:
        contrasts_df.to_csv(out / "paired_contrasts.csv", index=False)

    floor = stats.false_positive_floor(reps, outcomes.PRIMARY_OUTCOME)
    summary["false_positive_floor"] = floor
    if not contrasts_df.empty:
        plots.effect_sizes(contrasts_df, fig_dir,
                           floor if floor.get("ci_low") is not None else None)

    # --- mixed effects, determinism, convergence --------------------------
    for arm in sorted(reps["arm"].unique()):
        me = stats.mixed_effects(reps[reps.arm == arm], outcomes.PRIMARY_OUTCOME)
        summary.setdefault("mixed_effects", {})[arm] = me

    det = stats.determinism_report(reps, outcomes.PRIMARY_OUTCOME)
    det.to_csv(out / "determinism_report.csv", index=False)
    summary["max_cv_across_seeds"] = (
        float(det["cv_across_seeds"].max()) if not det.empty else None)

    conv_arm = next((a for a in reps["arm"].unique() if a.startswith("mono_")),
                    reps["arm"].iloc[0])
    conv_pop = sorted(reps[reps.arm == conv_arm]["population"].unique())[0]
    conv = stats.convergence(reps, conv_arm, conv_pop, outcomes.PRIMARY_OUTCOME)
    if not conv.empty:
        conv.to_csv(out / "convergence.csv", index=False)
        plots.convergence(conv, fig_dir)
        summary["convergence_arm"] = f"{conv_arm}/{conv_pop}"
        summary["final_relative_half_width"] = float(
            conv["relative_half_width"].iloc[-1])

    # --- competition ------------------------------------------------------
    comp = outcomes.competition_coefficients(reps)
    if not comp.empty:
        comp.to_csv(out / "competition_coefficients.csv", index=False)
        summary["competition"] = {
            f"{a}/{p}": stats.bootstrap_mean(
                g["competition_coefficient"].to_numpy()).as_dict()
            for (a, p), g in comp.groupby(["arm", "population"])
        }

    # --- figures ----------------------------------------------------------
    for arm in sorted(traj["arm"].unique()):
        plots.trajectories(traj, model, fig_dir, arm)

    # --- sweep, if sweep run directories are present ----------------------
    sweep_rows = []
    for d in sorted(Path(args.runs).glob("sweep_factor_*")):
        try:
            st = collect.read_all(d, exclude=())
        except ValueError:
            continue
        sr = outcomes.per_replicate(st)
        for pop, g in sr.groupby("population"):
            sweep_rows.append({
                "apoptosis_factor": float(d.name.split("_")[-1]),
                "population": pop,
                "mean": float(g["n_end"].mean()),
                "sd": float(g["n_end"].std(ddof=1)) if len(g) > 1 else float("nan"),
                "n_replicates": len(g),
            })
    if sweep_rows:
        sweep = pd.DataFrame(sweep_rows)
        sweep.to_csv(out / "sensitivity_sweep.csv", index=False)
        plots.sensitivity(sweep, fig_dir, model)

    # --- manuscript artefacts --------------------------------------------
    prov = (runner.physicell_provenance(_physicell_root(args))
            if (args.physicell or os.environ.get("PHYSICELL_ROOT"))
            else {"physicell_version": model.meta["physicell_version"],
                  "physicell_commit": "unknown"})
    report.write_all(model, out, prov, int(reps["seed"].nunique()))

    # --- the eight manuscript figures ------------------------------------
    try:
        made = manuscript_figures.build_all(
            traj, reps, arm_reps, contrasts_df, floor, model,
            Path(args.runs), out / "manuscript_figures",
            arm=getattr(args, "arm", "co_separated"))
        summary["manuscript_figures"] = [p.name for p in made]
    except Exception as exc:  # noqa: BLE001
        summary["manuscript_figures_error"] = str(exc)
        print(f"  manuscript figures FAILED: {exc}")

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote analysis to {out}")
    print(f"  runs analysed        {summary['n_runs']}")
    print(f"  timepoints           expected {model.n_timepoints}, "
          f"observed {summary['n_timepoints_observed']}")
    if "expectation_check_passes" in summary:
        print(f"  matches closed form  {summary['expectation_check_passes']}")
    if summary.get("max_cv_across_seeds") is not None:
        print(f"  max CV across seeds  {summary['max_cv_across_seeds']:.3%}"
              "   (non-zero => a single run was never sufficient)")
    print(f"  figures              {fig_dir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="glioma", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--params", default=str(DEFAULT_PARAMS))
    p.add_argument("--runs", default=str(DEFAULT_RUNS))
    p.add_argument("--physicell", default=None,
                   help="path to the PhysiCell source tree (or set PHYSICELL_ROOT)")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("check", help="validate the model, run nothing")
    c.set_defaults(func=cmd_check)

    t = sub.add_parser("table", help="write Table 1, Methods and reviewer map")
    t.add_argument("--out", default=str(DEFAULT_RESULTS))
    t.add_argument("--replicates", type=int, default=None)
    t.set_defaults(func=cmd_table)

    r = sub.add_parser("run", help="run seeded replicates")
    r.add_argument("--arms", nargs="*", default=None)
    r.add_argument("--replicates", type=int, default=None)
    r.add_argument("--apoptosis-factor", type=float, default=None)
    r.add_argument("--force", action="store_true", help="rerun completed runs")
    r.add_argument("--max-seconds", type=float, default=None,
                   help="wall-clock budget; stops between runs so a long set can "
                        "be worked through in chunks. Rerun to continue.")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("sweep", help="run the apoptosis-factor sensitivity sweep")
    s.add_argument("--factors", nargs="*", type=float, default=None)
    s.add_argument("--arms", nargs="*", default=["mono_TP53wt", "mono_TP53mut"])
    s.add_argument("--replicates", type=int, default=None)
    s.add_argument("--force", action="store_true")
    s.add_argument("--max-seconds", type=float, default=None,
                   help="wall-clock budget per sweep factor; stops between runs.")
    s.set_defaults(func=cmd_sweep)

    a = sub.add_parser("analyze", help="outcomes, statistics, figures, report")
    a.add_argument("--out", default=str(DEFAULT_RESULTS))
    a.add_argument("--arm", default="co_separated",
                   help="arm used for the manuscript figures (default co_separated)")
    a.set_defaults(func=cmd_analyze)
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
