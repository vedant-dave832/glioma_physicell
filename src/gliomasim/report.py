"""Emit the manuscript artefacts: Table 1, the Methods text, and a reviewer map.

These are generated from the same registry object the simulation is built from,
so a parameter cannot appear one way in the XML and another way in the table.
That coupling is the whole point -- reviewer major issue #2 is a list of places
where the manuscript and the model disagreed.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd

from . import convert
from .registry import Model

PROVENANCE_ORDER = {"measured": 0, "calculated": 1, "default": 2, "assumed": 3}


def table1(model: Model) -> pd.DataFrame:
    """The corrected parameter table the reviewer asked for, column by column.

    "a corrected parameter table listing each parameter, its unit, exact
    PhysiCell field, biological meaning, source, and whether it was measured,
    calculated, assumed, or retained as a default."
    """
    rows = []
    for p in model.params:
        src = model.sources.get(p.source or "", {})
        citation = src.get("citation", "").strip() if p.source else ""
        if src.get("type") == "none":
            citation = "NO SOURCE — assumption, swept"
        rows.append(
            {
                "Parameter": p.name,
                "Value": p.value,
                "Unit": p.units,
                "PhysiCell field": p.field_path,
                "Biological meaning": " ".join(str(p.meaning).split()),
                "Provenance": p.provenance,
                "Source": citation or ("PhysiCell default" if p.provenance == "default" else "—"),
                "Uncertainty": (
                    f"{p.uncertainty['low']}–{p.uncertainty['high']} ({p.uncertainty['kind']})"
                    if p.uncertainty
                    else ""
                ),
                "Swept over": ", ".join(str(s) for s in p.sweep) if p.sweep else "",
            }
        )
    df = pd.DataFrame(rows)
    df["_o"] = df["Provenance"].map(PROVENANCE_ORDER).fillna(9)
    return df.sort_values(["_o", "Parameter"]).drop(columns="_o").reset_index(drop=True)


def provenance_summary(model: Model) -> pd.DataFrame:
    t = table1(model)
    return (
        t.groupby("Provenance").size().rename("n_parameters").reset_index()
        .sort_values("n_parameters", ascending=False)
    )


def methods(model: Model, provenance: dict, n_replicates: int) -> str:
    """The Methods paragraph, with every number computed rather than typed."""
    n_tp = model.n_timepoints
    wt = model.genotypes["IDHmut_TP53wt"]
    mt = model.genotypes["IDHmut_TP53mut"]
    factor = wt.apoptosis_rate_per_min / mt.apoptosis_rate_per_min
    naive = convert.naive_phase_duration_error(convert.days_to_min(wt.doubling_time_days))
    nx = int((model.domain["x_max"] - model.domain["x_min"]) / model.domain["dx"])
    ny = int((model.domain["y_max"] - model.domain["y_min"]) / model.domain["dx"])

    return f"""\
## Methods (generated {date.today().isoformat()})

### Simulation platform

Agent-based simulations were performed with PhysiCell version \
{provenance.get('physicell_version', '?')} (git commit \
{provenance.get('physicell_commit', '?')[:12]}), compiled from source and run \
from the command line on a local installation rather than through the NanoHub \
web interface. The unmodified `template` sample project was used; no custom C++ \
was written. All model specification is contained in the XML configuration and \
the initial-condition CSV emitted by the accompanying harness, both of which are \
archived with the run outputs. Runs executed with \
{model.run['omp_num_threads']} OpenMP threads; because PhysiCell seeds one \
random-number stream per thread from the base seed, the thread count is part of \
the reproducibility specification and was held constant across all runs.

### Domain and discretisation

The domain spans {model.domain['x_min']:.0f} to {model.domain['x_max']:.0f} µm \
in x and y, discretised at {model.domain['dx']:.0f} µm, giving a {nx} × {ny} \
mesh of {model.n_voxels} voxels, restricted to the z = 0 plane. Total simulated \
time was {model.time['max_time']:.0f} min ({model.max_time_days:.2f} days) with \
full output saved every {model.time['save_interval']:.0f} min. This yields \
{int(model.time['max_time'] / model.time['save_interval'])} intervals and \
therefore **{n_tp} saved timepoints** including t = 0 and the endpoint. Diffusion, \
mechanics and phenotype timesteps were {model.time['dt_diffusion']}, \
{model.time['dt_mechanics']} and {model.time['dt_phenotype']} min respectively.

### Cell definitions and the doubling-time conversion

All agents carry the IDH1 R132H mutation. The comparison is TP53-wild-type \
versus TP53-mutant on that shared background; ATRX and CDKN2A/B status are held \
constant and are not represented, and CNS WHO grade 4 is assumed for the model \
rather than demonstrated by it.

Proliferation uses the PhysiCell Live cycle model (code 5), a single phase with \
an exponentially distributed duration and transition rate *r* (min⁻¹). Because \
each transition yields one additional cell and apoptosis removes cells at hazard \
*d*, the population satisfies dN/dt = (r − d)N, so the population doubling time \
is ln(2)/(r − d) — **not** the mean phase duration. Published doubling times \
were converted to cycle transition rates by inverting that relation, \
r = ln(2)/T_double + d, since a doubling time measured in culture already nets \
out cell loss. Setting the phase duration equal to the published doubling time, \
as is sometimes done, would encode a doubling of \
{naive['naive_actual_doubling_days']:.2f} days rather than the intended \
{naive['intended_doubling_days']:.2f} days, an error of a factor of ln(2).

The TP53-wild-type agent was given a population doubling time of \
{wt.doubling_time_days:.1f} days, the mean reported for {'' }\
patient-derived IDH1-mutant glioma cultures, yielding r = \
{wt.cycle.cycle_rate_per_min:.4e} min⁻¹ (mean phase duration \
{wt.cycle.mean_phase_duration_min:.0f} min) at the default apoptosis hazard of \
{wt.apoptosis_rate_per_min:.4e} min⁻¹. The TP53-mutant agent was given the \
**same cycle transition rate** and an apoptosis hazard reduced \
{factor:.0f}-fold to {mt.apoptosis_rate_per_min:.4e} min⁻¹. Its population \
doubling time of {mt.doubling_time_days:.2f} days is therefore an output of \
those two settings, not an input.

The magnitude of the apoptosis reduction is an assumption. No study was \
identified that measures a basal per-minute apoptosis hazard for a named TP53 \
substitution in IDH1-mutant glioma, so this parameter is reported as unsourced \
and is varied across a sensitivity sweep rather than asserted; results are \
presented as a function of it. No variant-specific mechanism is claimed.

Motility was disabled and mechanical parameters were identical across all cell \
definitions, so the model contains no genotype-specific spatial phenotype and \
spatial differences between populations reflect only their seeding positions.

### Replication and analysis

Each condition was run as {n_replicates} independent replicates differing only \
in random seed, with seeds recorded in each run's manifest. Replicate *i* of \
every arm shares a seed, so arms are compared within seed. **The unit of \
analysis is the run, not the timepoint**: each run was reduced to \
replicate-level outcomes — the fitted log growth rate (primary), endpoint count, \
log₂ fold change, area under the log-count curve, final proportion, Shannon \
diversity change and final radius of gyration — before any comparison was made. \
Successive timepoints within a run are not treated as independent observations.

Effects are reported as within-seed paired mean differences with percentile \
bootstrap 95% confidence intervals, corroborated by a mixed-effects model with \
seed as a random intercept. A replicate-count convergence check reports the \
relative CI half-width as a function of *n*. Monoculture arms were run for every \
genotype and checked against the closed-form expectation N(t) = N₀·exp((r−d)t). \
Position-swapped and well-mixed layouts, and a label control in which two \
identically parameterised populations are seeded and compared, bound the \
false-positive floor of the analysis; competition coefficients are defined as \
the difference between a genotype's growth rate in co-culture and in \
monoculture at the same seed. Runs that terminated before producing all \
{n_tp} timepoints were rerun, not trimmed.

### Sources cited by parameters

{_source_list(model)}
"""


def _source_list(model: Model) -> str:
    """Only the sources actually referenced by a parameter, with what each supports."""
    used = {p.source for p in model.params if p.source}
    lines = []
    for key in sorted(used):
        src = model.sources.get(key, {})
        citation = " ".join(str(src.get("citation", "")).split())
        doi = src.get("doi")
        lines.append(f"- **{key}** — {citation}" + (f" doi:{doi}" if doi else ""))
        for s in src.get("supports", []):
            lines.append(f"    - supports: {' '.join(str(s).split())}")
        for c in src.get("caveats", []):
            lines.append(f"    - caveat: {' '.join(str(c).split())}")
    return "\n".join(lines)


REVIEWER_MAP = {
    "Major 1 — tumour-type framing": (
        "params/model.yaml: all genotypes carry IDH1 R132H; comparison is "
        "TP53-wt vs TP53-mut on that background. ATRX/CDKN2A/B held constant "
        "and declared. IDH-wildtype arm marked include_in_primary: false and "
        "labelled a reference, not a subclone. Grade 4 declared as assumed."
    ),
    "Major 2 — parameter table and provenance": (
        "report.table1() emits the unit / PhysiCell field / meaning / source / "
        "provenance table from the same registry the XML is built from. "
        "registry._validate() refuses to load if a 'measured' parameter points "
        "at a non-primary source, or an unsourced assumption is not swept."
    ),
    "Major 2b — doubling-time conversion": (
        "convert.doubling_time_to_cycle_rate: r = ln(2)/T + d, with the "
        "ln(2) error quantified by convert.naive_phase_duration_error."
    ),
    "Major 3 — counts must reconcile": (
        "convert.n_saved_timepoints computes the timepoint count (23, not 33) "
        "and raises if max_time is not a whole multiple of the save interval. "
        "outcomes.observed_doubling_time_days is fitted from the data and kept "
        "in a separate column from the entered parameters."
    ),
    "Major 4 — replication and statistics": (
        "run.plan generates seeded replicates; outcomes.per_replicate reduces "
        "each run to one row; stats works only on replicate-level rows; "
        "stats.convergence justifies n; stats.determinism_report tests the "
        "determinism claim; incomplete runs are excluded and reported by "
        "collect.read_all rather than trimmed."
    ),
    "Major 5 — competition and spatial claims": (
        "Monoculture, well-mixed, position-swapped and label-control arms in "
        "params/model.yaml; outcomes.competition_coefficients measures growth "
        "against monoculture; motility and mechanics reported in Table 1 and "
        "identical across genotypes; minimum oxygen reached is recorded so the "
        "non-limiting claim is an observation."
    ),
    "Major 8 — missing figures": "plots.py emits every figure with a companion CSV.",
    "Recommended 2 — sensitivity analysis": (
        "The apoptosis suppression factor is a swept axis "
        "(sweep: [1, 2, 5, 10, 20, 100]); plots.sensitivity renders the result "
        "as a function of the assumption."
    ),
}


def write_all(
    model: Model, out_dir: Path, provenance: dict, n_replicates: int
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    t = table1(model)
    paths = {
        "table1_csv": out_dir / "table1_parameters.csv",
        "table1_md": out_dir / "table1_parameters.md",
        "methods_md": out_dir / "methods.md",
        "reviewer_map": out_dir / "reviewer_map.json",
        "provenance_csv": out_dir / "provenance_summary.csv",
    }
    t.to_csv(paths["table1_csv"], index=False)
    paths["table1_md"].write_text(t.to_markdown(index=False))
    paths["methods_md"].write_text(methods(model, provenance, n_replicates))
    paths["reviewer_map"].write_text(json.dumps(REVIEWER_MAP, indent=2))
    provenance_summary(model).to_csv(paths["provenance_csv"], index=False)
    return paths
