"""Loads params/model.yaml + params/sources.yaml and resolves them into settings.

The point of this layer: the PhysiCell XML and the manuscript's Table 1 are
generated from the same object. They cannot disagree, which is the class of
error reviewer major issue #2 is about.

It also enforces the provenance discipline. A parameter marked ``measured`` must
name a source of type ``primary``. A parameter marked ``assumed`` must either be
swept or be explicitly acknowledged. Violations raise at load time rather than
turning into a sentence in a manuscript.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from . import convert

VALID_PROVENANCE = {"measured", "calculated", "assumed", "default"}


class RegistryError(RuntimeError):
    pass


@dataclass
class Param:
    """One parameter, with everything Table 1 needs to display it."""

    name: str
    value: Any
    units: str
    field_path: str
    meaning: str
    provenance: str
    source: str | None
    uncertainty: dict | None = None
    sweep: list | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        if self.provenance not in VALID_PROVENANCE:
            raise RegistryError(
                f"{self.name}: provenance {self.provenance!r} not in {VALID_PROVENANCE}"
            )


@dataclass
class Genotype:
    key: str
    display_name: str
    colour: str
    include_in_primary: bool
    doubling_time_days: float
    apoptosis_rate_per_min: float
    cycle: convert.CycleSetting
    params: list[Param] = field(default_factory=list)
    note: str | None = None


@dataclass
class Population:
    genotype: str
    centre: tuple[float, float]
    label_as: str | None = None


@dataclass
class Arm:
    key: str
    kind: str
    layout: str | None
    populations: list[Population]
    purpose: str
    identical_parameters: bool = False


@dataclass
class Model:
    meta: dict
    domain: dict
    time: dict
    substrate: dict
    shared: dict
    genotypes: dict[str, Genotype]
    arms: dict[str, Arm]
    seeding: dict
    run: dict
    sources: dict
    raw: dict
    params: list[Param]

    # -- computed, never asserted -------------------------------------------
    @property
    def n_timepoints(self) -> int:
        return convert.n_saved_timepoints(
            self.time["max_time"], self.time["save_interval"]
        )

    @property
    def max_time_days(self) -> float:
        return convert.min_to_days(self.time["max_time"])

    @property
    def n_voxels(self) -> int:
        nx = (self.domain["x_max"] - self.domain["x_min"]) / self.domain["dx"]
        ny = (self.domain["y_max"] - self.domain["y_min"]) / self.domain["dx"]
        return int(round(nx * ny))

    def primary_genotypes(self) -> list[Genotype]:
        return [g for g in self.genotypes.values() if g.include_in_primary]


def _leaf(name: str, node: dict) -> Param:
    return Param(
        name=name,
        value=node.get("value"),
        units=node.get("units", ""),
        field_path=node.get("field", ""),
        meaning=node.get("meaning", ""),
        provenance=node.get("provenance", "assumed"),
        source=node.get("source"),
        uncertainty=node.get("uncertainty"),
        sweep=node.get("sweep"),
        notes=node.get("sweep_note") or node.get("note"),
    )


def _is_leaf(node: Any) -> bool:
    return isinstance(node, dict) and "provenance" in node


def _collect(prefix: str, node: Any, out: list[Param]) -> None:
    if _is_leaf(node):
        out.append(_leaf(prefix, node))
        return
    if isinstance(node, dict):
        for k, v in node.items():
            _collect(f"{prefix}.{k}" if prefix else k, v, out)


def load(params_dir: str | Path, apoptosis_factor: float | None = None) -> Model:
    """Load and resolve the model.

    Parameters
    ----------
    apoptosis_factor
        Override for the TP53 apoptosis suppression factor. This is the swept
        axis; passing a value here is how the sensitivity analysis walks it.
    """
    params_dir = Path(params_dir)
    raw = yaml.safe_load((params_dir / "model.yaml").read_text())
    sources = yaml.safe_load((params_dir / "sources.yaml").read_text())["sources"]

    all_params: list[Param] = []
    _collect("", raw.get("domain", {}), all_params)
    _collect("", raw.get("time", {}), all_params)
    _collect("", raw.get("substrate", {}), all_params)
    _collect("", raw.get("shared_phenotype", {}), all_params)

    scalars = lambda section: {  # noqa: E731
        k: v["value"] for k, v in raw[section].items() if _is_leaf(v)
    }
    domain = scalars("domain")
    time = scalars("time")
    substrate = scalars("substrate")
    substrate["name"] = raw["substrate"]["name"]
    shared = scalars("shared_phenotype")

    # --- genotypes, with derived values resolved --------------------------
    genotypes: dict[str, Genotype] = {}
    resolved_apoptosis: dict[str, float] = {}

    # first pass: direct apoptosis rates
    for key, g in raw["genotypes"].items():
        node = g["apoptosis_rate"]
        if "derived_from" not in node:
            resolved_apoptosis[key] = float(node["value"])

    # second pass: derived apoptosis rates
    for key, g in raw["genotypes"].items():
        node = g["apoptosis_rate"]
        if "derived_from" in node:
            parent = node["derived_from"]
            if parent not in resolved_apoptosis:
                raise RegistryError(
                    f"{key}: derived_from {parent!r}, which is itself derived or missing"
                )
            factor = float(
                apoptosis_factor
                if apoptosis_factor is not None
                else node["factor"]["value"]
            )
            if factor <= 0:
                raise RegistryError(f"{key}: apoptosis factor must be positive")
            op = node.get("operation", "divide")
            base = resolved_apoptosis[parent]
            resolved_apoptosis[key] = base / factor if op == "divide" else base * factor

    # Genotypes that state a doubling time are resolved first; genotypes that
    # match another genotype's cycle rate are resolved second, from that result.
    order = sorted(raw["genotypes"], key=lambda k: "cycle" in raw["genotypes"][k])
    cycle_rates: dict[str, float] = {}

    for key in order:
        g = raw["genotypes"][key]
        d = resolved_apoptosis[key]
        gparams: list[Param] = []

        if "cycle" in g:
            # Cycle rate is inherited. Doubling time becomes an OUTPUT.
            parent = g["cycle"]["match_to"]
            if parent not in cycle_rates:
                raise RegistryError(
                    f"{key}: cycle matched to {parent!r}, which is not resolved first"
                )
            r = cycle_rates[parent]
            cycle = convert.CycleSetting(
                cycle_rate_per_min=r,
                mean_phase_duration_min=1.0 / r,
                apoptosis_rate_per_min=d,
                target_doubling_time_min=convert.doubling_time_from_rates(r, d),
                net_growth_rate_per_min=r - d,
            )
            dt_days = convert.min_to_days(cycle.implied_doubling_time_min)
            dt_node = {"source": g["cycle"].get("source")}
            gparams.append(
                Param(
                    name=f"genotype.{key}.doubling_time",
                    value=round(dt_days, 4),
                    units="days",
                    field_path="OUTPUT (not entered into the model)",
                    meaning=(
                        f"population doubling time implied by the inherited cycle "
                        f"rate and this genotype's apoptosis rate; matched to "
                        f"{parent} on cycle rate, so any growth difference is a "
                        f"direct arithmetic consequence of the apoptosis setting"
                    ),
                    provenance="calculated",
                    source=g["cycle"].get("source"),
                )
            )
            gparams.append(_leaf(f"genotype.{key}.cycle_match", g["cycle"] | {"value": parent, "units": "genotype key", "field": "cell_definition/phenotype/cycle[code=5]/phase_transition_rates/rate[0->0]"}))
        else:
            dt_node = g["doubling_time"]
            dt_days = float(dt_node["value"])
            cycle = convert.doubling_time_to_cycle_rate(convert.days_to_min(dt_days), d)
            _collect(f"genotype.{key}", {"doubling_time": dt_node}, gparams)

        cycle_rates[key] = cycle.cycle_rate_per_min
        ap_node = g["apoptosis_rate"]
        if "derived_from" in ap_node:
            p = _leaf(f"genotype.{key}.apoptosis_factor", ap_node["factor"])
            gparams.append(p)
            gparams.append(
                Param(
                    name=f"genotype.{key}.apoptosis_rate",
                    value=d,
                    units="1/min",
                    field_path=ap_node["factor"]["field"],
                    meaning=(
                        f"basal apoptosis hazard = {ap_node['derived_from']} rate "
                        f"/ {ap_node['factor']['value']}"
                    ),
                    provenance="calculated",
                    source=ap_node["factor"]["source"],
                )
            )
        else:
            gparams.append(_leaf(f"genotype.{key}.apoptosis_rate", ap_node))

        gparams.append(
            Param(
                name=f"genotype.{key}.cycle_transition_rate",
                value=cycle.cycle_rate_per_min,
                units="1/min",
                field_path=(
                    "cell_definition/phenotype/cycle[code=5]/"
                    "phase_transition_rates/rate[0->0]"
                ),
                meaning=(
                    (
                        f"Live-cycle 0->0 transition rate inherited from "
                        f"{g['cycle']['match_to']}"
                        if "cycle" in g
                        else f"Live-cycle 0->0 transition rate solved from a "
                        f"{dt_days}-day population doubling time at an apoptosis "
                        f"rate of {d:.6g} /min, via r = ln(2)/T + d"
                    )
                    + f" (mean phase duration {cycle.mean_phase_duration_min:.1f} min)"
                ),
                provenance="calculated",
                source=dt_node.get("source"),
            )
        )
        all_params.extend(gparams)

        genotypes[key] = Genotype(
            key=key,
            display_name=g["display_name"],
            colour=g.get("colour", "black"),
            include_in_primary=g.get("include_in_primary", True),
            doubling_time_days=dt_days,
            apoptosis_rate_per_min=d,
            cycle=cycle,
            params=gparams,
            note=g.get("note"),
        )

    arms = {
        key: Arm(
            key=key,
            kind=a["kind"],
            layout=a.get("layout"),
            populations=[
                Population(
                    genotype=p["genotype"],
                    centre=tuple(p["centre"]),
                    label_as=p.get("label_as"),
                )
                for p in a["populations"]
            ],
            purpose=a["purpose"],
            identical_parameters=a.get("identical_parameters", False),
        )
        for key, a in raw["arms"].items()
    }

    model = Model(
        meta=raw["meta"],
        domain=domain,
        time=time,
        substrate=substrate,
        shared=shared,
        genotypes=genotypes,
        arms=arms,
        seeding=raw["seeding"],
        run=raw["run"],
        sources=sources,
        raw=raw,
        params=all_params,
    )
    _validate(model)
    return model


def _validate(model: Model) -> None:
    problems: list[str] = []

    for p in model.params:
        if p.provenance == "measured":
            src = model.sources.get(p.source or "")
            if src is None:
                problems.append(f"{p.name}: provenance 'measured' but source {p.source!r} is not in sources.yaml")
            elif src.get("type") != "primary":
                problems.append(
                    f"{p.name}: provenance 'measured' but source {p.source!r} has "
                    f"type {src.get('type')!r}; a measured value needs a study that measured it"
                )
        if p.provenance == "assumed" and p.source:
            src = model.sources.get(p.source, {})
            if src.get("type") == "none" and not p.sweep:
                problems.append(
                    f"{p.name}: assumed, unsourced, and not swept. Either give it a "
                    f"sweep range or find a source."
                )

    # arms must reference real genotypes
    for arm in model.arms.values():
        for pop in arm.populations:
            if pop.genotype not in model.genotypes:
                problems.append(f"arm {arm.key}: unknown genotype {pop.genotype!r}")

    # timepoint arithmetic must be exact
    try:
        model.n_timepoints
    except ValueError as exc:
        problems.append(f"time: {exc}")

    # every genotype's derived setting must round-trip
    for g in model.genotypes.values():
        implied = convert.min_to_days(g.cycle.implied_doubling_time_min)
        if not math.isclose(implied, g.doubling_time_days, rel_tol=1e-9):
            problems.append(
                f"genotype {g.key}: derived cycle rate implies a {implied:.4f}-day "
                f"doubling but the target was {g.doubling_time_days}"
            )

    if problems:
        raise RegistryError(
            "model.yaml failed validation:\n  - " + "\n  - ".join(problems)
        )
