"""Emit a PhysiCell settings XML and a cells.csv for one (arm, seed) run.

Rather than templating XML by hand, this mutates a copy of the stock
``sample_projects/template/config/PhysiCell_settings.xml`` shipped with the
pinned PhysiCell version. That way the emitted config is guaranteed to contain
every element the binary expects, including any this harness never touches, and
it stays valid if the schema gains elements in a later version.

No C++ is modified. The stock ``template`` project is used unchanged, with
``number_of_cells`` set to 0 so that every agent comes from cells.csv. This is
what makes the runs reproducible from this repository alone.
"""

from __future__ import annotations

import copy
import math
import random
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from .registry import Arm, Model


@dataclass
class RunSpec:
    """Everything that distinguishes one PhysiCell invocation from another."""

    arm: str
    seed: int
    apoptosis_factor: float
    out_dir: Path
    config_path: Path
    cells_csv_path: Path
    omp_num_threads: int
    populations: list[dict]  # resolved: label, genotype, centre, n, type_id


def _set(node: ET.Element, path: str, value) -> None:
    target = node.find(path)
    if target is None:
        raise KeyError(f"XML path not found in template: {path}")
    target.text = str(value)


def _require(node: ET.Element, path: str) -> ET.Element:
    target = node.find(path)
    if target is None:
        raise KeyError(f"XML path not found in template: {path}")
    return target


def _build_cell_definition(
    base: ET.Element,
    *,
    name: str,
    type_id: int,
    cycle_rate: float,
    apoptosis_rate: float,
    shared: dict,
    substrate_name: str,
    all_names: list[str],
) -> ET.Element:
    cd = copy.deepcopy(base)
    cd.set("name", name)
    cd.set("ID", str(type_id))

    # --- cycle: swap the template's flow-cytometry model for Live (code 5) ---
    phenotype = _require(cd, "phenotype")
    old_cycle = _require(phenotype, "cycle")
    idx = list(phenotype).index(old_cycle)
    phenotype.remove(old_cycle)

    cycle = ET.Element("cycle", {"code": "5", "name": "Live"})
    rates = ET.SubElement(cycle, "phase_transition_rates", {"units": "1/min"})
    rate = ET.SubElement(
        rates,
        "rate",
        {"start_index": "0", "end_index": "0", "fixed_duration": "false"},
    )
    rate.text = f"{cycle_rate:.10g}"
    phenotype.insert(idx, cycle)

    # --- death ------------------------------------------------------------
    _set(cd, "phenotype/death/model[@code='100']/death_rate", f"{apoptosis_rate:.10g}")
    _set(cd, "phenotype/death/model[@code='101']/death_rate", f"{shared['necrosis_rate']:.10g}")

    # --- mechanics and motility, identical across genotypes ---------------
    _set(cd, "phenotype/mechanics/cell_cell_adhesion_strength", shared["cell_cell_adhesion"])
    _set(cd, "phenotype/mechanics/cell_cell_repulsion_strength", shared["cell_cell_repulsion"])
    _set(cd, "phenotype/motility/options/enabled", str(bool(shared["motility_enabled"])).lower())

    # --- rename the substrate everywhere it appears -----------------------
    sec = _require(cd, "phenotype/secretion/substrate")
    sec.set("name", substrate_name)
    chemo = cd.find("phenotype/motility/options/chemotaxis/substrate")
    if chemo is not None:
        chemo.text = substrate_name
    adv = cd.find(
        "phenotype/motility/options/advanced_chemotaxis/chemotactic_sensitivities/"
        "chemotactic_sensitivity"
    )
    if adv is not None:
        adv.set("substrate", substrate_name)

    # --- per-cell-type lists must name every cell definition --------------
    def _fill(list_path: str, item_tag: str, value: str = "0") -> None:
        container = cd.find(list_path)
        if container is None:
            return
        for child in list(container):
            container.remove(child)
        for other in all_names:
            el = ET.SubElement(container, item_tag, {"name": other})
            el.text = value

    _fill("phenotype/mechanics/cell_adhesion_affinities", "cell_adhesion_affinity", "1")
    _fill("phenotype/cell_interactions/live_phagocytosis_rates", "phagocytosis_rate")
    _fill("phenotype/cell_interactions/attack_rates", "attack_rate")
    _fill("phenotype/cell_interactions/fusion_rates", "fusion_rate")
    _fill("phenotype/cell_transformations/transformation_rates", "transformation_rate")

    for el in cd.findall("phenotype/cell_interactions/live_phagocytosis_rates/phagocytosis_rate"):
        el.set("units", "1/min")
    for el in cd.findall("phenotype/cell_interactions/attack_rates/attack_rate"):
        el.set("units", "1/min")
    for el in cd.findall("phenotype/cell_transformations/transformation_rates/transformation_rate"):
        el.set("units", "1/min")

    return cd


def resolve_populations(model: Model, arm: Arm) -> list[dict]:
    """Assign a distinct cell-definition name and type id to each population.

    The label control needs two definitions with identical parameters and
    different names, so populations -- not genotypes -- own the naming.
    """
    pops = []
    for i, pop in enumerate(arm.populations):
        label = pop.label_as or pop.genotype
        g = model.genotypes[pop.genotype]
        pops.append(
            {
                "label": label,
                "genotype": pop.genotype,
                "centre": pop.centre,
                "n": model.seeding["n_cells_per_population"],
                "type_id": i,
                "cycle_rate": g.cycle.cycle_rate_per_min,
                "apoptosis_rate": g.apoptosis_rate_per_min,
            }
        )
    return pops


def write_cells_csv(
    path: Path, pops: list[dict], radius: float, layout: str | None, seed: int
) -> None:
    """Seed each population as a uniform disk. v1 CSV: x,y,z,typeID (no header).

    The RNG for cell placement is seeded from the run seed so that the initial
    condition varies between replicates exactly as the simulation does. A
    replicate that reuses the same starting arrangement is not an independent
    replicate.
    """
    rng = random.Random(seed)
    lines: list[str] = []
    for pop in pops:
        cx, cy = pop["centre"]
        for _ in range(pop["n"]):
            # uniform over the disk: sqrt keeps density flat with radius
            r = radius * math.sqrt(rng.random())
            theta = 2.0 * math.pi * rng.random()
            lines.append(
                f"{cx + r * math.cos(theta):.4f},{cy + r * math.sin(theta):.4f},"
                f"0.0,{pop['type_id']}"
            )
    if layout == "mixed":
        rng.shuffle(lines)
    path.write_text("\n".join(lines) + "\n")


def build(
    model: Model,
    arm_key: str,
    seed: int,
    run_root: Path,
    physicell_root: Path,
    apoptosis_factor: float,
) -> RunSpec:
    arm = model.arms[arm_key]
    pops = resolve_populations(model, arm)

    run_dir = run_root / arm_key / f"seed_{seed}"
    config_dir = run_dir / "config"
    out_dir = run_dir / "output"
    config_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    template_path = (
        physicell_root / "sample_projects" / "template" / "config" / "PhysiCell_settings.xml"
    )
    tree = ET.parse(template_path)
    root = tree.getroot()

    d = model.domain
    for key in ("x_min", "x_max", "y_min", "y_max"):
        _set(root, f"domain/{key}", d[key])
    for key in ("dx", "dy", "dz"):
        _set(root, f"domain/{key}", d["dx"])
    _set(root, "domain/use_2D", str(bool(d["use_2D"])).lower())

    t = model.time
    _set(root, "overall/max_time", t["max_time"])
    _set(root, "overall/dt_diffusion", t["dt_diffusion"])
    _set(root, "overall/dt_mechanics", t["dt_mechanics"])
    _set(root, "overall/dt_phenotype", t["dt_phenotype"])

    _set(root, "parallel/omp_num_threads", model.run["omp_num_threads"])

    _set(root, "save/folder", str(out_dir.resolve()))
    _set(root, "save/full_data/interval", t["save_interval"])
    _set(root, "save/SVG/interval", t["save_interval"])
    _set(root, "options/random_seed", seed)

    # --- substrate --------------------------------------------------------
    s = model.substrate
    var = _require(root, "microenvironment_setup/variable")
    var.set("name", s["name"])
    var.set("units", "mmHg")
    _set(var, "physical_parameter_set/diffusion_coefficient", s["diffusion"])
    _set(var, "physical_parameter_set/decay_rate", s["decay"])
    _set(var, "initial_condition", s["initial"])
    dirichlet = _require(var, "Dirichlet_boundary_condition")
    dirichlet.text = str(s["dirichlet"])
    dirichlet.set("enabled", "True")
    for bid in ("xmin", "xmax", "ymin", "ymax"):
        bv = _require(var, f"Dirichlet_options/boundary_value[@ID='{bid}']")
        bv.text = str(s["dirichlet"])
        bv.set("enabled", "True")

    # --- cell definitions -------------------------------------------------
    cds = _require(root, "cell_definitions")
    base = _require(cds, "cell_definition")
    for child in list(cds):
        cds.remove(child)
    names = [p["label"] for p in pops]
    for pop in pops:
        cds.append(
            _build_cell_definition(
                base,
                name=pop["label"],
                type_id=pop["type_id"],
                cycle_rate=pop["cycle_rate"],
                apoptosis_rate=pop["apoptosis_rate"],
                shared=model.shared,
                substrate_name=s["name"],
                all_names=names,
            )
        )
        # every cell type takes up oxygen at the shared default rate
        cd = cds[-1]
        _set(cd, f"phenotype/secretion/substrate[@name='{s['name']}']/uptake_rate", s["uptake"])

    # --- initial conditions from CSV, none placed at random ---------------
    cells_csv = config_dir / "cells.csv"
    cp = _require(root, "initial_conditions/cell_positions")
    cp.set("enabled", "true")
    cp.set("type", "csv")
    _set(cp, "folder", str(config_dir.resolve()))
    _set(cp, "filename", "cells.csv")
    _set(root, "user_parameters/number_of_cells", 0)

    write_cells_csv(
        cells_csv,
        pops,
        radius=model.seeding["disk_radius"]["value"],
        layout=arm.layout,
        seed=seed,
    )

    config_path = config_dir / "PhysiCell_settings.xml"
    ET.indent(tree, space="    ")
    tree.write(config_path, encoding="utf-8", xml_declaration=False)

    return RunSpec(
        arm=arm_key,
        seed=seed,
        apoptosis_factor=apoptosis_factor,
        out_dir=out_dir,
        config_path=config_path,
        cells_csv_path=cells_csv,
        omp_num_threads=model.run["omp_num_threads"],
        populations=pops,
    )
