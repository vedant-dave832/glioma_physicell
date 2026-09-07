"""Execute PhysiCell replicates and record a manifest for each one.

Every run writes a manifest.json next to its output. The manifest is what makes
the Methods section checkable: it records the PhysiCell version and git commit,
the seed, the thread count, the resolved parameter values, the wall time, and
the hash of the config actually handed to the binary. A result whose manifest is
missing is not a result.

Runs are resumable. A run whose output already contains the expected number of
saved timepoints is skipped, so a long replicate set can be executed in chunks.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from dataclasses import asdict
from pathlib import Path

from .build_config import RunSpec, build
from .registry import Model


class RunError(RuntimeError):
    pass


def physicell_provenance(physicell_root: Path) -> dict:
    version_file = physicell_root / "VERSION.txt"
    version = version_file.read_text().strip() if version_file.exists() else "unknown"
    commit = "unknown"
    try:
        commit = (
            subprocess.check_output(
                ["git", "-C", str(physicell_root), "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        pass
    return {"physicell_version": version, "physicell_commit": commit}


def _expected_outputs(model: Model) -> int:
    return model.n_timepoints


def is_complete(spec: RunSpec, model: Model) -> bool:
    n = len(list(spec.out_dir.glob("output*_cells.mat")))
    return n >= _expected_outputs(model)


def execute(
    spec: RunSpec,
    model: Model,
    physicell_root: Path,
    binary: Path | None = None,
    timeout_s: float | None = None,
) -> dict:
    binary = binary or (physicell_root / "project")
    if not binary.exists():
        raise RunError(
            f"PhysiCell binary not found at {binary}. Build it first:\n"
            f"  cd {physicell_root} && make template && make -j"
        )

    config_bytes = spec.config_path.read_bytes()
    manifest = {
        "arm": spec.arm,
        "seed": spec.seed,
        "apoptosis_factor": spec.apoptosis_factor,
        "omp_num_threads": spec.omp_num_threads,
        "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "cells_csv_sha256": hashlib.sha256(spec.cells_csv_path.read_bytes()).hexdigest(),
        "populations": spec.populations,
        "expected_timepoints": _expected_outputs(model),
        **physicell_provenance(physicell_root),
    }

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(spec.omp_num_threads)

    log_path = spec.out_dir.parent / "physicell.log"
    started = time.time()
    with log_path.open("w") as log:
        proc = subprocess.run(
            [str(binary.resolve()), str(spec.config_path.resolve())],
            cwd=str(physicell_root),
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            timeout=timeout_s,
        )
    manifest["wall_seconds"] = round(time.time() - started, 2)
    manifest["returncode"] = proc.returncode
    manifest["log"] = str(log_path)

    produced = len(list(spec.out_dir.glob("output*_cells.mat")))
    manifest["produced_timepoints"] = produced

    if proc.returncode != 0:
        manifest["status"] = "failed"
    elif produced < manifest["expected_timepoints"]:
        # This is the failure mode the manuscript hit: a trajectory that stopped
        # early. It is recorded as a failed run, not silently trimmed.
        manifest["status"] = "incomplete"
    else:
        manifest["status"] = "ok"

    (spec.out_dir.parent / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def plan(
    model: Model,
    run_root: Path,
    physicell_root: Path,
    arms: list[str] | None = None,
    n_replicates: int | None = None,
    apoptosis_factor: float | None = None,
) -> list[RunSpec]:
    arms = arms or list(model.arms)
    n = n_replicates if n_replicates is not None else model.run["n_replicates"]
    base = model.run["seed_base"]
    factor = (
        apoptosis_factor
        if apoptosis_factor is not None
        else model.raw["genotypes"]["IDHmut_TP53mut"]["apoptosis_rate"]["factor"]["value"]
    )
    specs = []
    for arm in arms:
        for i in range(n):
            # Same seed index across arms: replicate i of every arm shares a
            # seed, so arms are paired and can be compared within-seed.
            specs.append(
                build(model, arm, base + i, run_root, physicell_root, float(factor))
            )
    return specs


def run_all(
    specs: list[RunSpec],
    model: Model,
    physicell_root: Path,
    force: bool = False,
    on_progress=None,
) -> list[dict]:
    results = []
    for i, spec in enumerate(specs, 1):
        if not force and is_complete(spec, model):
            results.append({"arm": spec.arm, "seed": spec.seed, "status": "skipped"})
            if on_progress:
                on_progress(i, len(specs), spec, results[-1])
            continue
        m = execute(spec, model, physicell_root)
        results.append(m)
        if on_progress:
            on_progress(i, len(specs), spec, m)
    return results
