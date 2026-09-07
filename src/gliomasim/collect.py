"""Parse PhysiCell output into tidy per-timepoint records.

One row per (arm, seed, timepoint, population). Column labels are read from each
output XML rather than hard-coded, because the column layout of the cells matrix
changes between PhysiCell versions and a silently shifted column is the kind of
bug that produces a plausible-looking wrong answer.

Only LIVE cells are counted. The manuscript's counts do not say whether dead
agents were included; here it is explicit and recorded in the column name.
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

_OUT_RE = re.compile(r"output(\d+)\.xml$")


def _labels(xml_path: Path) -> dict[str, tuple[int, int]]:
    """Map label name -> (start index, size), read from the output itself."""
    root = ET.parse(xml_path).getroot()
    out: dict[str, tuple[int, int]] = {}
    for lab in root.findall(
        ".//cellular_information/cell_populations/cell_population/custom/"
        "simplified_data/labels/label"
    ):
        name = (lab.text or "").strip()
        idx, size = int(lab.get("index")), int(lab.get("size"))
        # 'elapsed_time_in_phase' appears twice; keep the first occurrence
        out.setdefault(name, (idx, size))
    return out


def _current_time(xml_path: Path) -> float:
    root = ET.parse(xml_path).getroot()
    node = root.find(".//current_time")
    return float(node.text) if node is not None else float("nan")


def _min_oxygen(me_path: Path) -> float:
    if not me_path.exists():
        return float("nan")
    m = loadmat(me_path)
    key = next(k for k in m if not k.startswith("__"))
    arr = m[key]
    # rows: x, y, z, voxel_volume, then one row per substrate
    if arr.shape[0] < 5:
        return float("nan")
    return float(np.min(arr[4, :]))


def read_run(run_dir: Path) -> pd.DataFrame:
    """Read one (arm, seed) run directory into a tidy frame."""
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"no manifest.json in {run_dir}; run is not trustworthy")
    manifest = json.loads(manifest_path.read_text())

    type_names = {p["type_id"]: p["label"] for p in manifest["populations"]}
    out_dir = run_dir / "output"
    xmls = sorted(p for p in out_dir.glob("output*.xml") if _OUT_RE.search(p.name))

    rows = []
    for xml_path in xmls:
        frame_idx = int(_OUT_RE.search(xml_path.name).group(1))
        mat_path = xml_path.with_name(xml_path.stem + "_cells.mat")
        if not mat_path.exists():
            continue
        labels = _labels(xml_path)
        cells = loadmat(mat_path)["cells"]

        t = _current_time(xml_path)
        ci, _ = labels["cell_type"]
        di, _ = labels["dead"]
        pi, _ = labels["position"]

        cell_type = cells[ci, :].astype(int)
        dead = cells[di, :].astype(int)
        x, y = cells[pi, :], cells[pi + 1, :]

        min_o2 = _min_oxygen(
            xml_path.with_name(xml_path.stem + "_microenvironment0.mat")
        )

        for tid, label in type_names.items():
            sel = (cell_type == tid) & (dead == 0)
            n = int(sel.sum())
            if n:
                px, py = x[sel], y[sel]
                cx, cy = px.mean(), py.mean()
                spread = float(np.sqrt(((px - cx) ** 2 + (py - cy) ** 2).mean()))
            else:
                cx = cy = spread = float("nan")
            rows.append(
                {
                    "arm": manifest["arm"],
                    "seed": manifest["seed"],
                    "apoptosis_factor": manifest["apoptosis_factor"],
                    "frame": frame_idx,
                    "time_min": t,
                    "time_days": t / 1440.0,
                    "population": label,
                    "n_live": n,
                    "n_dead": int(((cell_type == tid) & (dead == 1)).sum()),
                    "centroid_x": cx,
                    "centroid_y": cy,
                    "radius_of_gyration": spread,
                    "min_oxygen_mmHg": min_o2,
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"no parsable output frames in {out_dir}")
    return df.sort_values(["frame", "population"]).reset_index(drop=True)


def read_all(run_root: Path, exclude: tuple[str, ...] = ("sweep_factor_",)) -> pd.DataFrame:
    """Read every completed run under run_root.

    ``exclude`` drops runs whose path contains any of these fragments. It
    defaults to the sweep directories: those runs use a deliberately different
    apoptosis factor, so pooling them into the main replicate set would mix
    conditions and inflate the apparent variance. The sweep is read separately,
    with ``exclude=()``.
    """
    frames = []
    problems = []
    for manifest_path in sorted(run_root.rglob("manifest.json")):
        if any(frag in str(manifest_path) for frag in exclude):
            continue
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("status") not in ("ok",):
            problems.append(
                f"{manifest['arm']}/seed_{manifest['seed']}: status="
                f"{manifest.get('status')} "
                f"({manifest.get('produced_timepoints')}/"
                f"{manifest.get('expected_timepoints')} timepoints)"
            )
            continue
        frames.append(read_run(manifest_path.parent))
    if problems:
        # Surfaced, never silently dropped: an incomplete trajectory is the
        # exact failure the manuscript papered over by deleting a timepoint.
        print("Excluded runs (rerun these rather than trimming them):")
        for p in problems:
            print("  -", p)
    if not frames:
        raise ValueError(f"no completed runs under {run_root}")
    return pd.concat(frames, ignore_index=True)
