"""Tests for the parts that would fail quietly.

Run with: python -m pytest tests -q   (or: python tests/test_harness.py)

These are not tests of PhysiCell. They are tests of the arithmetic and the
bookkeeping between the parameter file, the XML, and the manuscript -- which is
where the manuscript's errors actually came from.
"""

from __future__ import annotations

import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pytest  # noqa: E402

from gliomasim import build_config, convert, registry, report  # noqa: E402

PARAMS = REPO / "params"


# --- the conversion the reviewer asked to see -------------------------------


def test_doubling_time_round_trips():
    for days in (1.0, 2.6, 4.9, 16.4):
        for d in (0.0, 5.31667e-05, 5.31667e-06):
            s = convert.doubling_time_to_cycle_rate(convert.days_to_min(days), d)
            assert math.isclose(
                convert.min_to_days(s.implied_doubling_time_min), days, rel_tol=1e-12
            )


def test_doubling_time_is_not_the_phase_duration():
    """The ln(2) error, pinned down."""
    e = convert.naive_phase_duration_error(convert.days_to_min(4.9))
    assert math.isclose(e["ratio"], math.log(2), rel_tol=1e-12)
    assert math.isclose(e["naive_actual_doubling_days"], 4.9 * math.log(2), rel_tol=1e-9)


def test_apoptosis_raises_the_required_cycle_rate():
    """A higher death rate needs a faster cycle for the same doubling time."""
    slow = convert.doubling_time_to_cycle_rate(7056, 0.0)
    fast = convert.doubling_time_to_cycle_rate(7056, 1e-4)
    assert fast.cycle_rate_per_min > slow.cycle_rate_per_min


def test_timepoint_count_matches_the_reviewers_arithmetic():
    assert convert.n_saved_timepoints(7920, 360) == 23
    assert convert.n_saved_timepoints(7920, 360) != 33


def test_non_divisible_save_interval_is_rejected():
    with pytest.raises(ValueError):
        convert.n_saved_timepoints(7920, 500)


def test_observed_doubling_time_reproduces_the_reviewers_numbers():
    """The reviewer recomputed three doubling times; we must get the same."""
    days = 5.5 * 1440
    assert math.isclose(
        convert.min_to_days(convert.observed_doubling_time(100, 13100, days)),
        0.782,
        abs_tol=0.002,
    )
    assert math.isclose(
        convert.min_to_days(convert.observed_doubling_time(100, 23600, days)),
        0.698,
        abs_tol=0.002,
    )
    assert math.isclose(
        convert.min_to_days(convert.observed_doubling_time(100, 180, days)),
        6.49,
        abs_tol=0.02,
    )


# --- the registry's guard rails ---------------------------------------------


def test_model_loads_and_validates():
    m = registry.load(PARAMS)
    assert m.n_timepoints == 23
    assert m.n_voxels == 2500
    assert set(m.genotypes) >= {"IDHmut_TP53wt", "IDHmut_TP53mut"}


def test_tp53_arms_share_a_cycle_rate():
    """The two primary arms must differ only in apoptosis."""
    m = registry.load(PARAMS)
    wt, mt = m.genotypes["IDHmut_TP53wt"], m.genotypes["IDHmut_TP53mut"]
    assert math.isclose(
        wt.cycle.cycle_rate_per_min, mt.cycle.cycle_rate_per_min, rel_tol=1e-12
    )
    assert mt.apoptosis_rate_per_min < wt.apoptosis_rate_per_min


def test_apoptosis_factor_override_flows_through():
    null = registry.load(PARAMS, apoptosis_factor=1.0)
    wt, mt = null.genotypes["IDHmut_TP53wt"], null.genotypes["IDHmut_TP53mut"]
    assert math.isclose(wt.apoptosis_rate_per_min, mt.apoptosis_rate_per_min)
    # at the null the two arms must be indistinguishable by construction
    assert math.isclose(wt.doubling_time_days, mt.doubling_time_days, rel_tol=1e-9)


def test_measured_parameters_point_at_primary_sources():
    m = registry.load(PARAMS)
    for p in m.params:
        if p.provenance == "measured":
            assert m.sources[p.source]["type"] == "primary", p.name


def test_unsourced_assumptions_must_be_swept():
    m = registry.load(PARAMS)
    for p in m.params:
        src = m.sources.get(p.source or "", {})
        if p.provenance == "assumed" and src.get("type") == "none":
            assert p.sweep, f"{p.name} is unsourced and not swept"


# --- the emitted config -----------------------------------------------------


@pytest.fixture(scope="module")
def physicell_root() -> Path:
    import os

    root = os.environ.get("PHYSICELL_ROOT")
    if not root:
        pytest.skip("PHYSICELL_ROOT not set")
    return Path(root)


def test_emitted_xml_carries_the_derived_rates(tmp_path, physicell_root):
    m = registry.load(PARAMS)
    spec = build_config.build(m, "co_separated", 42, tmp_path, physicell_root, 10.0)
    root = ET.parse(spec.config_path).getroot()

    assert root.find("options/random_seed").text == "42"

    for pop in spec.populations:
        cd = root.find(f"cell_definitions/cell_definition[@name='{pop['label']}']")
        assert cd is not None
        rate = float(cd.find("phenotype/cycle/phase_transition_rates/rate").text)
        assert math.isclose(rate, pop["cycle_rate"], rel_tol=1e-9)
        death = float(cd.find("phenotype/death/model[@code='100']/death_rate").text)
        assert math.isclose(death, pop["apoptosis_rate"], rel_tol=1e-9)
        assert cd.find("phenotype/cycle").get("code") == "5"

    # no agents may be placed at random; every cell comes from the CSV
    assert root.find("user_parameters/number_of_cells").text == "0"
    assert root.find("initial_conditions/cell_positions").get("enabled") == "true"


def test_seeded_initial_conditions_differ_between_replicates(tmp_path, physicell_root):
    m = registry.load(PARAMS)
    a = build_config.build(m, "mono_TP53wt", 1, tmp_path / "a", physicell_root, 10.0)
    b = build_config.build(m, "mono_TP53wt", 2, tmp_path / "b", physicell_root, 10.0)
    c = build_config.build(m, "mono_TP53wt", 1, tmp_path / "c", physicell_root, 10.0)
    assert a.cells_csv_path.read_text() != b.cells_csv_path.read_text()
    assert a.cells_csv_path.read_text() == c.cells_csv_path.read_text()


def test_cells_csv_seeds_the_requested_number(tmp_path, physicell_root):
    m = registry.load(PARAMS)
    spec = build_config.build(m, "co_separated", 7, tmp_path, physicell_root, 10.0)
    lines = [ln for ln in spec.cells_csv_path.read_text().splitlines() if ln.strip()]
    assert len(lines) == 2 * m.seeding["n_cells_per_population"]


def test_label_control_has_two_identical_definitions(tmp_path, physicell_root):
    m = registry.load(PARAMS)
    spec = build_config.build(m, "label_control", 3, tmp_path, physicell_root, 10.0)
    a, b = spec.populations
    assert a["label"] != b["label"]
    assert math.isclose(a["cycle_rate"], b["cycle_rate"], rel_tol=1e-12)
    assert math.isclose(a["apoptosis_rate"], b["apoptosis_rate"], rel_tol=1e-12)


# --- the manuscript artefacts -----------------------------------------------


def test_table1_has_every_column_the_reviewer_listed():
    m = registry.load(PARAMS)
    t = report.table1(m)
    for col in ("Parameter", "Value", "Unit", "PhysiCell field",
                "Biological meaning", "Provenance", "Source"):
        assert col in t.columns
    assert set(t["Provenance"]) <= {"measured", "calculated", "assumed", "default"}


def test_methods_states_the_computed_timepoint_count():
    m = registry.load(PARAMS)
    text = report.methods(m, {"physicell_version": "1.14.2", "physicell_commit": "abc"}, 10)
    assert "23 saved timepoints" in text
    assert "33" not in text


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
