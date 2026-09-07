"""Tests for the replicate-level reduction and the derived comparisons.

The competition-coefficient mapping shipped inverted once and returned an empty
frame, which reads exactly like "no competition" rather than like a bug. These
tests pin the shape of the answer so a silent empty result fails loudly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from gliomasim import outcomes  # noqa: E402


def _traj(arm: str, seed: int, pop: str, rate: float, n0: float = 100.0) -> pd.DataFrame:
    """A clean exponential trajectory with a known log growth rate."""
    t = np.arange(0, 7921, 360.0)
    return pd.DataFrame(
        {
            "arm": arm,
            "seed": seed,
            "apoptosis_factor": 10.0,
            "frame": np.arange(len(t)),
            "time_min": t,
            "time_days": t / 1440.0,
            "population": pop,
            "n_live": n0 * np.exp(rate * t),
            "n_dead": 0,
            "centroid_x": 0.0,
            "centroid_y": 0.0,
            "radius_of_gyration": 50.0,
            "min_oxygen_mmHg": 38.0,
        }
    )


@pytest.fixture
def reps() -> pd.DataFrame:
    frames = []
    for seed in (1, 2, 3):
        frames.append(_traj("mono_TP53wt", seed, "IDHmut_TP53wt", 1.0e-4))
        frames.append(_traj("mono_TP53mut", seed, "IDHmut_TP53mut", 1.5e-4))
        # in co-culture the wild-type grows more slowly than it does alone
        frames.append(_traj("co_mixed", seed, "IDHmut_TP53wt", 0.8e-4))
        frames.append(_traj("co_mixed", seed, "IDHmut_TP53mut", 1.5e-4))
    return outcomes.per_replicate(pd.concat(frames, ignore_index=True))


def test_fitted_rate_recovers_the_planted_rate(reps):
    wt = reps[(reps.arm == "mono_TP53wt")][outcomes.PRIMARY_OUTCOME]
    assert np.allclose(wt, 1.0e-4, rtol=1e-6)


def test_competition_coefficients_are_not_empty(reps):
    comp = outcomes.competition_coefficients(reps)
    assert not comp.empty, "monoculture lookup failed; empty here reads as 'no competition'"
    assert set(comp["population"]) == {"IDHmut_TP53wt", "IDHmut_TP53mut"}
    assert len(comp) == 6  # 2 populations x 3 seeds, co_mixed only


def test_competition_coefficient_sign_and_size(reps):
    comp = outcomes.competition_coefficients(reps)
    wt = comp[comp.population == "IDHmut_TP53wt"]["competition_coefficient"]
    mut = comp[comp.population == "IDHmut_TP53mut"]["competition_coefficient"]
    # suppressed in co-culture by exactly the planted 0.2e-4
    assert np.allclose(wt, -0.2e-4, atol=1e-9)
    # unaffected
    assert np.allclose(mut, 0.0, atol=1e-9)


def test_monoculture_arms_are_excluded_from_their_own_comparison(reps):
    comp = outcomes.competition_coefficients(reps)
    assert not comp["arm"].str.startswith("mono_").any()


def test_shannon_change_is_zero_for_a_single_population():
    traj = _traj("mono_TP53wt", 1, "IDHmut_TP53wt", 1e-4)
    arm = outcomes.per_replicate_arm(traj)
    assert arm["shannon_start"].iloc[0] == 0.0
    assert arm["shannon_change"].iloc[0] == 0.0


def test_expectation_check_flags_a_mismatch():
    """A monoculture that misses the closed form must not silently pass."""
    from gliomasim import registry

    model = registry.load(REPO / "params")
    # plant a rate far from what the parameters imply
    traj = _traj("mono_TP53wt", 1, "IDHmut_TP53wt", 5.0e-4)
    for seed in (2, 3):
        traj = pd.concat([traj, _traj("mono_TP53wt", seed, "IDHmut_TP53wt", 5.0e-4)])
    check = outcomes.expectation_check(outcomes.per_replicate(traj), model)
    assert not check["passes"].all()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
