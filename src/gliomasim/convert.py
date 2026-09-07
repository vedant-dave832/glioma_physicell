"""Doubling time <-> PhysiCell Live-cycle rate conversions.

This module exists because reviewer major issue #2 asks, explicitly, "how you
converted published doubling times into your PhysiCell settings", and reviewer
major issue #3 points out that the manuscript's fitted doubling times do not
agree with its endpoint counts or with its cycle settings.

The algebra, stated once, in one place.

PhysiCell's Live cycle model (code 5) has a single phase with one transition,
phase 0 -> phase 0, which is a cell division. With ``fixed_duration=false`` the
transition is a Poisson process with rate ``r`` (units 1/min), so the phase
duration is exponentially distributed with mean ``1/r``.

Each transition takes one cell to two cells: a net gain of one cell. Apoptosis
removes cells with hazard ``d`` (units 1/min). In the exponential regime, before
crowding or oxygen limitation bites:

    dN/dt = (r - d) N        =>       N(t) = N0 * exp((r - d) t)

so the POPULATION doubling time is

    T_double = ln(2) / (r - d)

Three consequences, each of which is a mistake the manuscript currently makes:

1. ``T_double`` is NOT the mean phase duration ``1/r``. It is ``ln(2)`` times it
   (when d = 0). Setting the phase duration to 7056 min because the published
   doubling time is 4.9 days encodes a 3.4-day doubling, not a 4.9-day one.

2. ``T_double`` depends on the apoptosis rate as well as the cycle rate. Two
   genotypes given the same cycle rate but different apoptosis rates have
   different population doubling times. So "identical division rates" and
   "different growth" are not in contradiction, and the Conclusion's claim of
   identical division rates has to be stated as identical *cycle transition
   rates*.

3. Published doubling times from culture are POPULATION doubling times: the
   measured cells were dying as well as dividing. So a published ``T_double``
   constrains ``r - d``, not ``r``. To set the XML you must add the apoptosis
   rate back in:

    r = ln(2) / T_double + d

That last line is the conversion. Everything else here is bookkeeping around it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

LN2 = math.log(2.0)
MIN_PER_DAY = 1440.0


def days_to_min(days: float) -> float:
    return days * MIN_PER_DAY


def min_to_days(minutes: float) -> float:
    return minutes / MIN_PER_DAY


@dataclass(frozen=True)
class CycleSetting:
    """A fully derived Live-cycle setting, with its inputs kept alongside it."""

    cycle_rate_per_min: float          # what goes in the XML
    mean_phase_duration_min: float     # 1 / cycle_rate, for human sanity checking
    apoptosis_rate_per_min: float      # the d used in the conversion
    target_doubling_time_min: float    # the published number we were matching
    net_growth_rate_per_min: float     # r - d

    @property
    def implied_doubling_time_min(self) -> float:
        """Round-trip check: what doubling time does this setting actually give?"""
        return doubling_time_from_rates(
            self.cycle_rate_per_min, self.apoptosis_rate_per_min
        )

    def as_dict(self) -> dict:
        return {
            "cycle_rate_per_min": self.cycle_rate_per_min,
            "mean_phase_duration_min": self.mean_phase_duration_min,
            "apoptosis_rate_per_min": self.apoptosis_rate_per_min,
            "target_doubling_time_min": self.target_doubling_time_min,
            "target_doubling_time_days": min_to_days(self.target_doubling_time_min),
            "net_growth_rate_per_min": self.net_growth_rate_per_min,
            "implied_doubling_time_days": min_to_days(self.implied_doubling_time_min),
        }


def doubling_time_to_cycle_rate(
    doubling_time_min: float, apoptosis_rate_per_min: float
) -> CycleSetting:
    """Solve ``r = ln(2)/T_double + d`` and return the whole derivation.

    Parameters
    ----------
    doubling_time_min
        Published POPULATION doubling time, in minutes.
    apoptosis_rate_per_min
        The apoptosis hazard ``d`` this genotype will be given in the XML.

    Raises
    ------
    ValueError
        If the requested doubling time is not achievable at this apoptosis rate
        (it always is for positive inputs, but a non-positive doubling time or a
        negative rate is a configuration error worth catching loudly).
    """
    if doubling_time_min <= 0:
        raise ValueError(f"doubling_time_min must be positive, got {doubling_time_min}")
    if apoptosis_rate_per_min < 0:
        raise ValueError(
            f"apoptosis_rate_per_min must be non-negative, got {apoptosis_rate_per_min}"
        )

    net_growth = LN2 / doubling_time_min
    cycle_rate = net_growth + apoptosis_rate_per_min

    return CycleSetting(
        cycle_rate_per_min=cycle_rate,
        mean_phase_duration_min=1.0 / cycle_rate,
        apoptosis_rate_per_min=apoptosis_rate_per_min,
        target_doubling_time_min=doubling_time_min,
        net_growth_rate_per_min=net_growth,
    )


def doubling_time_from_rates(
    cycle_rate_per_min: float, apoptosis_rate_per_min: float
) -> float:
    """``T_double = ln(2) / (r - d)``, in minutes. Inverse of the above."""
    net = cycle_rate_per_min - apoptosis_rate_per_min
    if net <= 0:
        return math.inf  # population shrinks or holds; no doubling time exists
    return LN2 / net


def observed_doubling_time(
    n_start: float, n_end: float, elapsed_min: float
) -> float:
    """Doubling time implied by two counts. This is a MEASUREMENT, not a setting.

    Reviewer major issue #3 turns on keeping these apart: what you entered into
    the model and what you observed after running it are different numbers and
    must be reported in different columns.
    """
    if n_start <= 0 or n_end <= 0:
        raise ValueError("counts must be positive to fit a doubling time")
    if n_end <= n_start:
        return math.inf
    return elapsed_min * LN2 / math.log(n_end / n_start)


def n_saved_timepoints(max_time_min: float, save_interval_min: float) -> int:
    """Number of saved outputs, counting t=0 and the endpoint.

    A 7920-minute run saved every 360 minutes gives 22 intervals and therefore
    23 saved timepoints. The manuscript reports 33. This function is why nothing
    downstream is allowed to state a timepoint count by hand.
    """
    if save_interval_min <= 0:
        raise ValueError("save_interval_min must be positive")
    n_intervals = max_time_min / save_interval_min
    if abs(n_intervals - round(n_intervals)) > 1e-9:
        raise ValueError(
            f"max_time ({max_time_min}) is not an exact multiple of save_interval "
            f"({save_interval_min}); the final save would not land on the endpoint"
        )
    return int(round(n_intervals)) + 1


def naive_phase_duration_error(doubling_time_min: float) -> dict:
    """Quantify the mistake of setting the phase duration equal to the doubling time.

    Kept as a callable rather than a comment so the Methods text can cite an
    actual number and the test suite can assert on it.
    """
    setting = doubling_time_to_cycle_rate(doubling_time_min, 0.0)
    naive_rate = 1.0 / doubling_time_min  # phase duration typed in as T_double
    naive_doubling = doubling_time_from_rates(naive_rate, 0.0)
    return {
        "intended_doubling_days": min_to_days(doubling_time_min),
        "naive_phase_duration_min": doubling_time_min,
        "naive_actual_doubling_days": min_to_days(naive_doubling),
        "correct_phase_duration_min": setting.mean_phase_duration_min,
        "ratio": naive_doubling / doubling_time_min,  # == ln(2)
    }
