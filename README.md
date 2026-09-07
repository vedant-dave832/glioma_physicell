# glioma_physicell

A PhysiCell replicate-run and analysis harness for the astrocytoma agent-based
model, built to answer the NHSJS reviewer report point by point.

The design principle throughout: **the parameter file, the simulation config,
and the manuscript's Table 1 are the same object.** They are generated from
`params/model.yaml`, so they cannot drift apart — which is where most of the
reviewer's major issue #2 came from.

---

## What it does

| Reviewer point | What the harness does about it |
|---|---|
| **Major 1** — the three modelled groups don't represent astrocytoma IDH-mutant | Every genotype carries IDH1 R132H. The comparison is TP53-wt vs TP53-mut on that shared background. ATRX and CDKN2A/B are held constant and declared. The IDH-wildtype arm is marked a reference, not a subclone. Grade 4 is declared as assumed. |
| **Major 2** — parameter table, provenance, doubling-time conversion | `glioma table` emits the table with unit / exact PhysiCell field / meaning / source / provenance for every parameter. The registry **refuses to load** if a `measured` parameter points at a non-primary source, or an unsourced assumption isn't swept. |
| **Major 3** — counts don't reconcile | The timepoint count is computed, never typed (7920 ÷ 360 → **23**, not 33). Fitted doubling times come from the data and live in different columns from the entered parameters. |
| **Major 4** — one run, timepoints treated as independent | Seeded replicates; each run reduced to one row *before* any statistics; within-seed paired contrasts with bootstrap CIs; mixed-effects corroboration; a convergence check that justifies *n*; a determinism report that tests the "one run was enough" claim. Incomplete runs are reported and rerun, never trimmed. |
| **Major 5** — competition and spatial claims | Monoculture, well-mixed, position-swapped and label-control arms. Competition coefficients measured against monoculture at the same seed. Motility and mechanics identical across genotypes and reported. Minimum oxygen reached is recorded. |
| **Major 8** — missing figures | Every figure is generated, with a companion CSV of the exact values plotted. |
| **Recommended 2** — sensitivity analysis | The apoptosis suppression factor is a swept axis, not an asserted value. |

### The correction at the centre of it

PhysiCell's Live cycle model gives population doubling time

```
T_double = ln(2) / (r - d)
```

where `r` is the cycle transition rate and `d` the apoptosis hazard — so the
doubling time is **not** the mean phase duration, and it depends on the death
rate too. Setting the phase duration to 7056 min because the published doubling
time is 4.9 days actually encodes a **3.40-day** doubling. Published doubling
times are population doubling times, so the conversion has to add the death rate
back in:

```
r = ln(2) / T_double + d
```

`convert.py` implements exactly that, and the test suite reproduces the
reviewer's own recomputed doubling times (0.782 d, 0.698 d, 6.49 d) to confirm
the arithmetic agrees with theirs.

---

## Install

**macOS** — one command, because there is a trap. PhysiCell needs OpenMP and
Apple's `g++` is clang, which ships without it; the build fails with an error
that doesn't clearly say so. The script installs Homebrew GCC, sets
`PHYSICELL_CPP`, builds PhysiCell, and verifies the harness:

```bash
bash scripts/setup_macos.sh
source .envrc
```

**Linux** — the system `g++` already has OpenMP:

```bash
git clone --depth 1 --branch 1.14.2 https://github.com/MathCancer/PhysiCell.git
cd PhysiCell && make template && make -j && cd ..
export PHYSICELL_ROOT=$PWD/PhysiCell
pip install -r requirements.txt
./bin/glioma check
```

No PhysiCell C++ is modified. The stock `template` sample project is used with
`number_of_cells` set to 0, so every agent comes from the generated `cells.csv`.

## Use

```bash
./bin/glioma check                       # validate the model; runs nothing
./bin/glioma run --replicates 10         # seeded replicates, all arms
./bin/glioma sweep                       # the apoptosis-factor sensitivity sweep
./bin/glioma analyze --out results       # outcomes, statistics, figures, Methods
./bin/glioma table --out results         # Table 1 + Methods on their own
```

`check` is the fast loop: it prints the derived cycle rates, the computed
timepoint count, and the endpoint counts each arm is expected to reach *by
construction*, before anything is simulated.

Runs are resumable — a completed run is skipped, so a long replicate set can be
worked through in chunks. Each run writes a `manifest.json` recording the seed,
the PhysiCell version and commit, the thread count, and a hash of the config
actually handed to the binary.

### Timing

Roughly 40 min per replicate at the 5.5-day default on 2 cores. Ten replicates
across six arms is a multi-hour job; run it in chunks or raise
`run.omp_num_threads` in `params/model.yaml` if you have the cores. **Changing
the thread count changes the results for a given seed** — PhysiCell seeds one
RNG stream per thread — so pick it once and leave it alone for a whole
experiment.

## Layout

```
params/model.yaml     single source of truth: every parameter with provenance
params/sources.yaml   citation registry; a 'measured' value needs a primary source
src/gliomasim/
  convert.py          doubling time <-> cycle rate; timepoint arithmetic
  registry.py         loads and validates params; resolves derived values
  build_config.py     emits PhysiCell XML + cells.csv per (arm, seed)
  run.py              replicate runner, manifests, resume
  collect.py          parses PhysiCell output into tidy frames
  outcomes.py         replicate-level outcomes; closed-form expectation check
  stats.py            paired contrasts, bootstrap CIs, convergence, determinism
  plots.py            figures, each with a companion CSV
  report.py           Table 1, Methods text, reviewer map
tests/                the arithmetic and bookkeeping, not PhysiCell itself
```

## What this harness will not do for you

It will not make the TP53 apoptosis factor sourced. No study was found measuring
a basal per-minute apoptosis hazard for a named TP53 substitution in IDH1-mutant
glioma, so the factor is declared unsourced and swept. If the sweep shows the
headline result is a function of that assumption, that is the finding, and the
manuscript should say so rather than claiming a mechanism.

It will also not make the growth difference between the TP53 arms a discovery.
The two arms share a cycle rate and differ only in apoptosis, so the endpoint
ratio is arithmetic — `glioma check` prints it before any simulation runs. What
the simulation adds is the effect of stochasticity, crowding and spatial
arrangement on top of that, and the label control tells you how much of any
apparent difference is noise.

---

## What is tracked in git

Tracked: the parameter registry, the source, the tests, and — for every run —
its `manifest.json` and the exact `config/PhysiCell_settings.xml` and
`config/cells.csv` handed to the binary. Those are small and they are what make
a result checkable.

Not tracked: the raw PhysiCell output frames, which are large and fully
regenerable from the tracked config plus the seed recorded in the manifest. Also
not tracked: PhysiCell itself, which is a pinned dependency rather than part of
this work.

## Citing

`CITATION.cff` carries the machine-readable form. Anything using this harness
should also cite PhysiCell (Ghaffarizadeh et al. 2018,
doi:10.1371/journal.pcbi.1005991) and, for the doubling times the cycle rates
are derived from, Verheul et al. 2021 (doi:10.1093/noajnl/vdab103).

## License

BSD 3-Clause — see `LICENSE`. This matches PhysiCell's own license, which the
harness depends on but does not include.
