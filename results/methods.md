## Methods (generated 2026-09-07)

### Simulation platform

Agent-based simulations were performed with PhysiCell version 1.14.2 (git commit dbd349925014), compiled from source and run from the command line on a local installation rather than through the NanoHub web interface. The unmodified `template` sample project was used; no custom C++ was written. All model specification is contained in the XML configuration and the initial-condition CSV emitted by the accompanying harness, both of which are archived with the run outputs. Runs executed with 2 OpenMP threads; because PhysiCell seeds one random-number stream per thread from the base seed, the thread count is part of the reproducibility specification and was held constant across all runs.

### Domain and discretisation

The domain spans -500 to 500 µm in x and y, discretised at 20 µm, giving a 50 × 50 mesh of 2500 voxels, restricted to the z = 0 plane. Total simulated time was 7920 min (5.50 days) with full output saved every 360 min. This yields 22 intervals and therefore **23 saved timepoints** including t = 0 and the endpoint. Diffusion, mechanics and phenotype timesteps were 0.01, 0.1 and 6 min respectively.

### Cell definitions and the doubling-time conversion

All agents carry the IDH1 R132H mutation. The comparison is TP53-wild-type versus TP53-mutant on that shared background; ATRX and CDKN2A/B status are held constant and are not represented, and CNS WHO grade 4 is assumed for the model rather than demonstrated by it.

Proliferation uses the PhysiCell Live cycle model (code 5), a single phase with an exponentially distributed duration and transition rate *r* (min⁻¹). Because each transition yields one additional cell and apoptosis removes cells at hazard *d*, the population satisfies dN/dt = (r − d)N, so the population doubling time is ln(2)/(r − d) — **not** the mean phase duration. Published doubling times were converted to cycle transition rates by inverting that relation, r = ln(2)/T_double + d, since a doubling time measured in culture already nets out cell loss. Setting the phase duration equal to the published doubling time, as is sometimes done, would encode a doubling of 3.40 days rather than the intended 4.90 days, an error of a factor of ln(2).

The TP53-wild-type agent was given a population doubling time of 4.9 days, the mean reported for patient-derived IDH1-mutant glioma cultures, yielding r = 1.5140e-04 min⁻¹ (mean phase duration 6605 min) at the default apoptosis hazard of 5.3167e-05 min⁻¹. The TP53-mutant agent was given the **same cycle transition rate** and an apoptosis hazard reduced 10-fold to 5.3167e-06 min⁻¹. Its population doubling time of 3.30 days is therefore an output of those two settings, not an input.

The magnitude of the apoptosis reduction is an assumption. No study was identified that measures a basal per-minute apoptosis hazard for a named TP53 substitution in IDH1-mutant glioma, so this parameter is reported as unsourced and is varied across a sensitivity sweep rather than asserted; results are presented as a function of it. No variant-specific mechanism is claimed.

Motility was disabled and mechanical parameters were identical across all cell definitions, so the model contains no genotype-specific spatial phenotype and spatial differences between populations reflect only their seeding positions.

### Replication and analysis

Each condition was run as 10 independent replicates differing only in random seed, with seeds recorded in each run's manifest. Replicate *i* of every arm shares a seed, so arms are compared within seed. **The unit of analysis is the run, not the timepoint**: each run was reduced to replicate-level outcomes — the fitted log growth rate (primary), endpoint count, log₂ fold change, area under the log-count curve, final proportion, Shannon diversity change and final radius of gyration — before any comparison was made. Successive timepoints within a run are not treated as independent observations.

Effects are reported as within-seed paired mean differences with percentile bootstrap 95% confidence intervals, corroborated by a mixed-effects model with seed as a random intercept. A replicate-count convergence check reports the relative CI half-width as a function of *n*. Monoculture arms were run for every genotype and checked against the closed-form expectation N(t) = N₀·exp((r−d)t). Position-swapped and well-mixed layouts, and a label control in which two identically parameterised populations are seeded and compared, bound the false-positive floor of the analysis; competition coefficients are defined as the difference between a genotype's growth rate in co-culture and in monoculture at the same seed. Runs that terminated before producing all 23 timepoints were rerun, not trimmed.

### Sources cited by parameters

- **physicell_default** — PhysiCell 1.14.2 default cell definition, as distributed in sample_projects/template/config/PhysiCell_settings.xml.
    - supports: apoptosis death_rate 5.31667e-05 1/min
    - supports: oxygen diffusion 100000 micron^2/min, decay 0.1 1/min
    - supports: volume, mechanics, and necrosis fluid/biomass change rates
    - caveat: These are software defaults, not measurements for astrocytoma. They are reported as 'default' in Table 1 and must not be described as literature-derived.
- **tp53_apoptosis_UNSOURCED** — No source. Value is an assumption swept over a range.
    - caveat: Any claim that TP53 mutation reduces the apoptosis rate by a specific factor is NOT supported. Report results as conditional on the swept factor.
- **verheul2021** — Verheul C, Ntafoulis I, Kers TV, Hoogstrate Y, Mastroberardino PG, Barnhoorn S, et al. Generation, characterization, and drug sensitivities of 12 patient-derived IDH1-mutant glioma cell cultures. Neuro-Oncology Advances. 2021;3(1):vdab103. doi:10.1093/noajnl/vdab103
    - supports: IDH1-mutant patient-derived glioma cultures: mean population doubling time 4.9 days (range 1.8-16.4 days), n = 12 cultures.
    - supports: IDH-wildtype comparator cultures: mean population doubling time 2.6 days (range 1.6-3.6 days).
    - caveat: In vitro monolayer/spheroid doubling times. They are population doubling times, so they already net out cell loss; the conversion in convert.py accounts for this when solving for the PhysiCell cycle rate.
    - caveat: Wide range (1.8-16.4 d) means the mean is a weak point estimate. The range is carried into the sensitivity sweep rather than discarded.
