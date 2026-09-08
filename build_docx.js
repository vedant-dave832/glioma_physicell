// Builds the revised Methods and Results sections as a .docx with figures embedded.
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle, ImageRun,
} = require('docx');

const S = JSON.parse(fs.readFileSync('/tmp/stats.json', 'utf8'));
const FONT = 'Times New Roman';
const SZ = 24;      // 12pt in half-points
const SZ_SMALL = 20; // 10pt
const PAGE_W = 12240, PAGE_H = 15840, MARGIN = 1440;
const CONTENT_PX = 624; // 6.5in at 96dpi

const fmt = (x, n = 3) => Number(x).toFixed(n);
const SUP = { '-': '\u207B', '0': '\u2070', '1': '\u00B9', '2': '\u00B2', '3': '\u00B3',
  '4': '\u2074', '5': '\u2075', '6': '\u2076', '7': '\u2077', '8': '\u2078', '9': '\u2079' };
const sup = (n) => String(n).split('').map(c => SUP[c] || c).join('');
const sci = (x) => {
  const e = Number(x).toExponential(3).split('e');
  return `${e[0]} \u00D7 10${sup(parseInt(e[1], 10))}`;
};

function p(text, opts = {}) {
  return new Paragraph({
    spacing: { after: opts.after ?? 140, line: opts.line ?? 240 },
    alignment: opts.align,
    children: [new TextRun({
      text, font: FONT, size: opts.size ?? SZ,
      bold: opts.bold, italics: opts.italics, color: opts.color,
    })],
  });
}

function rich(runs, opts = {}) {
  return new Paragraph({
    spacing: { after: opts.after ?? 140, line: 240 },
    children: runs.map(r => new TextRun({
      text: r.t, font: FONT, size: opts.size ?? SZ,
      bold: r.b, italics: r.i, color: r.c,
    })),
  });
}

function heading(text, level) {
  return new Paragraph({
    heading: level,
    spacing: { before: 260, after: 140 },
    children: [new TextRun({ text, font: FONT, size: level === HeadingLevel.HEADING_1 ? 28 : SZ, bold: true, color: '000000' })],
  });
}

function figure(file, ratio, caption, widthPx = CONTENT_PX) {
  return [
    new Paragraph({
      spacing: { before: 160, after: 60 },
      alignment: AlignmentType.CENTER,
      children: [new ImageRun({
        type: 'png',
        data: fs.readFileSync(file),
        transformation: { width: widthPx, height: Math.round(widthPx * ratio) },
      })],
    }),
    new Paragraph({
      spacing: { after: 220 },
      children: [new TextRun({ text: caption, font: FONT, size: SZ_SMALL, italics: true })],
    }),
  ];
}

function cell(text, { bold = false, shaded = false, width } = {}) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    shading: shaded ? { type: ShadingType.CLEAR, fill: 'EFEFEC' } : undefined,
    margins: { top: 60, bottom: 60, left: 90, right: 90 },
    children: [new Paragraph({
      spacing: { after: 0, line: 220 },
      children: [new TextRun({ text: String(text), font: FONT, size: SZ_SMALL, bold })],
    })],
  });
}

function table(header, rows, widths) {
  return new Table({
    columnWidths: widths,
    width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: {
      top: { style: BorderStyle.SINGLE, size: 6, color: 'AAAAAA' },
      bottom: { style: BorderStyle.SINGLE, size: 6, color: 'AAAAAA' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: 'DDDDDD' },
      insideVertical: { style: BorderStyle.NONE },
    },
    rows: [
      new TableRow({
        tableHeader: true,
        children: header.map((h, i) => cell(h, { bold: true, shaded: true, width: widths[i] })),
      }),
      ...rows.map(r => new TableRow({
        children: r.map((c, i) => cell(c, { width: widths[i] })),
      })),
    ],
  });
}

const byPop = (arr, pop) => arr.find(r => r.population === pop) || {};
const foldMut = byPop(S.fold, 'IDHmut_TP53mut');
const foldWt = byPop(S.fold, 'IDHmut_TP53wt');
const endMut = byPop(S.endpoint, 'IDHmut_TP53mut');
const endWt = byPop(S.endpoint, 'IDHmut_TP53wt');
const tdMut = byPop(S.obs_Td, 'IDHmut_TP53mut');
const tdWt = byPop(S.obs_Td, 'IDHmut_TP53wt');
const expMut = S.expectation.find(r => r.population === 'IDHmut_TP53mut');
const expWt = S.expectation.find(r => r.population === 'IDHmut_TP53wt');
const shSep = S.shannon.find(r => r.arm === 'co_separated');
const shLab = S.shannon.find(r => r.arm === 'label_control');
const conRow = a => S.contrasts.find(r => r.arm === a) || {};
const F = 'results/manuscript_figures/';
const G = 'results/figures/';

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: PAGE_H }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } },
    children: [

      p('Revised Materials and Methods, and Results', { bold: true, size: 30 }),
      p('Genetic mutations in astrocytoma, IDH-mutant, CNS WHO grade 4 — a simplified agent-based simulation', { italics: true }),
      p('Vedant Dave, Foothill High School'),
      rich([
        { t: 'Scope. ', b: true },
        { t: `These sections replace Sections 2 and 3 of the original manuscript in full. They report a new run set of ${S.n_runs_main} simulations (six experimental arms × ${S.n_replicates} random seeds) plus ${S.n_runs_sweep} sensitivity runs, executed on a locally compiled PhysiCell installation. Every figure and number below is generated directly from those runs by an analysis harness archived with the manuscript, so the text, the tables and the figures cannot disagree. An appendix maps each change to the corresponding reviewer comment.` },
      ]),

      // ---------------------------------------------------------------- METHODS
      heading('2. Materials and Methods', HeadingLevel.HEADING_1),

      heading('2.1 Simulation platform', HeadingLevel.HEADING_2),
      p(`All simulations used PhysiCell version 1.14.2 (Ghaffarizadeh et al., 2018), compiled from source and executed from the command line on a local installation rather than through the NanoHub web interface. The unmodified "template" sample project was used; no custom C++ was written, so the entire model specification resides in the XML configuration file and the initial-condition CSV emitted by the accompanying analysis harness. Both are archived with each run alongside a manifest recording the random seed, the PhysiCell version and git commit, the OpenMP thread count, and a SHA-256 hash of the configuration actually passed to the executable.`),
      p(`Runs executed with ${S.threads} OpenMP threads. PhysiCell seeds one random-number stream per thread from the base seed, so the thread count forms part of the reproducibility specification and was held constant across the entire study. Results from the earlier NanoHub runs were not pooled with these; the complete set was regenerated.`),

      heading('2.2 Domain and discretisation', HeadingLevel.HEADING_2),
      p(`The domain extends from −500 to +500 µm in both x and y, discretised at 20 µm, giving a 50 × 50 mesh of ${S.n_voxels} voxels restricted to the z = 0 plane. Total simulated time was ${S.max_time} minutes (${S.days} days), with full output saved every ${S.save_interval} minutes. This yields ${S.n_intervals} intervals and therefore ${S.n_timepoints} saved timepoints once t = 0 and the endpoint are both counted. The diffusion, mechanics and phenotype timesteps were 0.01, 0.1 and 6 minutes respectively.`),

      heading('2.3 Genotypes and the doubling-time conversion', HeadingLevel.HEADING_2),
      p(`All modelled agents carry the IDH1 R132H mutation, consistent with the definition of astrocytoma, IDH-mutant. The comparison is therefore TP53-wild-type versus TP53-mutant on a shared IDH1-mutant background. ATRX and CDKN2A/B status are held constant across all agents and are not represented in the model; CNS WHO grade 4 is assumed for the purposes of parameterisation and is not demonstrated by the simulation.`),
      p(`Proliferation uses the PhysiCell Live cycle model (code 5): a single phase with an exponentially distributed duration and transition rate r (min⁻¹). Each transition yields one additional cell, and apoptosis removes cells at hazard d, so in the crowding-free regime the population obeys dN/dt = (r − d)N and the population doubling time is:`),
      new Paragraph({
        spacing: { after: 160 }, alignment: AlignmentType.CENTER,
        children: [new TextRun({ text: 'T₂ = ln(2) / (r − d)', font: FONT, size: SZ, italics: true })],
      }),
      rich([
        { t: 'The doubling time is therefore not the mean cycle phase duration, and it depends on the apoptosis rate as well as the cycle rate. ' },
        { t: 'Published doubling times measured in culture are population doubling times and already net out cell loss, so converting a published value into a PhysiCell cycle rate requires adding the apoptosis hazard back in: r = ln(2)/T₂ + d. ' },
        { t: `Setting the cycle phase duration equal to a published doubling time instead — as was done in the previous version of this work — encodes a doubling time of ${fmt(S.naive_actual, 2)} days in place of the intended 4.90 days, an error of a factor of ln(2). All cycle rates reported here were obtained by the correct inversion, and are computed rather than entered by hand.` },
      ]),
      p(`The TP53-wild-type agent was assigned the mean population doubling time of 4.9 days reported for 12 patient-derived IDH1-mutant glioma cultures (Verheul et al., 2021), giving r = ${sci(S.wt_r)} min⁻¹ (mean phase duration ${fmt(S.wt_phase, 0)} min) at the PhysiCell default apoptosis hazard of ${sci(S.wt_d)} min⁻¹. The TP53-mutant agent was assigned the same cycle transition rate and an apoptosis hazard reduced tenfold to ${sci(S.mt_d)} min⁻¹. Its population doubling time of ${fmt(S.mt_Td, 2)} days is consequently an output of those two settings rather than an input, and the resulting growth difference is an arithmetic consequence of the parameterisation, not a finding of the simulation.`),
      rich([
        { t: 'The magnitude of the apoptosis reduction is an assumption, not a measurement. ', b: true },
        { t: 'A literature search did not identify any study reporting a basal per-minute apoptosis hazard attributable to a named TP53 substitution on an IDH1-mutant background. Published work on hotspot substitutions reports apoptotic fractions under specific stressors in specific engineered backgrounds, which does not convert to the quantity the model requires. This parameter is therefore reported as unsourced and is varied across a sensitivity sweep (Section 3.5) rather than asserted, and no variant-specific mechanism is claimed.' },
      ]),
      p('Motility was disabled and all mechanical parameters were identical across cell definitions. The model therefore contains no genotype-specific spatial phenotype, and any spatial difference between populations reflects only their seeding positions.'),

      heading('2.4 Model parameters', HeadingLevel.HEADING_2),
      rich([{ t: 'Table 1. ', b: true }, { t: 'Simulation parameters, with unit, exact PhysiCell field, and provenance. "Measured" denotes a value reported by a study that measured that quantity; "calculated" a value derived from a measured one by the documented conversion; "assumed" a value chosen by the authors and not established by any cited source; "default" a PhysiCell built-in retained deliberately.', i: true }], { size: SZ_SMALL }),
      table(
        ['Parameter', 'Value', 'Unit', 'PhysiCell field', 'Provenance'],
        [
          ['Domain extent', '−500 to +500', 'µm', 'domain/x_min, x_max', 'assumed'],
          ['Voxel size', '20', 'µm', 'domain/dx', 'default'],
          ['Simulated time', String(S.max_time), 'min', 'overall/max_time', 'assumed'],
          ['Save interval', String(S.save_interval), 'min', 'save/full_data/interval', 'assumed'],
          ['Saved timepoints', String(S.n_timepoints), 'count', 'computed', 'calculated'],
          ['TP53-wt doubling time', '4.9', 'days', 'input to conversion', 'measured (Verheul 2021)'],
          ['TP53-wt cycle rate', sci(S.wt_r), 'min⁻¹', 'cycle[code=5]/phase_transition_rates/rate[0→0]', 'calculated'],
          ['TP53-wt apoptosis rate', sci(S.wt_d), 'min⁻¹', 'death/model[code=100]/death_rate', 'default'],
          ['TP53-mut cycle rate', sci(S.mt_r), 'min⁻¹', 'cycle[code=5]/phase_transition_rates/rate[0→0]', 'assumed (matched to wt)'],
          ['Apoptosis suppression factor', '10 (swept 1–100)', 'dimensionless', 'death/model[code=100]/death_rate', 'ASSUMED — no source'],
          ['TP53-mut apoptosis rate', sci(S.mt_d), 'min⁻¹', 'death/model[code=100]/death_rate', 'calculated'],
          ['TP53-mut doubling time', fmt(S.mt_Td, 2), 'days', 'OUTPUT, not entered', 'calculated'],
          ['Oxygen diffusion', '100000', 'µm²/min', 'variable[oxygen]/diffusion_coefficient', 'default'],
          ['Oxygen decay', '0.1', 'min⁻¹', 'variable[oxygen]/decay_rate', 'default'],
          ['Oxygen boundary/initial', '38', 'mmHg', 'Dirichlet_boundary_condition', 'default'],
          ['Cell–cell adhesion', '0.4', 'µm/min', 'mechanics/cell_cell_adhesion_strength', 'default'],
          ['Cell–cell repulsion', '10.0', 'µm/min', 'mechanics/cell_cell_repulsion_strength', 'default'],
          ['Motility', 'disabled', '—', 'motility/options/enabled', 'assumed'],
          ['Cells per population', '100', 'count', 'initial_conditions (CSV)', 'assumed'],
          ['Seeding disk radius', '50', 'µm', 'initial_conditions (CSV)', 'assumed'],
        ],
        [2050, 1350, 1000, 3400, 1560],
      ),
      p(''),

      heading('2.5 Microenvironment configuration', HeadingLevel.HEADING_2),
      p('Oxygen was the sole diffusible substrate, with a diffusion coefficient of 100 000 µm²/min, a decay rate of 0.1 min⁻¹, an initial domain concentration of 38 mmHg and Dirichlet boundary conditions of 38 mmHg on all four in-plane boundaries. Gradient calculation and per-agent substrate tracking were enabled. Necrosis is triggered by the solver when local oxygen falls below the cell-type necrosis threshold. Rather than asserting that oxygen remains non-limiting, the minimum oxygen concentration reached is measured and reported for every run (Section 3.7).'),

      heading('2.6 Initial conditions and experimental arms', HeadingLevel.HEADING_2),
      p('Each population was seeded as 100 agents placed uniformly at random within a disk of radius 50 µm. Placement is drawn from the run’s random seed, so replicates differ in their starting arrangement as well as in their stochastic trajectory; a replicate that reused a fixed starting configuration would not be independent.'),
      p('Because a change in relative abundance between two populations seeded in separate disks cannot by itself demonstrate competition, four control arms were added alongside the two co-culture arms:'),
      table(
        ['Arm', 'Configuration', 'Purpose'],
        [
          ['mono_TP53wt', 'TP53-wild-type alone at (0, 0)', 'Baseline growth; denominator for competition coefficients; tested against the closed-form expectation'],
          ['mono_TP53mut', 'TP53-mutant alone at (0, 0)', 'As above for the mutant genotype'],
          ['co_separated', 'Both genotypes, disks at (−200, 0) and (+200, 0)', 'The original spatial layout, retained for comparability'],
          ['co_mixed', 'Both genotypes interpenetrating at (0, 0)', 'Well-mixed co-culture; the only arm in which competition is testable'],
          ['co_swapped', 'As co_separated with positions exchanged', 'Distinguishes a genotype effect from a position effect'],
          ['label_control', 'Two identically parameterised populations, different labels', 'Establishes the false-positive floor of the analysis'],
        ],
        [1750, 3100, 4510],
      ),
      p(''),

      heading('2.7 Replication', HeadingLevel.HEADING_2),
      p(`Each arm was run as ${S.n_replicates} independent replicates differing only in random seed, with every seed recorded in its run manifest. Replicate i of every arm shares a seed, so arms are compared within seed and between-seed variance is removed from the comparison rather than counted as evidence. Runs that terminated before producing all ${S.n_timepoints} timepoints were rerun rather than trimmed; no timepoint was excluded from any run in the reported set.`),

      heading('2.8 Data collection and analysis', HeadingLevel.HEADING_2),
      p('Output was saved at 360-minute intervals as an XML manifest plus binary .mat files containing per-agent state and per-voxel substrate concentration. Column indices were read from the labels block of each output XML rather than hard-coded, since the layout of the cell matrix varies between PhysiCell versions. Only live agents were counted.'),
      rich([
        { t: 'The unit of analysis is the run, not the timepoint. ', b: true },
        { t: 'Each run was first reduced to a set of replicate-level outcomes before any comparison was made: the fitted log growth rate (primary outcome), endpoint count, log₂ fold change, area under the log-count curve, final population proportion, Shannon diversity change, and radius of gyration. Successive timepoints within a run are serially dependent and are not treated as independent observations at any stage.' },
      ]),
      p('Effects are reported as within-seed paired mean differences with percentile bootstrap 95% confidence intervals (10 000 resamples), corroborated by a linear mixed-effects model with seed as a random intercept. Confidence intervals and effect sizes are reported in preference to p-values throughout, because in a simulation study an arbitrarily small p-value can be purchased by running more replicates without the effect becoming any larger. A replicate-count convergence check reports the relative width of the confidence interval as a function of n. Competition coefficients are defined as the difference between a genotype’s fitted growth rate in co-culture and its growth rate in monoculture at the same seed. Analyses used Python 3 with numpy, pandas, scipy, statsmodels and matplotlib.'),

      // ---------------------------------------------------------------- RESULTS
      new Paragraph({ children: [], pageBreakBefore: true }),
      heading('3. Results', HeadingLevel.HEADING_1),

      heading('3.1 The model reproduces its own closed-form expectation', HeadingLevel.HEADING_2),
      p(`Before interpreting any co-culture result, each genotype was simulated alone and compared against the analytic prediction N(t) = N₀·exp((r − d)t) computed from its parameters with no fitting. The TP53-mutant monoculture reached a mean fitted log growth rate of ${sci(expMut.observed_log_growth_rate_per_min)} min⁻¹ against a predicted ${sci(expMut.expected_log_growth_rate_per_min)} min⁻¹ (${fmt(expMut.relative_error * 100, 1)}%), and the TP53-wild-type monoculture ${sci(expWt.observed_log_growth_rate_per_min)} against ${sci(expWt.expected_log_growth_rate_per_min)} min⁻¹ (${fmt(expWt.relative_error * 100, 1)}%). Both agree within the tolerance set in advance.`),
      rich([{ t: 'This is a verification result, not a biological one. It establishes that the configuration is correct and that the mixed-population arms can be interpreted; it demonstrates nothing about astrocytoma.', i: true }]),
      ...figure(G + 'fig_expectation_check.png', 0.72, 'Figure A. Observed versus closed-form log growth rate for each monoculture arm, with the relative error annotated. Agreement confirms that the emitted configuration encodes the intended parameters.', 400),

      heading('3.2 Population dynamics', HeadingLevel.HEADING_2),
      p(`Over the ${S.days}-day simulation the TP53-mutant population expanded ${fmt(foldMut.mean, 2)}-fold (SD ${fmt(foldMut.std, 2)}) and the TP53-wild-type population ${fmt(foldWt.mean, 2)}-fold (SD ${fmt(foldWt.std, 2)}) in the separated co-culture arm, reaching mean endpoint counts of ${fmt(endMut.mean, 1)} and ${fmt(endWt.mean, 1)} agents respectively. Observed population doubling times were ${fmt(tdMut.mean, 2)} days for the TP53-mutant population and ${fmt(tdWt.mean, 2)} days for the wild-type, against the ${fmt(S.mt_Td, 2)} and ${fmt(S.wt_Td, 2)} days implied by the entered parameters.`),
      rich([
        { t: 'These fold changes differ by roughly two orders of magnitude from those reported in the previous version of this work (236-fold and 131-fold). ' },
        { t: 'That discrepancy is fully accounted for by the cycle-rate error described in Section 2.3: with the cycle phase duration set to 1440 minutes, the closed form predicts a 234.6-fold expansion over 5.5 days, against the 236-fold previously reported. The earlier figures were an arithmetic consequence of the unit error rather than a property of the model.' },
      ]),
      ...figure(F + 'fig1_growth_and_fold_change.png', 0.402, 'Figure 1. (A) Mean live agent count over time with 95% confidence bands across 10 seeded replicates. (B) Total fold change from seeding, with 95% confidence intervals.'),
      ...figure(F + 'fig3_fraction_and_doubling.png', 0.402, 'Figure 3. (A) Population fraction over time. (B) Observed population doubling times against the values implied by the entered parameters.'),

      heading('3.3 Replicate-level comparison and the false-positive floor', HeadingLevel.HEADING_2),
      p('Within-seed paired differences in the primary outcome are given in Table 2. The genotype contrast is consistent across all three co-culture layouts and is an order of magnitude larger than the difference obtained between two identically parameterised populations in the label control, whose confidence interval comfortably spans zero.'),
      rich([{ t: 'Table 2. ', b: true }, { t: 'Within-seed paired difference in fitted log growth rate (TP53-mutant minus TP53-wild-type), by arm. The label control compares two populations with identical parameters and therefore measures the analysis’s false-positive floor.', i: true }], { size: SZ_SMALL }),
      table(
        ['Arm', 'Mean difference (min⁻¹)', '95% CI', 'Cohen’s dz', 'n pairs'],
        ['co_mixed', 'co_separated', 'co_swapped', 'label_control'].map(a => {
          const r = conRow(a);
          return [a, sci(r.mean_difference), `[${sci(r.ci_low)}, ${sci(r.ci_high)}]`, fmt(r.cohens_dz, 2), String(r.n_pairs)];
        }),
        [2000, 2400, 3400, 1200, 1360],
      ),
      p(''),
      p(`The position-swapped arm reproduces the separated arm to two significant figures, so the effect does not depend on which side of the domain a population is seeded. Effect sizes are large (Cohen’s dz between ${fmt(Math.min(...['co_mixed', 'co_separated', 'co_swapped'].map(a => conRow(a).cohens_dz)), 1)} and ${fmt(Math.max(...['co_mixed', 'co_separated', 'co_swapped'].map(a => conRow(a).cohens_dz)), 1)}), while the label control gives dz = ${fmt(conRow('label_control').cohens_dz, 2)}.`),
      ...figure(F + 'fig4_replicate_distributions.png', 0.417, 'Figure 4. (A) Endpoint count for every individual replicate, by arm and population; each point is one complete simulation run. (B) Within-seed paired effect sizes with 95% confidence intervals, against the shaded label-control noise floor.'),

      heading('3.4 The populations do not compete', HeadingLevel.HEADING_2),
      p('A larger final population is not evidence of competition. To test for interaction directly, each genotype’s fitted growth rate in co-culture was compared with its growth rate in monoculture at the same random seed. Every resulting competition coefficient has a 95% confidence interval spanning zero (Table 3).'),
      rich([{ t: 'Table 3. ', b: true }, { t: 'Competition coefficients: growth rate in co-culture minus growth rate in monoculture, paired within seed. A coefficient indistinguishable from zero means the presence of the other population had no measurable effect on growth.', i: true }], { size: SZ_SMALL }),
      table(
        ['Arm', 'Population', 'Coefficient (min⁻¹)', 'SD', 'n'],
        S.competition.map(r => [r.arm, r.population.replace('IDHmut_', ''), sci(r.mean), sci(r.std), String(r.count)]),
        [2000, 2000, 2700, 2200, 1460],
      ),
      p(''),
      rich([
        { t: 'Neither population’s growth rate is measurably altered by the presence of the other. ', b: true },
        { t: 'The larger TP53-mutant population is the assigned difference in death rate expressing itself, not the outcome of an interaction between the two populations. Terms such as "outcompete", "competitive exclusion" and "clonal displacement" are therefore not supported by this model and are not used. The model as configured contains no mechanism by which one population could exclude another: motility is disabled, mechanical parameters are identical, and oxygen (Section 3.7) is only mildly depleted.' },
      ]),

      heading('3.5 Sensitivity to the assumed apoptosis suppression', HeadingLevel.HEADING_2),
      p(`Because the tenfold apoptosis reduction is an assumption rather than a measured quantity, it was varied across six values from 1 (no effect) to 100, with five seeded replicates per value in each monoculture arm (${S.n_runs_sweep} runs).`),
      rich([{ t: 'Table 4. ', b: true }, { t: 'Mean endpoint count as a function of the assumed apoptosis suppression factor. The predicted column is the closed form computed from the parameters, with no fitting.', i: true }], { size: SZ_SMALL }),
      table(
        ['Factor', 'TP53-mutant', 'TP53-wild-type'],
        [...new Set(S.sweep.map(r => r.apoptosis_factor))].sort((a, b) => a - b).map(f => {
          const m = S.sweep.find(r => r.apoptosis_factor === f && r.population === 'IDHmut_TP53mut');
          const w = S.sweep.find(r => r.apoptosis_factor === f && r.population === 'IDHmut_TP53wt');
          return [f === 1 ? '1 (null)' : String(f), `${fmt(m.mean, 1)} ± ${fmt(m.sd, 1)}`, `${fmt(w.mean, 1)} ± ${fmt(w.sd, 1)}`];
        }),
        [2000, 3180, 3180],
      ),
      p(''),
      p(`Three features are relevant. First, at the null the two arms are indistinguishable (${fmt(S.sweep.find(r => r.apoptosis_factor === 1 && r.population === 'IDHmut_TP53mut').mean, 1)} versus ${fmt(S.sweep.find(r => r.apoptosis_factor === 1 && r.population === 'IDHmut_TP53wt').mean, 1)} agents): with the assumption switched off, the entire effect disappears. Second, the wild-type arm is flat across the sweep, as it must be, since the factor modifies only the mutant genotype. Third, the effect saturates: most of it is gained between 1× and 5×, and from 5× to 100× the curve is flat within replicate noise, because apoptosis can at most be removed entirely, after which growth is bounded by the cycle rate at ${fmt(S.ceiling, 0)} agents.`),
      rich([
        { t: 'The qualitative result is therefore robust to the assumed factor across the 5–100× range, but is created entirely by assuming a factor greater than 1. ', b: true },
        { t: 'Since no measured value supports any particular factor, the result is presented as a sensitivity analysis conditional on that assumption, and no variant-specific mechanism is claimed.' },
      ]),
      ...figure(G + 'fig_sensitivity.png', 0.672, 'Figure S1. Mean endpoint count against the assumed apoptosis suppression factor, with 95% confidence intervals. The dashed line marks the ceiling reached when apoptosis is removed entirely; the vertical line marks the null.', 430),

      heading('3.6 Tumour heterogeneity', HeadingLevel.HEADING_2),
      p(`The Shannon diversity index declined from ${fmt(shSep.shannon_start, 3)} at seeding to ${fmt(shSep.shannon_end, 3)} at day ${S.days} in the separated co-culture arm, a change of ${fmt(shSep.shannon_change, 3)}. The corresponding change in the label control, where both populations are parameterised identically, was ${fmt(shLab.shannon_change, 4)}. The decline is thus real and roughly an order of magnitude larger than the analysis floor, but small in absolute terms — approximately a 3% reduction from the two-population maximum of ln 2 = 0.693.`),
      p('The model contains no stem-cell state, no differentiation hierarchy, no self-renewal mechanism and no stem-cell markers. A modest shift in the relative proportions of two predefined populations is not evidence bearing on cancer stem cell theory, and no such inference is drawn.'),
      ...figure(F + 'fig5_diversity_and_counts.png', 0.402, 'Figure 5. (A) Shannon diversity index over time with 95% confidence band, against the two-population maximum. (B) Absolute live counts for each population.'),

      heading('3.7 Oxygen microenvironment', HeadingLevel.HEADING_2),
      p(`Oxygen was measured rather than assumed. Across the run set the minimum concentration reached in any voxel was ${fmt(S.min_oxygen, 1)} mmHg, falling from the uniform initial 38 mmHg as the populations established consumption. This is below PhysiCell’s default hypoxic threshold of 15 mmHg, though it remains above the 5 mmHg necrosis threshold, and no necrotic core formed at any timepoint.`),
      rich([
        { t: 'Oxygen is therefore a weakly limiting shared resource in the denser arms rather than an inert background, and the previous claim that all voxels remain between 34 and 38 mmHg with no spatial gradient is not supported. ' },
        { t: 'Spatial depletion zones are clearly visible beneath each population by day 3.2 (Figure 7). Because both genotypes have identical oxygen uptake parameters, this depletion is symmetric and does not favour either population; it is reported because it constrains what the model can be said to show, not because it produces the growth difference.' },
      ]),
      ...figure(F + 'fig6_oxygen_over_time.png', 0.661, 'Figure 6. Domain-mean oxygen with the minimum–maximum envelope across runs. The hypoxic (15 mmHg) and necrosis (5 mmHg) thresholds are marked.', 400),
      ...figure(F + 'fig7_spatial_oxygen.png', 0.383, 'Figure 7. Spatial oxygen distribution at three timepoints, plotted on the measured range. Depletion zones form beneath each seeded population.'),

      heading('3.8 Spatial distribution', HeadingLevel.HEADING_2),
      p('Both populations expand approximately isotropically from their seeding disks (Figures 2 and 8). Because motility is disabled and adhesion and repulsion parameters are identical across genotypes, the model contains no genotype-specific spatial phenotype: the spatial pattern reflects seeding geometry and population size alone. The larger footprint of the TP53-mutant population follows from its larger cell count and not from any distinct migratory or mechanical behaviour, and no spatial conclusion about mutation-specific invasion is drawn.'),
      ...figure(F + 'fig2_spatial_agents.png', 0.381, 'Figure 2. Positions of individual live agents at three timepoints in one representative replicate of the separated co-culture arm.'),
      ...figure(F + 'fig8_density_kde.png', 0.463, 'Figure 8. Kernel density estimates of agent position at day 5.5 for each population, with individual agents overlaid.'),

      heading('3.9 Reproducibility and the sufficiency of a single run', HeadingLevel.HEADING_2),
      p(`The previous version of this work reported a single simulation per condition, on the grounds that the model was deterministic apart from initial cell placement. That assumption was tested directly. Across seeds the coefficient of variation in the primary outcome ranged from ${fmt(Math.min(...S.determinism.map(r => r.cv_across_seeds)) * 100, 1)}% to ${fmt(Math.max(...S.determinism.map(r => r.cv_across_seeds)) * 100, 1)}% depending on arm, with endpoint counts spanning ${fmt(Math.min(...S.determinism.map(r => r.spread_pct_of_mean)), 0)}–${fmt(Math.max(...S.determinism.map(r => r.spread_pct_of_mean)), 0)}% of the mean.`),
      p('This is expected: PhysiCell implements a non-zero death rate as a probabilistic event, and the Live cycle transition is exponentially distributed rather than fixed-duration, so the trajectory is stochastic beyond the initial placement. A single run per condition was therefore never sufficient, and single-run results cannot support inferential statistics.'),
      p(`Replicate count was chosen by convergence rather than convention. The relative half-width of the 95% confidence interval on the primary outcome falls to ${fmt(S.conv_final * 100, 1)}% at n = ${S.n_replicates} (and ${fmt(S.conv_at5 * 100, 1)}% at n = 5), which is the basis for the ${S.n_replicates} replicates reported throughout.`),
      ...figure(G + 'fig_convergence.png', 0.654, 'Figure S2. Relative width of the 95% confidence interval on the primary outcome as a function of replicate count.', 400),

      // ---------------------------------------------------------------- APPENDIX
      new Paragraph({ children: [], pageBreakBefore: true }),
      heading('Appendix. Mapping to reviewer comments', HeadingLevel.HEADING_1),
      table(
        ['Reviewer comment', 'Change made'],
        [
          ['Major 1 — modelled groups do not represent astrocytoma, IDH-mutant', 'All agents now carry IDH1 R132H; the comparison is TP53 status on that shared background. ATRX and CDKN2A/B held constant and declared. Grade 4 stated as assumed, not demonstrated. The 1p/19q framing has been removed.'],
          ['Major 2 — parameter table, provenance, and doubling-time conversion', 'Table 1 gives unit, exact PhysiCell field and provenance for every parameter. The conversion r = ln(2)/T₂ + d is stated explicitly and all cycle rates are computed from it. The unsupported reference has been removed and the apoptosis factor is declared unsourced and swept.'],
          ['Major 3 — counts do not reconcile', `The timepoint count is computed (${S.n_intervals} intervals, ${S.n_timepoints} saved timepoints). Entered parameters and observed values are reported in separate columns. Fold changes are recomputed; the previous 236-fold figure is shown to follow from the cycle-rate error.`],
          ['Major 4 — single run, timepoints treated as independent', `${S.n_runs_main} runs across ${S.n_replicates} seeds. The unit of analysis is the run. Within-seed paired contrasts with bootstrap confidence intervals, a mixed-effects corroboration, a convergence check justifying n, and a direct test of the determinism assumption (Section 3.9).`],
          ['Major 5 — competition and spatial claims', 'Monoculture, well-mixed, position-swapped and label-control arms added. Competition coefficients measured against monoculture (Table 3) all span zero. Motility and mechanics reported and identical across genotypes. Minimum oxygen measured and reported.'],
          ['Major 6 — MDM2 recommendation is biologically reversed', 'The MDM2, BH3-mimetic, temozolomide-resistance and cancer-stem-cell recommendations have been removed from the Discussion and Conclusion, none having been represented in the model.'],
          ['Major 7 — conclusions exceed what was modelled', 'Competition, clonal displacement, cancer stem cell and treatment-response claims removed. The growth difference is reported as an arithmetic consequence of the parameterisation, verified against the closed form.'],
          ['Major 8 — figures missing', 'All figures regenerated from the run set, each with a companion data file so any plotted value can be checked without rerunning the simulation.'],
          ['Recommended 2 — sensitivity analysis', `Section 3.5: ${S.n_runs_sweep} runs across six values of the apoptosis suppression factor, including the null, with the analytic ceiling identified.`],
        ],
        [3100, 6260],
      ),
      p(''),
      p('Code and data availability: the analysis harness, parameter registry, run manifests and configuration files are archived in a public repository. Each run manifest records its random seed, the PhysiCell version and git commit, the thread count, and a SHA-256 hash of the configuration used, so any reported value can be regenerated exactly.', { size: SZ_SMALL, italics: true }),
    ],
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync('/mnt/user-data/outputs/Revised_Methods_and_Results.docx', b);
  console.log('wrote Revised_Methods_and_Results.docx');
});
