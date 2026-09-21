# Computational Optimization of pH/ROS Dual-Responsive Nanoparticle Systems for Targeted Atherosclerosis Therapy

This repository contains the computational model, release-benchmark data, convergence analyses, uncertainty results, parameter provenance tables, and source files supporting the manuscript **“Computational Optimization of pH/ROS Dual-Responsive Nanoparticle Systems for Targeted Atherosclerosis Therapy.”**

## Study scope

This is an entirely **in silico** study. The framework is intended for computational screening, hypothesis generation, and experimental prioritization. The integrated transport model is numerically verified and uncertainty-tested, but it has not been quantitatively validated against a formulation-matched in vivo biodistribution dataset. The results should not be interpreted as evidence of nanoparticle safety, therapeutic efficacy, clinically validated plaque targeting, or clinical feasibility.

## Methodological framework

The final framework links stimulus-responsive release with tissue-level exposure in a shared transport structure. It includes:

- dimensionally consistent pH- and ROS-response functions;
- nonresponsive, pH-only, ROS-only, dual-additive, and dual-interaction release mechanisms;
- nanoparticle-bound and free-drug states in blood, mononuclear phagocyte system (MPS), healthy arterial wall, plaque cap, and plaque core;
- Pareto screening using plaque AUC and off-target AUC as competing objectives;
- Pareto-balanced representative-design selection by normalized distance to the ideal objective point;
- paired quasi-Monte Carlo uncertainty propagation;
- Jansen/Saltelli global sensitivity analysis;
- candidate-count, uncertainty-sample, and Sobol-sample convergence analyses;
- internal held-out and leave-one-concentration-out benchmarking of the reduced ROS-release model.

The individual methods are established. The contribution of the work is their integration into one exposure-based framework that compares responsive mechanisms under the same transport assumptions and connects sensitivity results to experimentally measurable quantities.

## Repository structure

### Primary analysis
- `nanoparticle_optimization_submission_ready_v2.m` — primary MATLAB analysis used for the revised manuscript.
- `independent_reproduction_and_convergence.py` — independent reproduction of the final equations and convergence analyses.
- `nanoparticle_optimization.m` — compatibility entry point that redirects to the current primary MATLAB analysis.

### Parameter and provenance data
- `data/parameter_provenance_complete.csv` — complete machine-readable parameter table.
- `data/final_parameter_reference_audit.csv` — parameter/value/unit/uncertainty/reference audit used for the revised manuscript.

### Reproducibility and uncertainty outputs
- `results/convergence_candidate_screening.csv` — 1,000/2,000/3,000-candidate convergence results.
- `results/convergence_uncertainty.csv` — 500/1,000/2,500/5,000-draw uncertainty convergence results.
- `results/convergence_sobol.csv` — 512/1,024/2,048-base-sample sensitivity convergence results.
- `results/principal_uncertainty_effects.csv` — principal paired outcome intervals and differences.
- `results/paired_effect_sizes_python_check.csv` — paired relative effect sizes and probabilities.
- `results/external_circulation_benchmark.csv` — limited cross-formulation circulation comparison.
- `results/external_benchmark_model_timecourse.csv` — model time course used for the circulation comparison.

## Software requirements

### MATLAB workflow
MATLAB R2020b or later is recommended. The primary script uses base MATLAB functions, including `expm`, `fminsearch`, `table`, `writetable`, `ode45`, and `erfinv`. No Statistics, Optimization, or Global Optimization Toolbox is required.

### Independent reproduction
Python 3 with NumPy, SciPy, and pandas.

## Reproducibility settings

Primary publication settings:

- analysis horizon: 72 h;
- initial normalized dose: 1.0;
- candidate designs: 3,000 per mechanism;
- paired uncertainty draws: 5,000;
- Jansen/Saltelli base samples: 2,048;
- candidate Halton offset: 20;
- uncertainty Halton offset: 100;
- Sobol A offset: 50;
- Sobol B offset: 5,000;
- auxiliary pseudorandom seed: 42.

The principal sampling sequences are deterministic, so the reported screening, uncertainty, and sensitivity results do not depend on pseudorandom draws.

## Release model

The release model uses dimensionless activation functions:

```text
f_pH = 1 / [1 + 10^(n_pH*(pH-pH50))]
f_ROS = ROS^n_ROS / (K_ROS^n_ROS + ROS^n_ROS)
k_rel = k_leak + k_pH*f_pH + k_ROS*f_ROS + k_interaction*f_pH*f_ROS
```

Because `f_pH` and `f_ROS` are dimensionless, each release coefficient retains units of h^-1.

## Transport model

Nanoparticle-bound and free drug are represented in five biological compartments:

1. blood;
2. MPS;
3. healthy arterial wall;
4. plaque cap;
5. plaque core.

For a fixed parameter vector, the linear state system is solved by the matrix exponential and independently checked with `ode45`. The implementation verifies mass balance, nonnegative states, and agreement between the two numerical solutions.

The cap/core representation is a screening-level approximation of plaque heterogeneity. It does not resolve continuous spatial diffusion, extracellular-matrix structure, individual cells, receptor-mediated transport, dynamic protein-corona evolution, or intracellular pharmacodynamics.

## Exposure metrics

- **Plaque AUC:** integrated free-drug exposure in plaque cap + core.
- **Off-target AUC:** integrated free-drug exposure in blood + MPS + healthy arterial wall.
- **ESR:** `AUC_plaque / AUC_healthy_wall`.
- **SER:** `ESR_responsive / ESR_matched_nonresponsive`.
- **Exposure efficiency:** `AUC_plaque / (AUC_plaque + AUC_off_target)`.

All AUC values are normalized exposure proxies per unit initial payload.

## Multi-objective screening

Candidate designs are evaluated using two objectives:

1. maximize plaque free-drug AUC;
2. minimize off-target free-drug AUC.

No weighted composite score is used. Nondominated candidates define the Pareto front. The representative design is selected from the nondominated set by minimum normalized Euclidean distance to the ideal objective point. The manuscript refers to this procedure as **Pareto-balanced selection**.

## Convergence analyses

### Candidate screening
The selected ROS-only, dual-additive, and dual-interaction representative designs were unchanged between 2,000 and 3,000 candidates. The selected dual-interaction values were plaque AUC 1.5975, off-target AUC 3.5446, and SER 5.6081 at both sample sizes. The pH-only representative point was more sensitive to candidate count and is interpreted less uniquely.

### Paired uncertainty
Across 500, 1,000, 2,500, and 5,000 paired draws, the probability that dual plaque AUC exceeded ROS-only stayed between 96.9% and 97.1%. The probability of higher off-target exposure stayed between 82.8% and 83.4%, and the probability of higher SER stayed between 64.1% and 65.3%.

### Global sensitivity
Across 512, 1,024, and 2,048 base samples, plaque-cap pH and carrier pH50 remained the two dominant total-order SER determinants, followed by baseline leakage. At 2,048 samples, the total-order indices were approximately 0.4248 for cap pH, 0.4129 for pH50, and 0.1578 for baseline leakage.

## Principal paired uncertainty results

At 5,000 paired draws:

- dual plaque AUC median: 1.5169 (95% interval 0.7326-2.9588);
- ROS-only plaque AUC median: 1.0538 (0.4975-2.1239);
- paired plaque-AUC difference median: 0.4407 (-0.0125 to 1.3097);
- median relative plaque-AUC change: +43.2% (-1.1% to +129.8%);
- P(dual plaque AUC > ROS-only): 96.96%;
- dual off-target AUC median: 3.5045 (2.2781-5.8047);
- ROS-only off-target AUC median: 2.9210 (1.9959-4.1456);
- median relative off-target change: +18.3% (-10.3% to +92.5%);
- P(dual off-target AUC > ROS-only): 82.90%;
- dual SER median: 5.3968 (3.7672-7.2853);
- ROS-only SER median: 5.0065 (3.3598-7.0863);
- median relative SER change: +7.27% (-27.7% to +63.3%);
- P(dual SER > ROS-only): 64.14%.

These results show a consistent modeled advantage in absolute plaque exposure but substantially greater overlap in selectivity-based outcomes.

## Calibration, verification, and validation scope

The reduced ROS-release model was calibrated against digitized release data from Yang et al. A complete 5 mM H2O2 curve was withheld for the primary internal benchmark, and leave-one-concentration-out analysis was also performed. A published pH-responsive dataset was used only as a cross-formulation response-shape check.

The full transport model was numerically verified but is not independently validated against a formulation-matched in vivo biodistribution dataset. A limited external circulation comparison is included as context rather than validation. The selected model predicts a shorter blood-bound nanoparticle half-life than reported for published comparator systems, which is treated as evidence of remaining model uncertainty rather than as a successful fit.

## Reproducing the analysis

1. Clone or download this repository.
2. Open MATLAB in the repository root.
3. Run `nanoparticle_optimization_submission_ready_v2.m` using the `publication` profile.
4. Preserve the generated `publication_outputs` directory and record the MATLAB version.
5. Run `independent_reproduction_and_convergence.py` in Python.
6. Compare the generated outputs with the CSV files in `results/`.

## Data interpretation

The model provides relative, normalized exposure predictions. It does not establish safety, efficacy, clinical feasibility, or validated plaque targeting. Formulation-specific release, blood stability, MPS uptake, plaque penetration, pharmacodynamic, toxicology, and biodistribution measurements are required before preclinical or clinical conclusions can be made.

## CRediT author contribution

**Ayur Anchan:** Conceptualization, methodology, software, validation, formal analysis, investigation, data curation, visualization, writing – original draft, writing – review and editing, and project administration.

## Generative-AI disclosure

Generative AI was **not used to generate experimental data, perform the reported numerical simulations, choose the reported numerical results, or determine the scientific conclusions**. AI-assisted tools were used during manuscript and repository preparation for language editing, organization, and scientific schematic development. All scientific content, equations, code, numerical outputs, references, and conclusions were reviewed and verified by the author.

## Citation

If using this repository, please cite the associated manuscript/preprint and the archived repository release corresponding to the version used.

## License and disclaimer

This repository is provided for research and educational use. The computational outputs are exploratory and should not be interpreted as clinical recommendations or validated therapeutic predictions.