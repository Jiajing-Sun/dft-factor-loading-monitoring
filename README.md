# Self-Normalized Online Monitoring of Factor-Loading Instability

R replication materials for the paper by **Jiajing Sun, Yongmiao Hong and Meiting Zhu**,
prepared for *Econometric Theory*. This branch contains the theory experiments and corresponding replication draws. The earlier manuscript package is preserved on main at commit f200d38.

The main method is affine-equivariant adjusted-range-based self-normalized
monitoring (RSMS-AE, called RSMS in the main text). RSMS-LDL is the
coordinate-based construction in the online appendix. SSMS is the quadratic
self-normalized benchmark. **aersn is the R package name.** The historical
array key AERSN denotes RSMS-AE; all manuscript-facing tables and plots use RSMS.

## Run

Use R 4.5.1 (the replication environment) or a compatible R version >= 4.1.
From the repository root:

```sh
Rscript install_dependencies.R
Rscript run_all.R --quick
```

The quick command uses the supplied simulation draws, **recomputes** the
two-sample paired bootstrap and all numerical tables, and rebuilds the FHFA
panel and observed score paths from the frozen CSV. It does not regenerate the
Monte Carlo or FHFA null draws.

For all simulations and the empirical null reference from scratch:

```sh
Rscript run_all.R --full
```

The full run is substantially longer. DFT_CORES controls FHFA workers (default
2); Windows uses a serial loop. Every FHFA replicate sets its own seed, so the
worker count does not change the draws. DFT_RESULTS can select another output
directory; the default is results/. R_LIBS_USER can select a package library.
No author-specific file paths are needed.

## Flat layout

All R scripts are at the top level. There are three single-level folders:

- data/: the frozen public FHFA input.
- reference/: saved reference/calibration draws, evaluation stopping times
  and check values used for the paper.
- results/: generated tables, figures, simulation files and validation report.

| Script | Purpose |
| --- | --- |
| run_all.R | Execute the complete quick or full workflow |
| install_dependencies.R, setup.R | Install/check packages and initialize paths/RNG |
| monitoring_core.R | Exact q=2/q=4 increment-hull geometry and SSMS |
| rsms_monitor.R | Public-API RSMS-AE and RSMS-LDL monitoring functions |
| profile_run.R | q=2 study, profile index 1--5 |
| profile_confirm.R | q=4 study, profile index 1--2 |
| profile_confirm_summarize.R | Resample calibration and evaluation panels independently |
| fhfa_helpers.R | Read the balanced panel and fit/simulate the factor VAR(1) |
| fhfa_aersn_revision.R | FHFA calibration, independent evaluation and figures |
| build_tables.R | Turn stopping times into manuscript tables |
| validate.R | Gauge equivalence, affine invariance and saved-result comparisons |
| simulate_training.R | Sparse and diffuse training contamination, with flat and front-cycle variance profiles |
| simulate_local_power.R | Local Gaussian score law and fixed-alternative consistency with fitted PCA |
| simulate_delay.R | Jump and linear onset against the conditional Gaussian exit law |
| theory_table.R | Assemble the three experiments into main Table 2 |

Each script contains remarks explaining its statistical and computational
steps. Individual study commands, run from the root, are e.g.
Rscript profile_run.R 1 and Rscript profile_confirm.R 2.
Run the summarizer only after both q=4 profiles are available.

## Manuscript-to-output map

| Location | Output in results/ |
| --- | --- |
| Main Table 1: q=4 detection and restricted mean delay | three_methods.csv; three_method_rows.tex (RSMS/SSMS contrasts in table1_power_delay.csv) |
| Main Table 2: contamination, local power and delay | theory_table.tex; training_summary.csv; local_power_results.csv; delay_results.csv |
| Main Table 3: FHFA reference and alarms | table3_fhfa.csv |
| Figure 1: cumulative detection (three methods) | three_methods.pdf |
| Figure 2: bandwidth sensitivity under serial dependence | front_cycle_bandwidth.pdf |
| Figure 3: average house-price growth | fhfa_J1_cross_sectional_avg_growth.pdf |
| Figure 4: observed monitoring paths | fhfa_aersn_paths.pdf |
| Online Table S1: q=4 null frequencies | tableS1_size_q4.csv |
| Online Table S2: q=2 null frequencies | tableS2_size_q2.csv |
| Online Table S3: q=2 oscillation | tableS3_power_delay_q2.csv |
| Online Table S4: q=4 early alarms/delays | tableS4_early_delays_q4.csv |
| Online Table S5: pulse and reversal | tableS5_other_shapes_q4.csv |

The CSVs retain full precision. Main Table 1 rounds percentages to two decimals
and delays to one. PDF metadata and plot styling may differ between platforms;
the numeric series are the replication target.

## Designs, calibration and seeds

Simulations use N=40, m=120, T=2, one factor and Gaussian idiosyncratic errors.
Factor distributions are Gaussian and bounded (+/-1). Profiles 1--5 are flat,
front-cycle, back-cycle, rising and falling, normalized to training mean one.
The panel variance multiplier is v; the corresponding score variance is v^2.
The change starts after monitoring date 60. Shapes have equal squared energy,
and the loading-change scale is 0.12.

The q=2 study uses 8,000 Gaussian references and 3,000 panels for each independent
null-calibration/evaluation sample. The q=4 study uses 4,000 references and
2,000 panels per sample. Both methods share each panel. Critical values use
R quantile type 8 at 0.95. The known-profile reference uses the specified
full-horizon profile and the fitted-score correction. Size-adjusted power uses
an independent null sample from the same DGP.

For profile i and factor distribution d=1 (Gaussian), 2 (bounded):

| Stage | Seed |
| --- | --- |
| q=2 Gaussian reference | 83000+i |
| q=2 null calibration | 93000+i+100d |
| q=2 evaluation | 103000+i+100d |
| q=4 Gaussian reference | 183000+i |
| q=4 null calibration | 193000+i+100d |
| q=4 evaluation | 203000+i+100d |
| q=4 geometry validation | 213000+i |
| q=4 paired bootstrap (2,000 draws) | 223000+i+100d |
| FHFA calibration (2,000 draws) | 310001--312000 |
| FHFA independent null evaluation (1,000 draws) | 320001--321000 |

Detection is the probability of a first alarm after the change and by the
horizon. Early alarms are separate. Median delay conditions on detection.
Restricted mean delay conditions on no early alarm and assigns 181 to a
missed change; it therefore retains missed changes. Bootstrap intervals
resample calibration and paired evaluation panels separately and re-estimate
both critical values in each draw. The q=2 raw archives also preserve the
additional diagnostic shapes/calibration schemes; they are not added to the
manuscript tables.

## FHFA data and empirical reference

The input is the FHFA **All-Transactions metropolitan-area HPI** CSV retained
with the original project, containing observations through 2025Q4.
Source: https://www.fhfa.gov/data/hpi/datasets
Use this frozen vintage for exact replication; the live file is revised.
The input checksum is recorded in MANIFEST.sha256.

The code retains 405 complete metropolitan areas from 1994Q4 through 2025Q4,
computes quarterly log growth and standardizes using 1995Q1--2004Q4 moments.
Monitoring covers 2005Q1--2025Q4; m=40, K=3, J=1, q=6, gamma=0.

The empirical boundary is a training-fitted model reference: factor VAR(1)
innovations and whole cross-sectional idiosyncratic vectors are resampled,
with the latter independent across dates. Every replicated panel is
re-standardized and PCA is re-estimated on its own training window.
A 1e-8 ridge is used only for the small VAR coefficient regression;
**no ridge is used in the reported RSMS or SSMS normalizers**.
The independent null rates are 5.8% and 5.5%; observed alarms are 2007Q3 and
2007Q4 for RSMS and SSMS respectively.

## Software and numerical checks

Recorded environment: R 4.5.1 on macOS arm64, aersn 0.2.3, geometry 0.5.2,
data.table 1.17.8, sandwich 3.1-1 and lpSolve 5.6.23. aersn 0.2.3 is pinned
by the installer. The other packages are standard CRAN packages; versions
used on each run are recorded in results/sessionInfo.txt. No proprietary
software or restricted-access data are required.

The accelerated polygons/facets are checked against the exported
aersn_hull/aersn_gauge linear-program implementation. A common invertible
QR transformation conditions the gauge calculation without changing it;
this is not the LDL construction. Numerical checks also verify affine
invariance, bootstrap intervals, the source data panel and FHFA alarm dates.
See VALIDATION.md for checks actually performed for this release.

The reference/ files remain unchanged by both execution modes.
No journal correspondence or referee reports are included.

## Method reference

Hong, Y., Z. Lin, O. B. Linton, W. K. Newey and J. Sun (2026),
*Affine-equivariant adjusted-range self-normalization*, Cambridge Working
Papers in Economics 2678; Janeway Institute Working Paper 2637.
https://www.janeway.econ.cam.ac.uk/publication/affine-equivariant-adjusted-range-self-normalization

Package: https://cran.r-project.org/package=aersn

## Theory experiments in Table 2

All three studies use one factor, one Fourier pair, T=2, gamma=0 and a
5% KS critical value. An independent Gaussian calibration sample sets the
reference threshold in each study. Each finite panel re-estimates its
training PCA. The experiments answer three different probability questions.

- In panel A all monitoring observations are from the baseline regime.
  The Gaussian N=40 loading design changes on evenly spaced training dates.
  The training lengths are 100, 400 and 1600; the dispersed and concentrated
  alternatives have identical projected first-order training dose. A separate
  m=1600 experiment uses the front-cycle profile for both panel generation
  and Gaussian reference calibration. The CSV includes 999-draw paired
  bootstrap intervals that resample calibration and evaluation separately.
- Panel B uses the N=12 Gaussian loading design, m=80,240,720 and a unit-RMS
  oscillation after monitoring date m/2. The local amplitudes are h/sqrt(m)
  for h=0.7,1.6 and the fixed amplitude is 0.13. The CSV distinguishes any
  alarm, early alarms, first alarm after the break and detection among
  survivors. Both A and B use 4000 calibration and 2000 independent
  evaluation draws for panels and Gaussian references.
- Panel C uses the two-series Gaussian model defined in the manuscript.
  For each m=128,512,2048,8192 the same draws give a jump and a linear onset
  with amplitude 2. Delay is scaled by m*(sqrt(m)*2)^(-1/(nu+1)), with nu
  equal to zero or one. Early alarms are removed before summarizing delay,
  and survivors with a miss contribute positive infinity. The reference
  uses 6000 critical-value draws and 4000 separate exit draws on grid 4096.
  The order-statistic intervals for medians condition on the critical values.

Seed offsets for A are 840000 (calibration), 850000 (Gaussian evaluation),
860000 (panel evaluation) and 870000 (bootstrap), plus m and a further 50000
for the front-cycle profile. For B the three seeds are 453000+m, 463000+m
and 473000+m. For C each draw i uses 640000+i (calibration), 650000+i (exit)
or 660000+10000*j+i (panel, with j indexing the four sample sizes).
The scripts and the full-precision result CSVs provide all reported designs.

Quick mode rebuilds these tables from frozen raw draws. Full mode regenerates
the experiments from their seeds. Source scripts stay in the repository root.

compare_hac.R replays both four-coordinate profiles using the original seeds, two factor distributions and all three loading alternatives. The Bartlett bandwidth is floor(m^(1/3)); each method is separately calibrated. Use --cached after copying the reference draws to results/ to rebuild the tables and monochrome cumulative detection figure.

serial_panels.R generates the 15 Gaussian profile/dependence cells. serial_summary.R rebuilds power, delay and paired intervals from archived maxima and stopping times. figure_bandwidth.R plots the front-cycle oscillating design at rho=0.5. Bandwidths are 2, 4, 8 and 16. The calibration and evaluation samples are independent (2,000 and 1,500 panels per cell).
