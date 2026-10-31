# Release validation

Validated on 31 October 2026, R 4.5.1, macOS arm64.

- The complete quick workflow ran successfully: bootstrap intervals, both
  main tables, all five online tables, both FHFA figures and numerical checks.
  Recorded elapsed time was about 7 seconds in the local environment.
- Fresh full q=2 flat-profile experiment: all 8,000 Gaussian-reference draws,
  both 3,000-panel null samples and both 3,000-panel evaluation samples
  reproduce the archived calibration arrays and stopping times.
- Fresh full q=4 front-cycle experiment: all 4,000 reference draws, both
  2,000-panel null samples and both 2,000-panel evaluation samples reproduce
  the archived calibration arrays, stopping times and pre/post maxima.
- Both factor distributions and all scenario/reference combinations in these
  two selected experiments were retained in the comparisons.
- Two-sample paired bootstrap: all 2,000 draws for each q=4 profile/factor
  cell were regenerated; the reported pointwise and simultaneous intervals
  agree with the archived output.
- Accelerated q=2/q=4 gauges match the public aersn linear-program gauge
  within 1e-6. Affine invariance passes the same tolerance.
- RSMS-LDL agrees with the quadratic form from aersn_ldl_normalizer within
  1e-8.
- The FHFA CSV reconstructs the exact standardized panel. Recomputed observed
  paths, boundaries and alarm dates agree with the archived application.
- Eight newly generated FHFA panels (four from each independent seed block)
  reproduce the archived maxima within 1e-8.
- Every entry in the manuscript's expanded Table 1 matches table1_rows.tex.
- SHA-256 checksums cover the frozen input and reference files.

The release checks combine complete reruns of two representative simulation
configurations with all-table regeneration and empirical seed replay.
The remaining simulation profiles and the full 3,000-panel FHFA reference
are supplied as archived draws; run_all.R --full regenerates them.
Cross-platform numerical comparisons use tolerances; PDF metadata need not
be byte-identical.
