# Numerical checks of the released implementation and published results.  This is not
# an asymptotic size test: reported sizes use independent null panels.
source("setup.R")
source(file.path(ROOT, "monitoring_core.R"))
source(file.path(ROOT, "rsms_monitor.R"))
report <- character()
record <- function(x) {
    report <<- c(report, x)
    message(x)
}
set.seed(410031)
for (q in c(2L, 4L)) {
    psi <- matrix(rnorm(360L * q), 360L, q)
    accelerated <- seqstats(psi, 120L, check = TRUE)[, 1]
    public <- rsms_monitor(psi, 120L)$path
    stopifnot(max(abs(accelerated - public)) < 1e-06)
    # Translation and a nonsingular nonorthogonal mixing must preserve RSMS-AE.
    A <- diag(seq_len(q)) + matrix(0.15, q, q)
    shifted <- sweep(psi %*% A, 2, seq_len(q))
    transformed <- rsms_monitor(shifted, 120L)$path
    stopifnot(max(abs(public - transformed)) < 1e-06)
}
record("PASS: q=2/4 accelerated gauges equal public LP; affine invariance.")

# Compare regenerated table values against the frozen replication outputs.
for (id in 1:2) {
    new <- readRDS(file.path(OUTPUT, paste0("profile_confirm_raw_", id, ".rds")))
    ref <- readRDS(file.path(ROOT, "reference", paste0("profile_confirm_raw_", id, ".rds")))
    stopifnot(isTRUE(all.equal(new$raw, ref$raw, tolerance = 1e-08)))
    stopifnot(isTRUE(all.equal(new$critical, ref$critical, tolerance = 1e-08)))
}
for (id in 1:5) {
    new <- readRDS(file.path(OUTPUT, paste0("profile_raw_", id, ".rds")))
    ref <- readRDS(file.path(ROOT, "reference", paste0("profile_raw_", id, ".rds")))
    stopifnot(isTRUE(all.equal(new$raw, ref$raw, tolerance = 1e-08)))
}
record("PASS: stopping times and q=4 thresholds agree with frozen results.")
newboot <- read.csv(file.path(OUTPUT, "profile_confirm_bootstrap_pairs.csv"))
refboot <- read.csv(file.path(ROOT, "reference", "profile_confirm_bootstrap_pairs.csv"))
stopifnot(isTRUE(all.equal(newboot, refboot, tolerance = 1e-10)))
record("PASS: regenerated two-sample paired bootstrap reproduces reported intervals.")
new <- readRDS(file.path(OUTPUT, "fhfa_aersn_revision.rds"))
ref <- readRDS(file.path(ROOT, "reference", "fhfa_aersn_revision.rds"))
stopifnot(isTRUE(all.equal(new$panel$X_std, ref$panel$X_std, tolerance = 1e-12)), isTRUE(all.equal(new$actual,
    ref$actual, tolerance = 1e-08)), isTRUE(all.equal(new$critical, ref$critical, tolerance = 1e-08)))
tab <- read.csv(file.path(OUTPUT, "fhfa_aersn_results.csv"))
stopifnot(identical(tab$alarm_date, c("2007Q3", "2007Q4")), isTRUE(all.equal(tab$null_rate,
    c(0.058, 0.055))))
record("PASS: FHFA source panel, observed paths, boundaries and alarm dates.")
# Every table is emitted directly from stopping times; main Table 1 retains misses.
t1 <- read.csv(file.path(OUTPUT, "table1_power_delay.csv"))
stopifnot(nrow(t1) == 4L, all(t1$RSMS > t1$SSMS), all(t1$Delay_RSMS < t1$Delay_SSMS))
record("PASS: main Table 1 contains paired power and restricted mean delay.")
writeLines(report, file.path(OUTPUT, "validation.txt"))
