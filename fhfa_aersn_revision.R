# FHFA application: training-fitted model reference, re-estimated per panel.  --cached
# reuses the supplied independent calibration/evaluation draws, but always rebuilds the
# data panel and observed monitoring paths from CSV.
source("setup.R")
source(file.path(ROOT, "fhfa_helpers.R"))
library(aersn)
panel <- read_fhfa_growth_panel(file.path(ROOT, "data/hpi_at_metro.csv"))
setwd(OUTPUT)
cached <- "--cached" %in% commandArgs(TRUE)
m <- 40L
K <- 3L
T <- 2.1
n <- 124L
fastpca <- function(X) {
    ss <- svd(X, nu = 0, nv = K)
    U <- ss$v[, 1:K, drop = FALSE]
    ev <- ss$d[1:K]^2/nrow(X)
    F <- sweep(X %*% U, 2, sqrt(ev), "/")
    L <- sweep(U, 2, sqrt(ev), "*")
    list(U = U, ev = ev, F = F, L = L, E = X - F %*% t(L))
}
pc <- fastpca(panel$X_std[1:m, ])
vf <- fit_var1_bootstrap(pc$F)
P <- exp((0+2i) * pi * (1:n)/(n - m))
P <- P - mean(exp((0+2i) * pi * (1:(n - m))/(n - m)))
score <- function(X) {
    pp <- fastpca(X[1:m, ])
    xp <- X %*% pp$U
    z <- sweep(xp, 2, sqrt(pp$ev), "/") * as.numeric(rowMeans(X) - xp %*% colMeans(pp$U))
    cbind(z * Re(P), z * Im(P))
}
pathstats <- function(psi) {
    ap <- aersn_path(psi[1:m, , drop = FALSE])
    h <- aersn_hull(ap, n_directions = 0)
    stopifnot(h$full_rank)
    mon <- sweep(psi[(m + 1):n, ], 2, colMeans(psi[1:m, ]))
    Z <- apply(mon, 2, cumsum)/sqrt(m)/(1 + (1:(n - m))/m)
    cbind(AERSN = aersn_gauge(h, Z)^2, SSMS = rowSums((Z %*% solve(crossprod(ap$G)/m)) *
        Z))
}
actual <- pathstats(score(panel$X_std))
one <- function(i, base) {
    set.seed(base + i)
    F <- simulate_var1_bootstrap(vf, n)
    X <- F %*% t(pc$L) + pc$E[sample.int(m, n, TRUE), ]
    X <- sweep(sweep(X, 2, colMeans(X[1:m, ])), 2, apply(X[1:m, ], 2, sd), "/")
    apply(pathstats(score(X)), 2, max)
}
runbatch <- function(B, base) {
    # Each replication has its own seed, so serial and parallel runs agree.
    cores <- as.integer(Sys.getenv("DFT_CORES", "2"))
    stopifnot(is.finite(cores), cores >= 1L)
    out <- if (.Platform$OS.type == "windows" || cores == 1L)
        lapply(1:B, one, base = base) else parallel::mclapply(1:B, one, base = base, mc.cores = cores, mc.preschedule = TRUE,
        mc.set.seed = FALSE)
    stopifnot(all(vapply(out, is.numeric, logical(1))))
    do.call(rbind, out)
}
cat("data", nrow(panel$X_std), ncol(panel$X_std), "actual maxima", apply(actual, 2, max),
    "\n")
flush.console()
cal <- if (cached) readRDS("fhfa_aersn_nullcal.rds") else runbatch(2000, 310000)
saveRDS(cal, "fhfa_aersn_nullcal.rds")
crit <- apply(cal, 2, quantile, 0.95, type = 8)
cat("calibration", crit, "\n")
flush.console()
eval <- if (cached) readRDS("fhfa_aersn_nulleval.rds") else runbatch(1000, 320000)
saveRDS(eval, "fhfa_aersn_nulleval.rds")
rates <- colMeans(sweep(eval, 2, crit, ">"))
rows <- lapply(1:2, function(j) {
    hits <- which(actual[, j] > crit[j])
    k <- if (length(hits))
        hits[1] else NA_integer_
    ci <- binom.test(sum(eval[, j] > crit[j]), nrow(eval))$conf.int
    data.frame(method = colnames(actual)[j], q = 6, K = 3, J = 1, m = m, N = ncol(panel$X_std),
        scheme = "training_fitted_model", critical = crit[j], calibration_B = 2000, evaluation_B = 1000,
        null_rate = rates[j], null_lo = ci[1], null_hi = ci[2], stop_k = k, alarm_date = if (is.na(k))
            "none" else panel$time_growth[m + k], observed_max = max(actual[, j]))
})
tab <- do.call(rbind, rows)
write.csv(tab, "fhfa_aersn_results.csv", row.names = FALSE)
write.csv(data.frame(date = panel$time_growth[(m + 1):n], actual), "fhfa_aersn_paths.csv",
    row.names = FALSE)
saveRDS(list(panel = panel, pc = pc, vf = vf, actual = actual, critical = crit), "fhfa_aersn_revision.rds")
pdf("fhfa_aersn_paths.pdf", width = 8, height = 3.8)
par(mar = c(4, 4, 1, 1))
matplot(1:(n - m), sweep(actual, 2, crit, "/"), type = "l", lty = 1:2, col = c("#8B0000",
    "#303030"), log = "y", xaxt = "n", xlab = "Quarter", ylab = "Statistic / model-reference boundary")
at <- seq(1, n - m, by = 12)
axis(1, at = at, labels = panel$time_growth[m + at], cex.axis = 0.8)
abline(h = 1, lty = 3)
legend("topleft", c("RSMS", "SSMS"), lty = 1:2, col = c("#8B0000", "#303030"), bty = "n")
dev.off()

# Figure 1: raw quarterly growth, before training-window standardization.
pdf("fhfa_J1_cross_sectional_avg_growth.pdf", width = 8, height = 3.8)
par(mar = c(4, 4, 1, 1))
plot(seq_len(n), rowMeans(panel$X_growth), type = "l", col = "black", xaxt = "n", xlab = "Quarter",
    ylab = "Average quarterly log growth (%)")
at <- seq(1, n, by = 16)
axis(1, at = at, labels = panel$time_growth[at], cex.axis = 0.8)
abline(h = 0, lty = 3, col = "grey60")
dev.off()
print(tab)
