source("setup.R")
setwd(OUTPUT)
# Summarize fresh q=4 confirmation and propagate estimated-threshold uncertainty.
# Resample independent actual-null calibration panels and paired evaluation panels
# independently; resample methods and scenarios together within each.
out <- list()
boot_rows <- list()
null_rows <- list()
pairout <- list()
Rboot <- 2000L
for (id in 1:2) {
    saved <- readRDS(paste0("profile_confirm_raw_", id, ".rds"))
    out[[id]] <- read.csv(paste0("profile_confirm_results_", id, ".csv"))
    pairout[[id]] <- read.csv(paste0("profile_confirm_pairs_", id, ".csv"))
    for (dgp in c("Gaussian", "Bounded")) {
        maxima <- saved$maxima[[dgp]]
        nullcal <- readRDS(paste0("profile_confirm_nullcal_", id, "_", dgp, ".rds"))
        B <- dim(maxima)[1]
        scenarios <- dimnames(maxima)[[2]]
        observed_cv <- saved$critical[[dgp]][, "size_adjusted"]
        calc <- function(cvs, ii) {
            nul <- sapply(1:2, function(j) mean(pmax(maxima[ii, 1, j, 1], maxima[ii, 1, j,
                2]) > cvs[j]))
            det <- vapply(2:4, function(s) {
                p <- sapply(1:2, function(j) mean(maxima[ii, s, j, 1] <= cvs[j] & maxima[ii,
                  s, j, 2] > cvs[j]))
                p[1] - p[2]
            }, numeric(1))
            c(nul, det)
        }
        actual <- calc(observed_cv, 1:B)
        set.seed(223000 + id + 100 * match(dgp, c("Gaussian", "Bounded")))
        bs <- matrix(NA_real_, Rboot, 5)
        for (r in 1:Rboot) {
            ci <- sample.int(nrow(nullcal), nrow(nullcal), TRUE)
            cvs <- apply(nullcal[ci, , drop = FALSE], 2, quantile, 0.95, type = 8)
            ii <- sample.int(B, B, TRUE)
            bs[r, ] <- calc(cvs, ii)
        }
        for (j in 1:2) null_rows[[length(null_rows) + 1L]] <- data.frame(id = id, profile = saved$config$profile_name,
            dgp = dgp, method = c("AERSN", "SSMS")[j], scheme = "size_adjusted", rejection = actual[j],
            bootstrap_lo = unname(quantile(bs[, j], 0.025)), bootstrap_hi = unname(quantile(bs[,
                j], 0.975)), Rboot = Rboot)
        for (j in 1:3) boot_rows[[length(boot_rows) + 1L]] <- data.frame(id = id, profile = saved$config$profile_name,
            dgp = dgp, scenario = scenarios[j + 1], scheme = "size_adjusted", difference = actual[j +
                2], bootstrap_lo = unname(quantile(bs[, j + 2], 0.025)), bootstrap_hi = unname(quantile(bs[,
                j + 2], 0.975)), oscillation_family_lo = if (scenarios[j + 1] == "oscillation")
                unname(quantile(bs[, j + 2], 0.05/(2 * 4))) else NA_real_, oscillation_family_hi = if (scenarios[j + 1] == "oscillation")
                unname(quantile(bs[, j + 2], 1 - 0.05/(2 * 4))) else NA_real_, Rboot = Rboot)
        saveRDS(bs, paste0("profile_confirm_bootstrap_", id, "_", dgp, ".rds"))
    }
}
results <- do.call(rbind, out)
pairs <- do.call(rbind, pairout)
boot <- do.call(rbind, boot_rows)
nullboot <- do.call(rbind, null_rows)
write.csv(results, "profile_confirm_results.csv", row.names = FALSE)
write.csv(pairs, "profile_confirm_pairs.csv", row.names = FALSE)
write.csv(boot, "profile_confirm_bootstrap_pairs.csv", row.names = FALSE)
write.csv(nullboot, "profile_confirm_bootstrap_null.csv", row.names = FALSE)
cat("Two-sample paired bootstrap complete.\n")
