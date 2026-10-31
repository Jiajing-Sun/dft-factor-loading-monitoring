# q=2 variance-profile experiment (Online Appendix Tables S2 and S3).  Reference, null
# calibration and evaluation use independent seed blocks.  Additional shapes/reference
# schemes are retained as raw diagnostics; the manuscript uses oscillation and the
# matched/size-adjusted calibrations.

source("setup.R")
source(file.path(ROOT, "monitoring_core.R"))
setwd(OUTPUT)
id <- as.integer(commandArgs(TRUE)[1])
m <- 120
N <- 40
T <- 2
n <- 360
k0 <- 60
M <- 180
B <- 3000
C <- 8000
rr <- (1:n)/m
v <- switch(id, rep(1, n), 1 + 0.65 * sin(2 * pi * rr), 1 - 0.65 * sin(2 * pi * rr), exp(0.35 *
    (rr - 0.5)), exp(-0.35 * (rr - 0.5)))
v <- v/mean(v[1:m])
H <- cbind(cos(2 * pi * rr), sin(2 * pi * rr))
profname <- c("flat", "front_cycle", "back_cycle", "rising", "falling")[id]
# Baseline loading and orthogonal loading-change direction, RMS normalized.
L <- sin(2 * pi * (1:N)/N) + 0.5
U <- L/sqrt(sum(L^2))
a <- rep(1, N) - U * sum(U)
a <- a/sqrt(mean(a^2))
cc <- cos(4 * pi * (1:N)/N)
Z <- cbind(L, rep(1, N))
cc <- cc - Z %*% solve(crossprod(Z), crossprod(Z, cc))
cc <- as.numeric(cc)/sqrt(mean(cc^2))
d <- 1:M
# Shapes have equal total squared loading-change energy.
rawg <- cbind(step = rep(1, M), ramp = d/M, pulse = as.numeric(d <= 45), reversal = ifelse(d <=
    90, 1, -1), oscillation = cos(2 * pi * (m + k0 + d)/m), cancellation = rep(1, M))
g <- sweep(rawg, 2, sqrt(colMeans(rawg^2)), "/")
write.csv(data.frame(date = 1:n, panel_variance = v, score_variance = v^2), paste0("profile_grid_",
    id, ".csv"), row.names = FALSE)
write.csv(data.frame(shape = colnames(g), energy = colSums(g^2), peak = apply(abs(g), 2,
    max), delta = 0.12), paste0("profile_signals_", id, ".csv"), row.names = FALSE)
# Gaussian reference for fitted scores; calibration quantiles use R type 8.
set.seed(83000 + id)
cal <- array(NA_real_, c(C, 2, 3), dimnames = list(NULL, c("AERSN", "SSMS"), c("flat", "known_profile",
    "wrong_variance_clock")))
for (b in 1:C) {
    e <- rnorm(n)
    z <- v * e
    cal[b, , 1] <- apply(seqstats(H * (e - mean(e[1:m])), m), 2, max)
    cal[b, , 2] <- apply(seqstats(H * (z - v * sum(z[1:m])/sum(v[1:m])), m), 2, max)
    cal[b, , 3] <- apply(seqstats(H * (z - v^2 * sum(z[1:m])/sum(v[1:m]^2)), m), 2, max)
}
crit <- apply(cal, c(2, 3), quantile, 0.95, type = 8)
saveRDS(cal, paste0("profile_reference_", id, ".rds"))
cat("reference complete", id, "\n")
flush.console()
# Generate factors first, then errors; this order is part of the fixed seeds.
panel <- function(dgp) {
    f <- if (dgp == "Gaussian")
        rnorm(n) else sample(c(-1, 1), n, TRUE)
    list(f = f, X = (tcrossprod(f, L) + matrix(rnorm(n * N), n, N)) * sqrt(v))
}
# PCA uses training observations only and is held fixed during monitoring.
fit <- function(X) {
    ee <- eigen(crossprod(X[1:m, ])/m, symmetric = TRUE)
    list(U = ee$vectors[, 1], ev = ee$values[1])
}
# Cross-sectional residual projection multiplied by the estimated factor, then by
# real/imaginary Fourier weights.
scores <- function(X, pc) {
    xp <- as.numeric(X %*% pc$U)
    H * (xp/sqrt(pc$ev) * (rowMeans(X) - xp * mean(pc$U)))
}
# Stopping times are Inf for misses, never replaced before computing power.
results <- list()
pairrows <- list()
rawall <- list()
thresholds <- list()
for (dgp in c("Gaussian", "Bounded")) {
    set.seed(93000 + id + 100 * match(dgp, c("Gaussian", "Bounded")))
    nullcal <- t(replicate(B, {
        X <- panel(dgp)$X
        apply(seqstats(scores(X, fit(X)), m), 2, max)
    }))
    cv <- cbind(crit, size_adjusted = apply(nullcal, 2, quantile, 0.95, type = 8))
    thresholds[[dgp]] <- cv
    saveRDS(nullcal, paste0("profile_nullcal_", id, "_", dgp, ".rds"))
    scenarios <- c("null", colnames(g))
    raw <- array(Inf, c(B, 7, 2, 4), dimnames = list(NULL, scenarios, rownames(cv), colnames(cv)))
    set.seed(103000 + id + 100 * match(dgp, c("Gaussian", "Bounded")))
    for (b in 1:B) {
        dat <- panel(dgp)
        pc <- fit(dat$X)
        for (j in 1:7) {
            X <- dat$X
            if (j > 1) {
                tt <- (m + k0 + 1):n
                avec <- if (scenarios[j] == "cancellation")
                  cc else a
                X[tt, ] <- X[tt, ] + 0.12 * tcrossprod(sqrt(v[tt]) * dat$f[tt] * g[, j -
                  1], avec)
            }
            ss <- seqstats(scores(X, pc), m)
            for (sch in 1:4) for (met in 1:2) {
                hits <- which(ss[, met] > cv[met, sch])
                if (length(hits))
                  raw[b, j, met, sch] <- hits[1]
            }
        }
    }
    rawall[[dgp]] <- raw
    # Delay conditions on no early alarm, assigning horizon+1 to a miss.
    for (j in 1:7) for (met in 1:2) for (sch in 1:4) {
        st <- raw[, j, met, sch]
        kk <- if (j == 1)
            0 else k0
        early <- st <= kk
        post <- st > kk & is.finite(st)
        delay <- st[post] - kk
        ci <- binom.test(sum(post), B)$conf.int
        results[[length(results) + 1]] <- data.frame(profile = profname, id = id, dgp = dgp,
            scenario = scenarios[j], method = rownames(cv)[met], scheme = colnames(cv)[sch],
            B = B, early = mean(early), detection = mean(post), lo = ci[1], hi = ci[2], miss = mean(!is.finite(st)),
            median_delay = if (length(delay))
                median(delay) else NA, p90_delay = if (length(delay))
                quantile(delay, 0.9) else NA, restricted_delay = mean(pmin(st[!early] - kk, 240 - kk + 1)), detect40 = mean(post &
                st - kk <= 40))
    }
    for (j in 2:7) for (sch in 1:4) {
        st <- raw[, j, , sch]
        diff <- as.numeric(st[, 1] > k0 & is.finite(st[, 1])) - as.numeric(st[, 2] > k0 &
            is.finite(st[, 2]))
        se <- sd(diff)/sqrt(B)
        pairrows[[length(pairrows) + 1]] <- data.frame(profile = profname, dgp = dgp, scenario = scenarios[j],
            scheme = colnames(cv)[sch], difference = mean(diff), lo = mean(diff) - 1.96 *
                se, hi = mean(diff) + 1.96 * se)
    }
    write.csv(do.call(rbind, results), paste0("profile_results_", id, ".csv"), row.names = FALSE)
    write.csv(do.call(rbind, pairrows), paste0("profile_pairs_", id, ".csv"), row.names = FALSE)
    saveRDS(list(profile = v, g = g, critical = thresholds, raw = rawall), paste0("profile_raw_",
        id, ".rds"))
    cat("evaluation complete", id, dgp, "\n")
    flush.console()
}
