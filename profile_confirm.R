# q=4 experiment (main Table 1; Online Appendix Tables S1, S4 and S5).  Both statistics
# use identical fitted PCA-DFT scores and paired panel draws.  Independent null
# simulations set size-adjusted thresholds before evaluation.
source("setup.R")
source(file.path(ROOT, "monitoring_core.R"))
setwd(OUTPUT)

args <- commandArgs(TRUE)
id <- as.integer(args[1])
stopifnot(id %in% 1:2)
m <- 120L
N <- 40L
T <- 2L
n <- 360L
k0 <- 60L
M <- 180L
B <- 2000L
C <- 4000L
delta <- 0.12
rr <- (1:n)/m
v <- if (id == 1L) rep(1, n) else 1 + 0.65 * sin(2 * pi * rr)
v <- v/mean(v[1:m])
profile_name <- c("flat", "front_cycle")[id]
H <- do.call(cbind, lapply(1:2, function(j) cbind(cos(2 * pi * j * rr), sin(2 * pi * j *
    rr))))
L <- sin(2 * pi * (1:N)/N) + 0.5
U <- L/sqrt(sum(L^2))
a <- 1 - U * sum(U)
a <- a/sqrt(mean(a^2))
d <- 1:M
rawg <- cbind(oscillation = cos(2 * pi * (m + k0 + d)/m), reversal = ifelse(d <= 90, 1, -1),
    pulse = as.numeric(d <= 45))
g <- sweep(rawg, 2, sqrt(colMeans(rawg^2)), "/")
stopifnot(max(abs(colSums(g^2) - M)) < 1e-10, abs(sum(L * a)) < 1e-10, abs(sum(a^2) - N) <
    1e-10)
seeds <- c(reference = 183000 + id, Gaussian_null = 193100 + id, Bounded_null = 193200 +
    id, Gaussian_eval = 203100 + id, Bounded_eval = 203200 + id, validation = 213000 + id)

panel <- function(dgp) {
    f <- if (dgp == "Gaussian")
        rnorm(n) else sample(c(-1, 1), n, TRUE)
    list(f = f, X = (tcrossprod(f, L) + matrix(rnorm(n * N), n, N)) * sqrt(v))
}
fit <- function(X) {
    ee <- eigen(crossprod(X[1:m, , drop = FALSE])/m, symmetric = TRUE)
    list(U = ee$vectors[, 1], ev = ee$values[1])
}
scores <- function(X, pc) {
    xp <- as.numeric(X %*% pc$U)
    H * (xp/sqrt(pc$ev) * (rowMeans(X) - xp * mean(pc$U)))
}

# Independent validation paths: seqstats(check=TRUE) checks 12 dates against the
# exported aersn_hull/aersn_gauge LP implementation, tolerance 1e-6.
set.seed(seeds["validation"])
validation_count <- 0L
for (j in 1:3) {
    e <- rnorm(n)
    z <- v * e
    invisible(seqstats(H * (z - v * sum(z[1:m])/sum(v[1:m])), m, check = TRUE))
    validation_count <- validation_count + 1L
}
for (dgp in c("Gaussian", "Bounded")) for (j in 1:3) {
    dat <- panel(dgp)
    pc <- fit(dat$X)
    tt <- (m + k0 + 1):n
    dat$X[tt, ] <- dat$X[tt, ] + delta * tcrossprod(sqrt(v[tt]) * dat$f[tt] * g[, j], a)
    invisible(seqstats(scores(dat$X, pc), m, check = TRUE))
    validation_count <- validation_count + 1L
}
cat("LP validation complete:", id, validation_count, "paths\n")
flush.console()

set.seed(seeds["reference"])
cal <- matrix(NA_real_, C, 2, dimnames = list(NULL, c("AERSN", "SSMS")))
for (b in 1:C) {
    e <- rnorm(n)
    z <- v * e
    cal[b, ] <- apply(seqstats(H * (z - v * sum(z[1:m])/sum(v[1:m])), m), 2, max)
}
known <- apply(cal, 2, quantile, 0.95, type = 8)
saveRDS(cal, paste0("profile_confirm_reference_", id, ".rds"))
cat("reference complete:", id, "\n")
flush.console()

scenarios <- c("null", colnames(g))
result_rows <- list()
pair_rows <- list()
all_raw <- list()
all_maxima <- list()
all_critical <- list()
for (dgp in c("Gaussian", "Bounded")) {
    set.seed(seeds[paste0(dgp, "_null")])
    nullcal <- t(replicate(B, {
        X <- panel(dgp)$X
        apply(seqstats(scores(X, fit(X)), m), 2, max)
    }))
    cv <- cbind(known_profile = known, size_adjusted = apply(nullcal, 2, quantile, 0.95,
        type = 8))
    all_critical[[dgp]] <- cv
    saveRDS(nullcal, paste0("profile_confirm_nullcal_", id, "_", dgp, ".rds"))
    cat("null calibration complete:", id, dgp, "\n")
    flush.console()

    raw <- array(Inf, c(B, 4, 2, 2), dimnames = list(NULL, scenarios, rownames(cv), colnames(cv)))
    # These maxima allow threshold-resampling comparisons of early/post-change event
    # probabilities without rerunning PCA; delay resampling needs replay.
    maxima <- array(NA_real_, c(B, 4, 2, 2), dimnames = list(NULL, scenarios, rownames(cv),
        c("pre60", "post60")))
    set.seed(seeds[paste0(dgp, "_eval")])
    for (b in 1:B) {
        dat <- panel(dgp)
        pc <- fit(dat$X)
        for (j in 1:4) {
            X <- dat$X
            if (j > 1) {
                tt <- (m + k0 + 1):n
                X[tt, ] <- X[tt, ] + delta * tcrossprod(sqrt(v[tt]) * dat$f[tt] * g[, j -
                  1], a)
            }
            ss <- seqstats(scores(X, pc), m)
            maxima[b, j, , 1] <- apply(ss[1:k0, , drop = FALSE], 2, max)
            maxima[b, j, , 2] <- apply(ss[(k0 + 1):(m * T), , drop = FALSE], 2, max)
            for (sch in 1:2) for (met in 1:2) {
                hits <- which(ss[, met] > cv[met, sch])
                if (length(hits))
                  raw[b, j, met, sch] <- hits[1]
            }
        }
    }
    all_raw[[dgp]] <- raw
    all_maxima[[dgp]] <- maxima

    for (j in 1:4) for (met in 1:2) for (sch in 1:2) {
        st <- raw[, j, met, sch]
        kk <- if (j == 1)
            0L else k0
        early <- st <= kk
        post <- st > kk & is.finite(st)
        delay <- st[post] - kk
        ci <- binom.test(sum(post), B)$conf.int
        result_rows[[length(result_rows) + 1L]] <- data.frame(id = id, profile = profile_name,
            N = N, m = m, T = T, K = 1, J = 2, q = 4, dgp = dgp, scenario = scenarios[j],
            delta = if (j == 1)
                0 else delta, k0 = kk, method = rownames(cv)[met], scheme = colnames(cv)[sch],
            B = B, early = mean(early), detection = mean(post), lo = ci[1], hi = ci[2], miss = mean(!is.finite(st)),
            conditional_detection = mean(is.finite(st[!early])), median_delay = if (length(delay))
                median(delay) else NA_real_, p90_delay = if (length(delay))
                unname(quantile(delay, 0.9)) else NA_real_, restricted_delay = mean(pmin(st[!early] - kk, m * T - kk + 1)),
            detect40 = mean(post & st - kk <= 40))
    }
    for (j in 2:4) for (sch in 1:2) {
        st <- raw[, j, , sch]
        dd <- as.numeric(st[, 1] > k0 & is.finite(st[, 1])) - as.numeric(st[, 2] > k0 & is.finite(st[,
            2]))
        se <- sd(dd)/sqrt(B)
        pair_rows[[length(pair_rows) + 1L]] <- data.frame(id = id, profile = profile_name,
            dgp = dgp, scenario = scenarios[j], scheme = colnames(cv)[sch], difference = mean(dd),
            lo = mean(dd) - 1.96 * se, hi = mean(dd) + 1.96 * se, B = B)
    }
    write.csv(do.call(rbind, result_rows), paste0("profile_confirm_results_", id, ".csv"),
        row.names = FALSE)
    write.csv(do.call(rbind, pair_rows), paste0("profile_confirm_pairs_", id, ".csv"), row.names = FALSE)
    saveRDS(list(config = list(N = N, m = m, T = T, K = 1, J = 2, q = 4, B = B, C = C, delta = delta,
        k0 = k0, profile_name = profile_name, seeds = seeds, validation_paths = validation_count,
        LP_tolerance = 1e-06), panel_variance = v, score_variance = v^2, g = g, critical = all_critical,
        raw = all_raw, maxima = all_maxima), paste0("profile_confirm_raw_", id, ".rds"))
    cat("evaluation complete:", id, dgp, "\n")
    flush.console()
}
