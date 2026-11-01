# Section 4: fitted PCA panels, matched Gaussian local laws and consistency.
# Run from the replication root with monitoring_core.R alongside the scripts:
# Rscript simulate_local_power.R
# An integer argument runs one training size. --tables combines completed runs.
core <- Sys.getenv("DFT_CORE", "monitoring_core.R")
source(core)
out <- Sys.getenv("DFT_RESULTS", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
ms <- c(80L, 240L, 720L)
C <- 4000L
B <- 2000L
strengths <- c(null = 0, local_weak = 0.7, local_strong = 1.6)
methods <- c("RSMS", "SSMS")
fname <- function(m) file.path(out, sprintf("local_power_m%d.csv", m))

run_size <- function(m) {
    N <- 12L
    n <- 3L * m
    k0 <- m %/% 2L
    dates <- seq_len(n)
    times <- dates / m
    H <- cbind(cos(2*pi*times),sin(2*pi*times))
    post <- dates > m + k0
    # A unit-RMS oscillation starts immediately after the break.
    g <- numeric(n)
    g[post] <- cos(2*pi*times[post])
    g <- g / sqrt(mean(g[post]^2))
    L <- sin(2*pi*seq_len(N)/N) + 0.5
    U <- L / sqrt(sum(L^2))
    a <- 1 - U * sum(U)
    a <- a / sqrt(mean(a^2))
    b <- rep(1/N,N)
    # At the null the projected factor has variance one, independent of
    # the residual mean. The coefficient below is exact for this loading
    # perturbation because a is orthogonal to U and errors have covariance I.
    sigma_z <- sqrt(sum((b - U * sum(U*b))^2))
    lambda <- 1 + sum(L^2)
    effect <- sum(U * L) * sum(a * b) / (sqrt(lambda) * sigma_z)

    # Separate reference calibration, Gaussian local-law evaluation and
    # actual-panel evaluation. Both methods and all changes share each draw.
    set.seed(453000L + m)
    maxima <- replicate(C, {
        e <- rnorm(n)
        psi <- H * (e - mean(e[seq_len(m)]))
        apply(seqstats(psi,m),2,max)
    })
    cv <- apply(maxima,1,quantile,probs=.95,type=8,names=FALSE)
    first_stop <- function(psi) {
        paths <- seqstats(psi,m)
        vapply(1:2,function(j) {
            ix <- which(paths[,j] > cv[j])
            if (length(ix)) ix[1] else Inf
        },numeric(1))
    }
    # The fitted Gaussian null path uses the same training endpoint as the
    # weighted score integral. Under a local change its monitoring increment
    # acquires h H g * effect/sqrt(m); the training path stays identical.
    ref <- array(Inf,c(B,2,length(strengths)),
                 dimnames=list(NULL,methods,names(strengths)))
    set.seed(463000L + m)
    for (i in seq_len(B)) {
        e <- rnorm(n)
        e <- e - mean(e[seq_len(m)])
        for (d in names(strengths))
            ref[i,,d] <- first_stop(H * (e + strengths[d]*effect*g/sqrt(m)))
    }

    delta <- c(strengths/sqrt(m),fixed=0.13)
    raw <- array(Inf,c(B,2,length(delta)),
                 dimnames=list(NULL,methods,names(delta)))
    set.seed(473000L + m)
    for (i in seq_len(B)) {
        f <- rnorm(n)
        X <- tcrossprod(f,L) + matrix(rnorm(n*N),n,N)
        pc <- eigen(crossprod(X[seq_len(m),,drop=FALSE])/m,symmetric=TRUE)
        u <- pc$vectors[,1]
        # Apply the actual training PCA to each post-break perturbation.
        # X_delta = X + delta*(g*f) a'; projections below avoid storing
        # an extra n-by-N matrix for every strength.
        xp0 <- as.vector(X %*% u)
        xb0 <- rowMeans(X)
        for (d in names(delta)) {
            amp <- delta[d] * f * g
            xp <- xp0 + amp * sum(a*u)
            xb <- xb0 + amp * mean(a)
            z <- xp / sqrt(pc$values[1]) * (xb - xp * mean(u))
            raw[i,,d] <- first_stop(H * z)
        }
    }
    sim <- list(m=m,N=N,k0=k0,delta=delta,effect=effect,
                reference_maxima=maxima,critical=cv,panel_stops=raw,reference_stops=ref)
    saveRDS(sim,file.path(out,sprintf("local_power_paths_m%d.rds",m)))
    summarize_draws(sim)
    message("Completed m=",m)
}

# All estimates are rebuilt from the stored stopping times, including misses.
summarize_draws <- function(sim) {
    m <- sim$m
    N <- sim$N
    k0 <- sim$k0
    delta <- sim$delta
    effect <- sim$effect
    raw <- sim$panel_stops
    ref <- sim$reference_stops
    cv <- sim$critical
    B <- dim(raw)[1]
    C <- ncol(sim$reference_maxima)
    rows <- list()
    for (mode in c("fitted_panel","shifted_reference")) {
        draws <- if (mode == "fitted_panel") raw else ref
        for (d in dimnames(draws)[[3]]) for (j in 1:2) {
            st <- draws[,j,d]
            survive <- st > k0
            detected <- is.finite(st)
            rows[[length(rows)+1L]] <- data.frame(m=m,N=N,T=2,k0=k0,C=C,B=B,
                method=methods[j],model=mode,scenario=d,delta=unname(delta[d]),
                h=unname(sqrt(m)*delta[d]),standardized_effect=effect,
                critical=cv[j],any_alarm=mean(detected),early=mean(!survive),
                post_break=mean(survive & detected),conditional=mean(detected[survive]),
                at_risk=sum(survive))
        }
    }
    write.csv(do.call(rbind,rows),fname(m),row.names=FALSE)
    do.call(rbind,rows)
}

make_tables <- function() {
    x <- do.call(rbind,lapply(ms,function(m) {
        sim <- readRDS(file.path(out,sprintf("local_power_paths_m%d.rds",m)))
        summarize_draws(sim)
    }))
    write.csv(x,file.path(out,"local_power_results.csv"),row.names=FALSE)
    take <- function(m,j,mode,sc,v) {
        i <- x$m==m & x$method==j & x$model==mode & x$scenario==sc
        stopifnot(sum(i)==1)
        100*x[[v]][i]
    }
    lines <- c("\\begin{tabular}{rlrrrrrr}", "\\toprule",
        "$m$ & Method & Null & \\multicolumn{2}{c}{$h=0.7$} & \\multicolumn{2}{c}{$h=1.6$} & Fixed \\\\",
        " & & & Panel & Ref. & Panel & Ref. & Panel \\\\", "\\midrule")
    for(m in ms) for(j in methods) {
        p <- function(sc,v="post_break") take(m,j,"fitted_panel",sc,v)
        r <- function(sc) take(m,j,"shifted_reference",sc,"post_break")
        lines <- c(lines,sprintf("%d & %s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\",
           m,j,p("null","any_alarm"),p("local_weak"),r("local_weak"),
           p("local_strong"),r("local_strong"),p("fixed")))
    }
    writeLines(c(lines,"\\bottomrule","\\end{tabular}"),file.path(out,"power_table.tex"))
}
args <- commandArgs(TRUE)
if (!length(args)) {
    for (m in ms) run_size(m)
    make_tables()
} else if (identical(args[1],"--tables")) {
    make_tables()
} else {
    m <- as.integer(args[1])
    stopifnot(m %in% ms)
    run_size(m)
}
