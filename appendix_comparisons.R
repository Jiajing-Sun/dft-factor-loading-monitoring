# Rebuild the appendix from the independent evaluation samples used in the paper.
# Infinite stopping times are misses. A pre-break alarm is never a detection.
source("setup.R")
profiles <- c("Flat", "Front-cycle", "Back-cycle", "Rising", "Falling")
shapes <- c("Oscillation", "Pulse", "Reversal", "Step", "Ramp", "Cancellation")
methods <- c("RSMS", "SSMS", "HAC")
draws <- function(name) readRDS(file.path(OUTPUT, name))
rate <- function(event) 100 * mean(event)

# The conditional average describes successful detections. The restricted mean
# also retains survivors that miss the change, at one observation past the horizon.
timing <- function(stop) {
    hit <- stop > 60 & is.finite(stop)
    stopifnot(any(hit), all(stop[is.finite(stop)] >= 1 & stop[is.finite(stop)] <= 240))
    c(adjusted_power = rate(hit), early = rate(stop <= 60),
      ADD = mean(stop[hit] - 60), restricted = mean(pmin(stop[stop > 60] - 60, 181)))
}

independent <- function(q, id) {
    name <- if (q == 4) "profile_confirm_raw" else "profile_raw"
    sn <- draws(sprintf("%s_%d.rds", name, id))$raw
    hac <- draws(sprintf("hac_q%d_%d.rds", q, id))$cells
    rows <- list()
    for (dgp in names(sn)) for (shape in shapes[seq_len(if (q == 4) 3 else 6)]) {
        for (j in seq_along(methods)) {
            x <- if (j == 3) hac[[dgp]]$raw else sn[[dgp]][,,j,]
            s <- tolower(shape)
            fitted <- x[,s,"known_profile"]
            adj <- x[,s,"size_adjusted"]
            metrics <- c(size = rate(is.finite(x[,"null","known_profile"])),
                         adjusted_size = rate(is.finite(x[,"null","size_adjusted"])),
                         power = rate(fitted > 60 & is.finite(fitted)), timing(adj))
            rows[[length(rows)+1L]] <- data.frame(q, profile=profiles[id], factor=dgp,
                shape, method=methods[j], t(metrics), check.names=FALSE)
        }
    }
    do.call(rbind, rows)
}

serial <- function(id) {
    profile_id <- (id - 1L) %% 5L + 1L
    rho <- c(0, 0.5, 0.85)[(id - 1L) %/% 5L + 1L]
    x <- draws(sprintf("cell_%02d.rds", id))
    # Hold the profile-specific independence cutoff fixed to isolate the effect
    # of serial dependence. The adjusted cutoff matches both profile and rho.
    fixed <- draws(sprintf("cell_%02d.rds", profile_id))$cut
    methods <- c("RSMS", "SSMS", "HAC(2)", "HAC(4)", "HAC(8)", "HAC(16)")
    rows <- list()
    for (s in 1:3) for (j in seq_along(methods)) {
        null <- pmax(x$pre[,j,1], x$post[,j,1])
        before <- x$pre[,j,s+1L]
        after <- x$post[,j,s+1L]
        times <- timing(x$stop[,j,s+1L])
        stopifnot(isTRUE(all.equal(unname(times[1]), rate(before <= x$cut[j] & after > x$cut[j]))))
        metrics <- c(size = rate(null > fixed[j]), adjusted_size = rate(null > x$cut[j]),
                     power = rate(before <= fixed[j] & after > fixed[j]), times)
        rows[[length(rows)+1L]] <- data.frame(q=4, profile=profiles[profile_id],
            rho, shape=shapes[s], method=methods[j], t(metrics), check.names=FALSE)
    }
    do.call(rbind, rows)
}

# Each table repeats the calibration labels so it can be read on its own.
table_tex <- function(d, title, reference, by_factor) {
    end <- "\\\\"
    header <- if (by_factor) "Factor & " else ""
    cell <- function(x, previous) if (identical(x, previous)) "" else x
    rows <- character(nrow(d))
    for (i in seq_len(nrow(d))) {
        factor <- if (by_factor) c(cell(d$factor[i], d$factor[i-1L])) else NULL
        shape <- cell(d$shape[i], d$shape[i-1L])
        if (by_factor && i > 1 && d$factor[i] != d$factor[i-1L]) shape <- d$shape[i]
        values <- c(factor, shape, d$method[i], sprintf("%.2f", unlist(d[i,c("size",
            "adjusted_size","power","adjusted_power")])), sprintf("%.1f", unlist(d[i,c("ADD","restricted")])))
        rows[i] <- paste(paste(values, collapse=" & "), end)
    }
    c("\\begin{table}[!htbp]\\centering\\color{revisionred}",
      paste0("\\caption{",title,"}"), "\\begingroup\\fontsize{9.5}{11}\\selectfont\\setstretch{1.0}",
      "\\setlength{\\tabcolsep}{3pt}",
      paste0("\\begin{tabular}{",if(by_factor) "l" else "", "llrrrrrr}\\toprule"),
      paste0(header,"Change & Method & \\multicolumn{2}{c}{Size (\\%)} & \\multicolumn{2}{c}{Power (\\%)} & \\multicolumn{2}{c}{Adjusted delay} ",end),
      paste0(if(by_factor) "& " else "", "& & ",reference," & Adj. & ",reference," & Adj. & ADD & Restricted ",end,"\\midrule"),
      rows, "\\bottomrule\\end{tabular}\\endgroup\\end{table}")
}

pages <- character()
panels <- list()
for (q in c(4,2)) {
    for (id in seq_len(if (q == 4) 2 else 5)) {
        d <- independent(q,id)
        panels[[length(panels)+1L]] <- d
        pages <- c(pages, table_tex(d, sprintf("Independent panels, %s variance, $q=%d$.",
            tolower(profiles[id]),q),"Fitted",TRUE))
    }
}
grid <- lapply(1:15, serial)
for (d in grid) {
    pages <- c(pages, table_tex(d, sprintf("Gaussian AR(1) panels, %s variance, $\\rho=%s$, $q=4$.",
        tolower(d$profile[1]),d$rho[1]),"Indep.",FALSE))
}
write.csv(do.call(rbind,panels),file.path(OUTPUT,"appendix_independent.csv"),row.names=FALSE)
write.csv(do.call(rbind,grid),file.path(OUTPUT,"appendix_serial.csv"),row.names=FALSE)
writeLines(pages,file.path(OUTPUT,"appendix_comparisons.tex"))
