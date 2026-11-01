# Recompute manuscript entries from replication-level stopping times.  A miss is Inf.
# Detection excludes alarms at/before k0=60.  Restricted mean delay is conditional on
# no early alarm; a miss receives 181.
source("setup.R")
setwd(OUTPUT)
profile_names <- c("Flat", "Front-cycle", "Back-cycle", "Rising", "Falling")
metrics <- function(st, k0 = 60L) {
    early <- st <= k0
    detected <- st > k0 & is.finite(st)
    c(early = 100 * mean(early), detection = 100 * mean(detected), median_delay = if (any(detected)) median(st[detected] -
        k0) else NA_real_, restricted_delay = mean(pmin(st[!early] - k0, 241L - k0)))
}
boot <- read.csv("profile_confirm_bootstrap_pairs.csv")
power_row <- function(id, dgp, scenario = "oscillation", q = 4L) {
    saved <- readRDS(paste0(if (q == 4L)
        "profile_confirm_raw_" else "profile_raw_", id, ".rds"))
    st <- saved$raw[[dgp]][, scenario, , "size_adjusted"]
    a <- metrics(st[, 1])
    b <- metrics(st[, 2])
    out <- data.frame(Profile = profile_names[id], Factor = dgp, RSMS = a["detection"], SSMS = b["detection"],
        Gain = a["detection"] - b["detection"])
    if (q == 4L) {
    ci <- boot[boot$id == id & boot$dgp == dgp & boot$scenario == scenario, , drop = FALSE]
        stopifnot(nrow(ci) == 1L)
        out$CI_low <- 100 * ci$bootstrap_lo
        out$CI_high <- 100 * ci$bootstrap_hi
    }
    out$Delay_RSMS <- a["restricted_delay"]
    out$Delay_SSMS <- b["restricted_delay"]
    rownames(out) <- NULL
    out
}
rows <- function(ids, fun) {
    do.call(rbind, lapply(ids, function(i) do.call(rbind, lapply(c("Gaussian", "Bounded"),
        function(d) fun(i, d)))))
}
table1 <- rows(1:2, function(i, d) power_row(i, d))
write.csv(table1, "table1_power_delay.csv", row.names = FALSE)
# LaTeX body with the exact rounding used in main Table 1.
lines <- apply(table1, 1, function(x) {
    nums <- as.numeric(x[3:9])
    sprintf("%s & %s & %.2f & %.2f & %.2f & [%.2f, %.2f] & %.1f & %.1f \\\\", x[1], x[2],
        nums[1], nums[2], nums[3], nums[4], nums[5], nums[6], nums[7])
})
writeLines(lines, "table1_rows.tex")
null_rows <- function(q, ids) rows(ids, function(i, d) {
    saved <- readRDS(paste0(if (q == 4L)
        "profile_confirm_raw_" else "profile_raw_", i, ".rds"))
    raw <- saved$raw[[d]][, "null", , ]
    data.frame(Profile = profile_names[i], Factor = d, Fitted_RSMS = 100 * mean(is.finite(raw[,
        1, "known_profile"])), Fitted_SSMS = 100 * mean(is.finite(raw[, 2, "known_profile"])),
        Adjusted_RSMS = 100 * mean(is.finite(raw[, 1, "size_adjusted"])), Adjusted_SSMS = 100 *
            mean(is.finite(raw[, 2, "size_adjusted"])))
})
write.csv(null_rows(4L, 1:2), "tableS1_size_q4.csv", row.names = FALSE)
write.csv(null_rows(2L, 1:5), "tableS2_size_q2.csv", row.names = FALSE)
write.csv(rows(1:5, function(i, d) power_row(i, d, q = 2L)), "tableS3_power_delay_q2.csv",
    row.names = FALSE)
delays <- rows(1:2, function(i, d) {
    st <- readRDS(paste0("profile_confirm_raw_", i, ".rds"))$raw[[d]][, "oscillation", ,
        "size_adjusted"]
    a <- metrics(st[, 1])
    b <- metrics(st[, 2])
    data.frame(Profile = profile_names[i], Factor = d, Early_RSMS = a["early"], Early_SSMS = b["early"],
        Median_RSMS = a["median_delay"], Median_SSMS = b["median_delay"], Restricted_RSMS = a["restricted_delay"],
        Restricted_SSMS = b["restricted_delay"])
})
write.csv(delays, "tableS4_early_delays_q4.csv", row.names = FALSE)
other <- rows(1:2, function(i, d) do.call(rbind, lapply(c("pulse", "reversal"), function(s) {
    x <- power_row(i, d, s)
    x$Shape <- s
    x
})))
write.csv(other, "tableS5_other_shapes_q4.csv", row.names = FALSE)
# Main Table 3 is the fitted-reference FHFA summary, read after its re-evaluation.
fhfa <- read.csv("fhfa_aersn_results.csv")
fhfa$method[fhfa$method == "AERSN"] <- "RSMS"
write.csv(fhfa, "table3_fhfa.csv", row.names = FALSE)
cat("Built main Tables 1 and 3 and online Tables S1-S5 from stored draws.\n")

# CvM uses the same first-alarm definitions and censoring horizon as main Table 1.
cvm_rows <- list()
for (id in 1:2) {
    draw <- readRDS(sprintf("cvm_q4_%d.rds", id))
    for (factor in names(draw$cells)) {
        cell <- draw$cells[[factor]]
        # Recheck the stored integral thresholds against their own null sample.
        q95 <- function(x) apply(x, 1, quantile, .95, type=8, names=FALSE)
        stopifnot(isTRUE(all.equal(cell$critical,
            cbind(fitted=q95(draw$reference), adjusted=q95(cell$terminal_null)))))
        for (method in dimnames(cell$stops)[[3]]) {
            st <- cell$stops[,,method,]
            time <- st[,"oscillation","adjusted"]
            hit <- is.finite(time) & time > 60
            cvm_rows[[length(cvm_rows)+1L]] <- data.frame(
                profile=profile_names[id],factor,method,
                fitted_size=100*mean(is.finite(st[,"null","fitted"])),
                adjusted_size=100*mean(is.finite(st[,"null","adjusted"])),
                adjusted_power=100*mean(hit),ADD=mean(time[hit]-60),
                restricted=mean(pmin(time[time>60]-60,181)))
        }
    }
}
cvm <- do.call(rbind,cvm_rows)
write.csv(cvm,"cvm_q4_summary.csv",row.names=FALSE)
cvm_lines <- with(cvm,sprintf("%s & %s & %s & %.2f & %.2f & %.2f & %.1f & %.1f %s",
    profile,factor,method,fitted_size,adjusted_size,adjusted_power,ADD,restricted,"\\\\"))
writeLines(c("\\begin{tabular}{lllrrrrr}\\toprule",
    "Profile & Factor & Method & \\multicolumn{2}{c}{Size (\\%)} & Adjusted & \\multicolumn{2}{c}{Delay} \\\\",
    " & & & Fitted & Adj. & power (\\%) & ADD & Restricted \\\\\\midrule",
    cvm_lines,"\\bottomrule\\end{tabular}"),"cvm_table.tex")
