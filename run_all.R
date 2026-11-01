# One entry point. --quick rebuilds all tables/figures from frozen draws.  --full
# regenerates every reference, null-calibration and evaluation sample.  Run from this
# repository's root: Rscript run_all.R --quick
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "--quick"
if (!mode %in% c("--quick", "--full")) stop("Use --quick or --full.")
source("setup.R")
Sys.setenv(DFT_RESULTS = OUTPUT)  # Pass the absolute output path to each R process.
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, args = character()) {
    message("Running ", script, " ", paste(args, collapse = " "))
    status <- system2(rscript, c(shQuote(file.path(ROOT, script)), args))
    if (status != 0)
        stop("Failed: ", script)
}
start <- proc.time()[["elapsed"]]
if (mode == "--quick") {
    inputs <- list.files(file.path(ROOT, "reference"), full.names = TRUE)
    stopifnot(all(file.copy(inputs, OUTPUT, overwrite = TRUE)))
} else {
    for (id in 1:5) run("profile_run.R", as.character(id))
    for (id in 1:2) run("profile_confirm.R", as.character(id))
}
run("profile_confirm_summarize.R")
run("fhfa_aersn_revision.R", if (mode == "--quick") "--cached" else character())
run("simulate_training.R", if (mode == "--quick") "--cached" else character())
run("simulate_local_power.R", if (mode == "--quick") "--tables" else character())
run("simulate_delay.R", mode)
run("theory_table.R")
run("compare_hac.R", if (mode == "--quick") "--cached" else character())
run("build_tables.R")
run("validate.R")
writeLines(c(paste("Mode:", mode), paste("Elapsed seconds:", proc.time()[["elapsed"]] - start),
    capture.output(sessionInfo())), file.path(OUTPUT, "sessionInfo.txt"))
message("Complete. Tables, figures and validation are in ", OUTPUT)
