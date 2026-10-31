# Shared paths and dependency checks. Run scripts from the repository root.
# DFT_RESULTS selects a different one-level output directory when needed.
ROOT <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(ROOT, "run_all.R"))) {
    stop("Run this script with the repository root as the working directory.")
}
OUTPUT <- Sys.getenv("DFT_RESULTS", "results")
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)
OUTPUT <- normalizePath(OUTPUT, mustWork = TRUE)
for (pkg in c("aersn", "geometry", "data.table")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package ", pkg, ". Run Rscript install_dependencies.R first.")
    }
}
if (as.character(utils::packageVersion("aersn")) != "0.2.3") {
    stop("The reported calculations require aersn 0.2.3.")
}
# Fix all three RNG algorithms; individual studies reset their documented seeds.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
