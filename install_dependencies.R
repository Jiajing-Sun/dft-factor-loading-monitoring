# Install standard CRAN dependencies and pin the statistical engine to 0.2.3.
# R_LIBS_USER can select a user library; no machine-specific path is required.
options(repos = c(CRAN = "https://cloud.r-project.org"))
needed <- c("geometry", "data.table", "sandwich", "lpSolve")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
have_version <- if (requireNamespace("aersn", quietly = TRUE)) as.character(packageVersion("aersn")) else ""
if (have_version != "0.2.3") {
    archive <- "https://cran.r-project.org/src/contrib/Archive/aersn/aersn_0.2.3.tar.gz"
    current <- "https://cran.r-project.org/src/contrib/aersn_0.2.3.tar.gz"
    dest <- tempfile(fileext = ".tar.gz")
    ok <- tryCatch({
        download.file(current, dest, mode = "wb", quiet = TRUE)
        TRUE
    }, error = function(e) FALSE)
    if (!ok)
        download.file(archive, dest, mode = "wb")
    install.packages(dest, repos = NULL, type = "source")
    unlink(dest)
}
stopifnot(as.character(packageVersion("aersn")) == "0.2.3")
print(sessionInfo())
