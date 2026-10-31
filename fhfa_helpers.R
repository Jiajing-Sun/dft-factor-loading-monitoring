# FHFA panel preparation and factor VAR(1) reference model.
fhfa_yq_index <- function(year, quarter) {
    as.integer(year) * 4L + as.integer(quarter)
}

fhfa_yq_label <- function(year, quarter) {
    sprintf("%dQ%d", as.integer(year), as.integer(quarter))
}

# The 1e-8 ridge is used only in fitting the VAR coefficients, not in the RSMS or SSMS
# normalizers.
safe_solve <- function(A, ridge = 1e-10) {
    # Numerically stable inverse for small positive-definite matrices.
    stopifnot(is.matrix(A), nrow(A) == ncol(A))
    q <- nrow(A)
    A2 <- A + ridge * diag(q)

    out <- tryCatch({
        R <- chol(A2)
        chol2inv(R)
    }, error = function(e) {
        solve(A2)
    })
    return(out)
}

fit_var1_bootstrap <- function(F_mat, ridge = 1e-08) {
    F_mat <- as.matrix(F_mat)
    TT <- nrow(F_mat)
    K <- ncol(F_mat)
    if (TT < K + 3L)
        stop("Too few factor observations to fit VAR(1).")

    Y <- F_mat[2:TT, , drop = FALSE]
    X <- cbind(Intercept = 1, F_mat[1:(TT - 1L), , drop = FALSE])
    XtX <- crossprod(X)
    coef <- safe_solve(XtX, ridge = ridge) %*% crossprod(X, Y)
    fitted <- X %*% coef
    innovations <- Y - fitted

    list(intercept = as.numeric(coef[1, ]), A = t(coef[-1, , drop = FALSE]), innovations = innovations,
        F_initial_pool = F_mat, K = K)
}

simulate_var1_bootstrap <- function(var_fit, n, seed = NULL) {
    if (!is.null(seed))
        set.seed(seed)
    n <- as.integer(n)
    K <- as.integer(var_fit$K)
    F_sim <- matrix(0, nrow = n, ncol = K)

    init_id <- sample(seq_len(nrow(var_fit$F_initial_pool)), size = 1L)
    F_sim[1, ] <- var_fit$F_initial_pool[init_id, ]

    innov <- as.matrix(var_fit$innovations)
    for (tt in 2:n) {
        e_id <- sample(seq_len(nrow(innov)), size = 1L)
        F_sim[tt, ] <- var_fit$intercept + as.numeric(var_fit$A %*% F_sim[tt - 1L, ]) + innov[e_id,
            ]
    }

    colnames(F_sim) <- paste0("F", seq_len(K))
    F_sim
}

read_fhfa_growth_panel <- function(hpi_file = file.path("data", "hpi_at_metro.csv"), level_start = c(1994,
    4), level_end = c(2025, 4), train_start = c(1995, 1), train_end = c(2004, 4)) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required. Install via install.packages('data.table').")
    }
    if (!file.exists(hpi_file))
        stop("Cannot find FHFA file: ", hpi_file)

    dt <- data.table::fread(hpi_file, header = FALSE, col.names = c("metro", "cbsa", "year",
        "quarter", "hpi_raw", "se_raw"), na.strings = c("-", ""))

    suppressWarnings(dt[, hpi := as.numeric(hpi_raw)])
    dt[, t_idx := fhfa_yq_index(year, quarter)]
    dt[, yq := fhfa_yq_label(year, quarter)]

    start_idx <- fhfa_yq_index(level_start[1], level_start[2])
    end_idx <- fhfa_yq_index(level_end[1], level_end[2])
    dt <- dt[t_idx >= start_idx & t_idx <= end_idx]

    n_quarters_levels <- as.integer(end_idx - start_idx + 1L)
    keep_cbsa <- dt[, .(n_nonmiss = sum(!is.na(hpi))), by = cbsa][n_nonmiss == n_quarters_levels,
        cbsa]
    if (length(keep_cbsa) == 0L)
        stop("No complete FHFA metro series in the requested level window.")

    dt_bal <- dt[cbsa %in% keep_cbsa]
    data.table::setkey(dt_bal, t_idx, cbsa)
    wide <- data.table::dcast(dt_bal, t_idx + yq ~ cbsa, value.var = "hpi")
    data.table::setorder(wide, t_idx)

    H <- as.matrix(wide[, -(1:2), with = FALSE])
    col_cbsa <- colnames(H)
    logH <- log(H)
    X_growth <- 100 * (logH[-1, , drop = FALSE] - logH[-nrow(logH), , drop = FALSE])

    time_levels <- wide$yq
    time_growth <- time_levels[-1]

    train_start_label <- fhfa_yq_label(train_start[1], train_start[2])
    train_end_label <- fhfa_yq_label(train_end[1], train_end[2])
    train_start_pos <- match(train_start_label, time_growth)
    train_end_pos <- match(train_end_label, time_growth)
    if (is.na(train_start_pos))
        stop("Training start quarter not found: ", train_start_label)
    if (is.na(train_end_pos))
        stop("Training end quarter not found: ", train_end_label)
    if (train_end_pos < train_start_pos)
        stop("train_end precedes train_start.")

    train_idx <- seq.int(train_start_pos, train_end_pos)
    mu_train <- colMeans(X_growth[train_idx, , drop = FALSE])
    sd_train <- apply(X_growth[train_idx, , drop = FALSE], 2, stats::sd)
    sd_train[sd_train < 1e-08] <- 1

    X_std <- sweep(X_growth, 2, mu_train, FUN = "-")
    X_std <- sweep(X_std, 2, sd_train, FUN = "/")

    list(X_growth = X_growth, X_std = X_std, time_growth = time_growth, cbsa = col_cbsa,
        train_idx = train_idx, m_empirical = length(train_idx), mu_train = mu_train, sd_train = sd_train,
        level_start = level_start, level_end = level_end, train_start = train_start, train_end = train_end)
}
