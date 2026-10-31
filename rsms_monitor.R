# RSMS has affine-equivariant and LDL-based constructions.  The method names are
# RSMS-AE and RSMS-LDL; aersn is the implementation package.  All calibration must
# match this sequential statistic and its fitted Fourier scores.

rsms_monitor <- function(psi, m, construction = c("affine", "ldl"), gamma = 0) {
    construction <- match.arg(construction)
    stopifnot(is.matrix(psi), nrow(psi) > m, gamma >= 0, gamma < 0.5)
    # Construct the centered, scaled training partial-sum path.
    training <- psi[1:m, , drop = FALSE]
    path <- aersn::aersn_path(training)
    mon <- sweep(psi[(m + 1):nrow(psi), , drop = FALSE], 2, colMeans(training))
    Z <- apply(mon, 2, cumsum)/sqrt(m)
    if (ncol(psi) == 1)
        Z <- matrix(Z, ncol = 1)
    s <- seq_len(nrow(Z))/m
    Z <- Z/((1 + s) * (s/(1 + s))^gamma)
    if (construction == "affine") {
        h <- aersn::aersn_hull(path, n_directions = 0)
        stat <- aersn::aersn_gauge(h, Z)^2
    } else {
        # Appendix construction: lag-zero LDL, then coordinatewise path ranges.  This
        # depends on coordinate ordering; it is separate from the affine hull.
        Sigma0 <- crossprod(sweep(training, 2, colMeans(training)))/m
        fac <- aersn::aersn_ldl(Sigma0)
        Linv <- forwardsolve(fac$L, diag(ncol(psi)))
        Q <- path$G %*% t(Linv)
        ranges <- apply(Q, 2, function(x) diff(range(x)))
        stopifnot(all(ranges > 0))
        stat <- rowSums(sweep(Z %*% t(Linv), 2, ranges, "/")^2)
    }
    list(path = stat, KS = max(stat), construction = construction)
}
