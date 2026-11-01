# Affine-equivariant adjusted-range monitoring of fitted Fourier scores.

rsms_monitor <- function(psi, m, construction = "affine", gamma = 0) {
    construction <- match.arg(construction, "affine")
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
    h <- aersn::aersn_hull(path, n_directions = 0)
    stat <- aersn::aersn_gauge(h, Z)^2
    list(path = stat, KS = max(stat), construction = construction)
}
