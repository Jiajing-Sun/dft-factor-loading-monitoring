# Monitoring geometry used by both Monte Carlo studies (gamma = 0).  AERSN is a
# historical array key; it denotes affine-equivariant RSMS throughout.  The calculation
# uses the unregularized training path.  The increment hull is
# conv{G_i - G_j}. Constructing it from the vertices of conv{G_i} gives the same set.
# The 2D orientation is counterclockwise, so edge normals provide the gauge through
# normalized supporting lines.
increment_vertices_2d <- function(G) {
    idx <- grDevices::chull(G[, 1L], G[, 2L])
    V <- G[idx, , drop = FALSE]
    nv <- nrow(V)
    ii <- rep(seq_len(nv), times = nv)
    jj <- rep(seq_len(nv), each = nv)
    D <- unique(V[ii, , drop = FALSE] - V[jj, , drop = FALSE])
    idx <- grDevices::chull(D[, 1L], D[, 2L])
    out <- D[idx, , drop = FALSE]
    out[rev(seq_len(nrow(out))), , drop = FALSE]
}
# Center monitoring scores by the training mean. The same invertible QR transformation
# is applied to the numerator and hull, preserving its gauge.  For q=4 use Qhull facet
# equations; check=TRUE compares 12 dates against the public aersn linear-program gauge
# to tolerance 1e-6.
seqstats <- function(psi, m, check = FALSE) {
    p <- aersn::aersn_path(psi[1:m, , drop = FALSE])
    G <- p$G
    q <- ncol(G)
    S <- apply(sweep(psi[(m + 1):nrow(psi), , drop = FALSE], 2, colMeans(psi[1:m, , drop = FALSE])),
        2, cumsum)/sqrt(m)
    Z <- S/(1 + seq_len(nrow(S))/m)
    # Same invertible QR transformation applied to numerator and hull.
    qrG <- qr(G)
    R <- qr.R(qrG)
    GG <- G %*% solve(R)
    ZZ <- Z %*% solve(R)
    if (q == 2) {
        V <- increment_vertices_2d(GG)
        E <- V[c(2:nrow(V), 1), ] - V
        A <- cbind(E[, 2], -E[, 1])
        A <- A/rowSums(A * V)
    } else {
        idx <- unique(as.vector(geometry::convhulln(GG)))
        P <- GG[idx, , drop = FALSE]
        ij <- expand.grid(i = 1:nrow(P), j = 1:nrow(P))
        D <- unique(P[ij$i, ] - P[ij$j, ])
        hh <- geometry::convhulln(D, options = "Qt", output.options = "n")
        # normals contain q normal coordinates followed by signed offset
        A <- hh$normals[, 1:q, drop = FALSE]/(-hh$normals[, q + 1])
    }
    ag <- apply(ZZ %*% t(A), 1, max)^2
    sq <- rowSums((Z %*% solve(crossprod(G)/m)) * Z)
    if (check) {
        h <- aersn::aersn_hull(p, n_directions = 0)
        sel <- unique(round(seq(1, nrow(Z), length.out = 12)))
        stopifnot(max(abs(sqrt(ag[sel]) - aersn::aersn_gauge(h, Z[sel, , drop = FALSE]))) <
            1e-06)
    }
    cbind(AERSN = ag, SSMS = sq)
}
