# Sourced inside run_profile() to share every panel and score with Table 1.
# With unit weights the CvM integral is monotone, so the terminal integral
# calibrates the first-crossing rule for the whole monitoring horizon.
source(file.path(ROOT,"monitoring_core.R"), local=TRUE)
methods <- c("RSMS","SSMS","HAC")
pointwise <- function(psi) {
    paths <- cbind(seqstats(psi,m),hac(psi,m))
    colnames(paths) <- methods
    paths
}
limits <- function(paths) cbind(CvM=colSums(paths)/m, KS=apply(paths,2,max))
first <- function(x, cutoff) {
    hit <- which(x > cutoff)
    if (length(hit)) hit[1] else Inf
}
read_ks <- function(name) readRDS(file.path(OUTPUT,name))
sn <- read_ks(sprintf("profile_confirm_raw_%d.rds",id))
bartlett <- read_ks(sprintf("hac_q4_%d.rds",id))
match_ks <- function(actual,expected) {
    stopifnot(isTRUE(all.equal(unname(actual),unname(expected),tolerance=1e-8)))
}
set.seed(seeds["reference"])
ref <- replicate(C,{
    z <- v*rnorm(n)
    limits(pointwise(H*(z-v*sum(z[1:m])/sum(v[1:m]))))
})
match_ks(t(ref[,"KS",]),cbind(read_ks(sprintf("profile_confirm_reference_%d.rds",id)),bartlett$cal))
quantiles <- function(draws) apply(draws,1,quantile,.95,type=8,names=FALSE)
fitted <- quantiles(ref[,"CvM",])
cells <- list()
for (factor in c("Gaussian","Bounded")) {
    set.seed(seeds[paste0(factor,"_null")])
    calibration <- replicate(B,{
        X <- panel(factor)$X
        limits(pointwise(scores(X,fit(X))))
    })
    match_ks(t(calibration[,"KS",]),cbind(read_ks(sprintf("profile_confirm_nullcal_%d_%s.rds",id,factor)),
                                         bartlett$cells[[factor]]$nullcal))
    cuts <- cbind(fitted=fitted,adjusted=quantiles(calibration[,"CvM",]))
    ks_cut <- c(sn$critical[[factor]][,"size_adjusted"],bartlett$cells[[factor]]$critical[2])
    old <- sn$raw[[factor]]
    stops <- array(Inf,c(B,2,3,2),dimnames=list(NULL,c("null","oscillation"),methods,colnames(cuts)))
    set.seed(seeds[paste0(factor,"_eval")])
    for (i in seq_len(B)) {
        dat <- panel(factor)
        X <- dat$X
        pc <- fit(X)
        for (shape in c("null","oscillation")) {
            if (shape=="oscillation") {
                tt <- (m+k0+1):n
                X[tt,] <- X[tt,]+delta*tcrossprod(sqrt(v[tt])*dat$f[tt]*g[,"oscillation"],a)
            }
            paths <- pointwise(scores(X,pc))
            ks_time <- vapply(1:3,function(j) first(paths[,j],ks_cut[j]),numeric(1))
            expected <- c(old[i,shape,,"size_adjusted"],bartlett$cells[[factor]]$raw[i,shape,"size_adjusted"])
            match_ks(ks_time,expected)
            integral <- apply(paths,2,cumsum)/m
            for (j in 1:3) for (k in 1:2)
                stops[i,shape,j,k] <- first(integral[,j],cuts[j,k])
        }
    }
    cells[[factor]] <- list(terminal_null=calibration[,"CvM",],critical=cuts,stops=stops)
}
saveRDS(list(reference=ref[,"CvM",],cells=cells),file.path(OUTPUT,sprintf("cvm_q4_%d.rds",id)))
