# Paired kernel comparison on the original Figure 2 panel streams. Every
# bandwidth, kernel and stopping rule has its own independent null cutoff.
method <- c("RSMS","SSMS",paste0("Bartlett",hs),paste0("Parzen",hs))
functional <- c("KS","CvM")
original <- readRDS(sprintf("cell_%02d.rds",id))
set.seed(813000+id)
cal <- replicate(C,draw(FALSE))
dimnames(cal) <- list(method,functional,NULL)
stopifnot(isTRUE(all.equal(unname(t(cal[1:6,"KS",])),original$cal,tolerance=1e-8)))
cut <- apply(cal,1:2,quantile,probs=.95,type=8)
stop_time <- before <- after <- array(Inf,c(B,length(method),2,2),
    dimnames=list(NULL,method,functional,c("null","oscillation")))
set.seed(913000+id)
for (i in seq_len(B)) {
    pair <- draw(TRUE,i==1)
    for (shape in names(pair)) {
        instantaneous <- pair[[shape]]
        for (rule in functional) {
            path <- if (rule=="KS") instantaneous else apply(instantaneous,2,cumsum)/m
            before[i,,rule,shape] <- apply(path[1:k0,,drop=FALSE],2,max)
            after[i,,rule,shape] <- apply(path[-(1:k0),,drop=FALSE],2,max)
            for (j in seq_along(method)) {
                hit <- which(path[,j] > cut[j,rule])
                if (length(hit)) stop_time[i,j,rule,shape] <- hit[1]
            }
        }
    }
}
stopifnot(identical(unname(stop_time[,1:6,"KS",]),original$stop[,1:6,1:2]))
saveRDS(list(cal=cal,cut=cut,stop=stop_time,pre=before,post=after),"kernel_front_rho05.rds")
