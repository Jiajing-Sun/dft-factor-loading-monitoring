# Replay the q=4 study on its original PCA panel streams. All three methods
# use the same observations. HAC follows the original Bartlett implementation.
source("setup.R")
out <- OUTPUT
cached <- "--cached" %in% commandArgs(TRUE)
hac <- function(psi, m) {
  x <- sweep(psi,2,colMeans(psi[1:m,,drop=FALSE]))
  tr <- x[1:m,,drop=FALSE]
  h <- floor(m^(1/3))
  V <- crossprod(tr)/m
  for (lag in seq_len(h-1L)) {
    g <- crossprod(tr[(lag+1):m,],tr[1:(m-lag),])/m
    V <- V+(1-lag/h)*(g+t(g))
  }
  Z <- apply(x[-(1:m),],2,cumsum)/sqrt(m)/(1+seq_len(nrow(x)-m)/m)
  rowSums((Z %*% solve(V+diag(1e-10,ncol(V))))*Z)
}
run_profile <- function(id) {
m <- 120L
N <- 40L
T <- 2L
n <- 360L
k0 <- 60L
M <- 180L
B <- 2000L
C <- 4000L
delta <- 0.12
rr <- (1:n)/m
v <- if (id == 1L) rep(1, n) else 1 + 0.65 * sin(2 * pi * rr)
v <- v/mean(v[1:m])
profile_name <- c("flat", "front_cycle")[id]
H <- do.call(cbind, lapply(1:2, function(j) cbind(cos(2 * pi * j * rr), sin(2 * pi * j *
    rr))))
L <- sin(2 * pi * (1:N)/N) + 0.5
U <- L/sqrt(sum(L^2))
a <- 1 - U * sum(U)
a <- a/sqrt(mean(a^2))
d <- 1:M
rawg <- cbind(oscillation = cos(2 * pi * (m + k0 + d)/m), reversal = ifelse(d <= 90, 1, -1),
    pulse = as.numeric(d <= 45))
g <- sweep(rawg, 2, sqrt(colMeans(rawg^2)), "/")
stopifnot(max(abs(colSums(g^2) - M)) < 1e-10, abs(sum(L * a)) < 1e-10, abs(sum(a^2) - N) <
    1e-10)
seeds <- c(reference = 183000 + id, Gaussian_null = 193100 + id, Bounded_null = 193200 +
    id, Gaussian_eval = 203100 + id, Bounded_eval = 203200 + id, validation = 213000 + id)

panel <- function(dgp) {
    f <- if (dgp == "Gaussian")
        rnorm(n) else sample(c(-1, 1), n, TRUE)
    list(f = f, X = (tcrossprod(f, L) + matrix(rnorm(n * N), n, N)) * sqrt(v))
}
fit <- function(X) {
    ee <- eigen(crossprod(X[1:m, , drop = FALSE])/m, symmetric = TRUE)
    list(U = ee$vectors[, 1], ev = ee$values[1])
}
scores <- function(X, pc) {
    xp <- as.numeric(X %*% pc$U)
    H * (xp/sqrt(pc$ev) * (rowMeans(X) - xp * mean(pc$U)))
}


  set.seed(seeds["reference"])
  cal <- replicate(C, {z <- v*rnorm(n); max(hac(H*(z-v*sum(z[1:m])/sum(v[1:m])),m))})
  saved <- list()
  for (dgp in c("Gaussian","Bounded")) {
    set.seed(seeds[paste0(dgp,"_null")])
    nul <- replicate(B,{ X<-panel(dgp)$X; max(hac(scores(X,fit(X)),m)) })
    cuts <- c(quantile(cal,.95,type=8),quantile(nul,.95,type=8))
    raw <- array(Inf,c(B,4,2),dimnames=list(NULL,c("null",colnames(g)),c("known_profile","size_adjusted")))
    set.seed(seeds[paste0(dgp,"_eval")])
    for(i in seq_len(B)){
      dat <- panel(dgp); pc <- fit(dat$X)
      for(j in seq_len(4)){
        X <- dat$X
        if(j>1){tt <- (m+k0+1):n; X[tt,] <- X[tt,]+delta*tcrossprod(sqrt(v[tt])*dat$f[tt]*g[,j-1],a)}
        path <- hac(scores(X,pc),m)
        for(k in 1:2){hit<-which(path>cuts[k]);if(length(hit)) raw[i,j,k]<-hit[1]}
      }
    }
    saved[[dgp]] <- list(raw=raw,nullcal=nul,critical=cuts)
  }
  saveRDS(list(cal=cal,cells=saved),file.path(out,paste0("hac_q4_",id,".rds")))
}
if(!cached) for(id in 1:2) run_profile(id)
# Form summary rows from first stopping times; early alarms are never detections.
rows <- list()
pdf(file.path(out,"three_methods.pdf"),width=7,height=5.4)
par(mfrow=c(2,2),mar=c(3.5,3.8,2,0.6),mgp=c(2.1,.65,0))
for(id in 1:2) for(dgp in c("Gaussian","Bounded")){
  getdraw <- function(name){
    path <- file.path(out,name)
    readRDS(if(file.exists(path)) path else file.path(ROOT,"reference",name))
  }
  original <- getdraw(paste0("profile_confirm_raw_",id,".rds"))$raw[[dgp]]
  hr <- getdraw(paste0("hac_q4_",id,".rds"))$cells[[dgp]]$raw
  profile <- c("Flat","Front-cycle")[id]
  plot(0,type="n",xlim=c(0,180),ylim=c(0,100),xlab="Observations after break",ylab="Detected (%)",main=paste(profile,dgp,sep=", "),las=1)
  for(ix in 1:3){
    met <- c("RSMS","SSMS","HAC")[ix]
    raw <- if(ix==3) hr else original[,,ix,]
    for(sc in c("oscillation","pulse","reversal")){
      st <- raw[,sc,"size_adjusted"]
      rows[[length(rows)+1]] <- data.frame(profile,dgp,method=met,shape=sc,
        size=100*mean(is.finite(raw[,"null","known_profile"])),adjusted_size=100*mean(is.finite(raw[,"null","size_adjusted"])),
        power=100*mean(is.finite(st)&st>60),early=100*mean(st<=60),
        median=median(st[is.finite(st)&st>60]-60),restricted=mean(pmin(st[st>60]-60,181)))
      if(sc=="oscillation"){
        y <- vapply(0:180,function(k)100*mean(st>60 & st<=60+k),numeric(1))
        lines(0:180,y,type="s",lty=c(1,2,3)[ix],lwd=1.5)
        at <- c(42,91,142)+4*ix
        points(at,y[at+1],pch=c(16,1,2)[ix],cex=.65)
      }
    }
  }
  legend("topleft",c("RSMS","SSMS","HAC"),lty=1:3,pch=c(16,1,2),cex=.72,bty="n")
}
dev.off()
write.csv(do.call(rbind,rows),file.path(out,"three_methods.csv"),row.names=FALSE)

# Output the manuscript rows directly from the simulated stopping times.
tab <- do.call(rbind,rows)
osc <- subset(tab,shape=="oscillation")
endrow <- paste0(intToUtf8(92),intToUtf8(92))
main_rows <- sprintf("%s & %s & %s & %.2f & %.2f & %.2f & %.1f & %.1f %s",
  osc$profile,osc$dgp,osc$method,osc$size,osc$adjusted_size,osc$power,osc$median,osc$restricted,endrow)
writeLines(main_rows,file.path(out,"three_method_rows.tex"))
more <- subset(tab,shape!="oscillation")
appendix_rows <- character()
for(profile in c("Flat","Front-cycle")) for(dgp in c("Gaussian","Bounded")) for(shape in c("pulse","reversal")){
  d <- more[more$profile==profile & more$dgp==dgp & more$shape==shape,]
  stopifnot(identical(d$method,c("RSMS","SSMS","HAC")))
  appendix_rows <- c(appendix_rows,paste0(paste(c(profile,dgp,shape,
    sprintf("%.2f",d$power),sprintf("%.1f",d$restricted)),collapse=" & ")," ",endrow))
}
writeLines(appendix_rows,file.path(out,"three_shape_rows.tex"))
