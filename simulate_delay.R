# Finite-panel stopping times and the Gaussian ray-exit distribution.
# Run from the directory in which delay_*.csv and delay_table.tex should be written.
# Calibration, local-ray simulation and finite-panel paths have separate seeds.
stopifnot(as.character(packageVersion("aersn")) == "0.2.3")
RNGkind("Mersenne-Twister","Inversion","Rejection")
B <- 2000L
C <- 6000L
grid <- 4096L
sizes <- c(128L,512L,2048L,8192L)
amplitude <- 2
workers <- if (.Platform$OS.type=="windows") 1L else as.integer(Sys.getenv("DFT_CORES","3"))

geometry <- function(z) {
  path <- aersn::aersn_path(z)
  G <- path$G
  # conv(G-G) is the difference body of the polygon conv(G).
  V <- G[chull(G),,drop=FALSE]
  pairs <- expand.grid(i=seq_len(nrow(V)),j=seq_len(nrow(V)))
  D <- V[pairs$i,,drop=FALSE] - V[pairs$j,,drop=FALSE]
  W <- D[chull(D),,drop=FALSE]
  edge <- W[c(2:nrow(W),1),,drop=FALSE]-W
  normals <- cbind(edge[,2],-edge[,1])
  list(normals=normals/rowSums(normals*W), inv=solve(crossprod(G)/nrow(z)),
       mean=colMeans(z), path=path)
}
# Evaluate both KS paths on the same centered partial sums.
mon_path <- function(z,geo,m) {
  y <- cbind(cumsum(z[,1]-geo$mean[1]),cumsum(z[,2]-geo$mean[2]))/sqrt(m)
  g <- rep(0,nrow(y))
  for(j in seq_len(nrow(geo$normals))) {
    a <- geo$normals[j,]
    g <- pmax(g,a[1]*y[,1]+a[2]*y[,2])
  }
  stat <- cbind(RSMS=g^2,SSMS=rowSums((y%*%geo$inv)*y))/(1+seq_len(nrow(y))/m)^2
  list(y=y,stat=stat)
}
# Exact drift normalization: standardized score direction at the break is
# (-1/sqrt(2),0). Both onsets use the same fixed loading amplitude.
make_h <- function(m) {
  t <- seq_len(3*m)/m
  cbind(cos(2*pi*t),sin(2*pi*t))
}
gaussian_run <- function(i,m,H,critical=NULL) {
  set.seed(if(is.null(critical)) 640000L+i else 650000L+i)
  pre <- seq_len(m)
  npost <- if(is.null(critical)) 2*m else m%/%2L
  e <- rnorm(m+npost)
  z <- H[seq_along(e),,drop=FALSE]*(e-mean(e[pre]))
  geo <- geometry(z[pre,,drop=FALSE])
  mon <- mon_path(z[-pre,,drop=FALSE],geo,m)
  if (i==1L) {
    idx <- unique(round(seq(1,nrow(mon$y),length.out=8)))
    hull <- aersn::aersn_hull(geo$path,n_directions=0)
    reported <- sqrt(mon$stat[idx,1])*(1+idx/m)
    stopifnot(max(abs(aersn::aersn_gauge(hull,mon$y[idx,,drop=FALSE])-reported))<1e-7)
  }
  maxima <- apply(mon$stat,2,max)
  if(is.null(critical)) return(maxima)
  # Distance along the drift ray until the training norm hits its boundary.
  pos <- mon$y[npost,]
  direction <- c(-1/sqrt(2),0)
  radius <- 1.5*sqrt(critical)
  A <- geo$normals
  forward <- as.numeric(A%*%direction)>0
  tR <- min((radius[1]-A[forward,,drop=FALSE]%*%pos)/(A[forward,,drop=FALSE]%*%direction))
  D <- geo$inv
  aa <- drop(crossprod(direction,D%*%direction))
  bb <- 2*drop(crossprod(direction,D%*%pos))
  cc <- drop(crossprod(pos,D%*%pos))-radius[2]^2
  tS <- (-bb+sqrt(max(0,bb^2-4*aa*cc)))/(2*aa)
  distance <- c(tR,tS)
  distance[maxima>=critical] <- NA_real_
  c(distance,sqrt(2*pmax(distance,0)))
}
finite_run <- function(i,m,H,critical) {
  set.seed(660000L+10000L*match(m,sizes)+i)
  n <- 3*m
  pre <- seq_len(m)
  monitor <- seq.int(m+1L,n)
  f <- rnorm(n)
  x1 <- f+rnorm(n)
  e2 <- rnorm(n)
  # The reference eigenvector is fitted on the actual Gaussian training panel.
  pc <- eigen(crossprod(cbind(x1[pre],e2[pre]))/m,symmetric=TRUE)
  u <- pc$vectors[,1]*sign(pc$vectors[1,1])
  v <- c(-u[2],u[1])
  scores <- function(ids,x2) {
    (x1[ids]*u[1]+x2*u[2])*(x1[ids]*v[1]+x2*v[2])*mean(v)/sqrt(pc$values[1])
  }
  geo <- geometry(H[pre,,drop=FALSE]*scores(pre,e2[pre]))
  post_age <- (seq_along(monitor)-m/2)/m
  out <- numeric(4)
  for(nu in 0:1) {
    shape <- if(nu==0) as.numeric(post_age>0) else pmax(0,post_age)
    x2 <- e2[monitor]+amplitude*shape*f[monitor]
    stat <- mon_path(H[monitor,,drop=FALSE]*scores(monitor,x2),geo,m)$stat
    for(j in 1:2) {
      hits <- which(stat[,j]>critical[j])
      k <- if(length(hits)) hits[1L] else Inf
      out[nu*2L+j] <- k
    }
  }
  out
}
# Independent seed per draw allows platform-independent batching.
batch <- function(B,fun) {
  groups <- split(seq_len(B),rep(seq_len(min(B,workers)),length.out=B))
  work <- function(ids) do.call(rbind,lapply(ids,fun))
  out <- if(workers==1L) lapply(groups,work) else parallel::mclapply(groups,work,
      mc.cores=workers,mc.set.seed=FALSE)
  stopifnot(all(vapply(out,is.matrix,logical(1))))
  do.call(rbind,out)[order(unlist(groups)),,drop=FALSE]
}
# All full-mode streams are independent across stages. Within each finite
# replication the methods and both shape alternatives use the same panel.
args <- commandArgs(trailingOnly=TRUE)
quick <- "--quick" %in% args
source_root <- dirname(normalizePath(sub("^--file=","",
    grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)[1])))
dest <- Sys.getenv("DFT_RESULTS",unset="results")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
dest <- normalizePath(dest)
archive_name <- "delay_draws.rds"

if (quick) {
  candidates <- c(file.path(dest,archive_name),
                  file.path(source_root,"reference",archive_name),
                  file.path(source_root,archive_name))
  archive <- candidates[file.exists(candidates)][1]
  if (is.na(archive)) stop("Cannot find the archived delay_draws.rds.")
  saved <- readRDS(archive)
  stopifnot(identical(saved$schema,1L))
} else {
  H <- make_h(grid)
  cal <- batch(C,function(i) gaussian_run(i,grid,H))
  crit <- apply(cal,2,quantile,.95,type=8)
  ray <- batch(2L*B,function(i) gaussian_run(i,grid,H,crit))
  stops <- list()
  for (m in sizes) {
    H <- make_h(m)
    stops[[as.character(m)]] <- batch(B,function(i) finite_run(i,m,H,crit))
    cat("Finished finite m =",m,"\n"); flush.console()
  }
  saved <- list(schema=1L,critical=crit,calibration=cal,ray=ray,stops=stops)
  saveRDS(saved,file.path(dest,archive_name))
}
crit <- saved$critical
write.csv(data.frame(method=names(crit),critical=unname(crit),draws=C,grid=grid),
          file.path(dest,"delay_critical.csv"),row.names=FALSE)
runs <- list(Limit=saved$ray)
for (m in sizes) {
  k <- saved$stops[[as.character(m)]]
  delay <- k - m/2
  delay[delay<=0] <- NA_real_
  for (nu in 0:1) {
    n_m <- m*(sqrt(m)*amplitude)^(-1/(nu+1))
    cols <- 2L*nu + (1:2)
    delay[,cols] <- delay[,cols]/n_m
  }
  runs[[as.character(m)]] <- delay
}
rows <- list()
for(level in names(runs)) for(nu in 0:1) for(j in 1:2) {
  z <- runs[[level]][,nu*2L+j]
  surv <- z[!is.na(z)]
  # Inf is retained for a missed change. Distribution-free order statistic
  # confidence interval for the survivor median, conditional on frozen boundaries.
  n <- length(surv)
  q <- quantile(surv,c(.25,.5,.75),type=1,names=FALSE)
  ci <- sort(surv)[pmin(n,pmax(1,c(qbinom(.025,n,.5),qbinom(.975,n,.5)+1)))]
  rows[[length(rows)+1L]] <- data.frame(grid=level,onset=c("Jump","Linear")[nu+1],
       method=c("RSMS","SSMS")[j],B=length(z),survivors=n,early=mean(is.na(z)),
       miss_among_survivors=mean(is.infinite(surv)),q25=q[1],median=q[2],q75=q[3],
       median_lo=ci[1],median_hi=ci[2])
}
tab <- do.call(rbind,rows)
write.csv(tab,file.path(dest,"delay_results.csv"),row.names=FALSE)

# The four columns of the display use conditional medians and their Monte
# Carlo intervals. Full-precision quartiles, early and miss rates are in CSV.
row <- function(level) {
  x <- tab[tab$grid==level,]
  stopifnot(nrow(x)==4L)
  cells <- sprintf("%.2f [%.2f, %.2f]",x$median,x$median_lo,x$median_hi)
  paste0(paste(c(level,cells),collapse=" & ")," \\\\")
}
tex <- c("\\begin{table}[!ht]\\centering\\color{revisionred}",
  "\\caption{Median scaled KS detection delay conditional on survival to the break.}\\label{tab:et_delay_sim}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular}{@{}rcccc@{}}\\toprule",
  "$m$ & \\multicolumn{2}{c}{Jump} & \\multicolumn{2}{c}{Linear onset} \\\\",
  "& RSMS & SSMS & RSMS & SSMS \\\\ \\midrule",
  vapply(c(as.character(sizes),"Limit"),row,character(1)),
  "\\bottomrule\\end{tabular}",
  "\\par\\smallskip\\begin{minipage}{0.97\\linewidth}Intervals give 95\\% Monte Carlo order-statistic bounds for the median at the simulated critical values. Delay is divided by $m(\\sqrt m\\,\\delta)^{-1/(\\nu+1)}$, with $\\nu=0$ for a jump, $\\nu=1$ for a linear onset and $\\delta=2$. A path that survives to the break and has no alarm by the horizon contributes $+\\infty$.\\end{minipage}",
  "\\end{table}")
writeLines(tex,file.path(dest,"delay_table.tex"))
