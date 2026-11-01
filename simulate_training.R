# Training replacements in a one-factor panel. Run with Rscript simulate_training.R.
# The same innovations generate the clean and contaminated panels in each draw.
# Use RSMS_CORE to locate the released monitoring_core.R if it is elsewhere.
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])))
source(Sys.getenv("RSMS_CORE", file.path(here, "monitoring_core.R")))
stopifnot(as.character(utils::packageVersion("aersn")) == "0.2.3")
# Published layout: all R sources are at root; generated files go in results/
# and saved calibration/evaluation arrays in reference/ support --cached.
args <- commandArgs(TRUE)
output <- Sys.getenv("DFT_RESULTS", file.path(here,"results"))
dir.create(output, recursive=TRUE, showWarnings=FALSE)
reference <- Sys.getenv("DFT_TRAINING_REFERENCE", file.path(here,"reference"))
RNGkind("Mersenne-Twister","Inversion","Rejection")

N <- 40L
B <- 2000L
C <- 4000L
grid <- c(100L,400L,1600L)
designs <- c("clean","subcritical","boundary")
methods <- c("RSMS","SSMS")
loading <- sin(2*pi*seq_len(N)/N) + .5
U <- loading/sqrt(sum(loading^2))
b <- rep(1/N,N)
bperp <- b-U*sum(U*b)
direction <- bperp/sqrt(mean(bperp^2))
Gamma <- sum(loading^2) + 1
# Under Gaussian factors/errors the baseline score variance is ||(I-P)b||^2.
# Choose standardized P(C1-C0)=pi: for H=(cos,sin)(2*pi*r), the
# largest norm of the boundary drift over the first cycle is exactly 1.
sd_score <- sqrt(sum(bperp^2))
delta <- pi*sd_score/(sqrt(sum(loading^2)/Gamma)*sum(direction*b))
theta <- delta*sqrt(sum(loading^2)/Gamma)*sum(direction*b)/sd_score

run_training <- function(m, front=FALSE) {
  tag <- paste0("m",m,if(front) "_front" else "")
  stream <- m + if(front) 50000L else 0L
  cached <- "--cached" %in% args
  raw_name <- paste0("training_draws_",tag,".rds")
  n <- 3L*m
  H <- cbind(cos(2*pi*seq_len(n)/m), sin(2*pi*seq_len(n)/m))
  v <- if(front) 1+.65*sin(2*pi*seq_len(n)/m) else rep(1,n)
  v <- v/mean(v[seq_len(m)])
  # Fitted Gaussian reference for this specific variance profile. The
  # Gaussian increments and training-PCA correction both carry v(t).
  fitted_reference <- function(e) {
    if (!front) return(e-mean(e[seq_len(m)]))
    v*(e-sum(v[seq_len(m)]*e[seq_len(m)])/sum(v[seq_len(m)]))
  }
  # Under a flat profile, add a diffuse shift over every training date with
  # the same total projected covariance dose as the sparse boundary design.
  if (!front) designs <- c(designs, "diffuse")
  count <- c(0L,floor(m^.25),floor(sqrt(m)), if(!front) m else integer())
  dose <- count
  if(!front) dose[4] <- dose[3]
  replace_dates <- function(j) if(j == 0L) integer() else floor((seq_len(j)-.5)*m/j)+1L
  dates <- lapply(count, replace_dates)
  maximum <- function(psi) apply(seqstats(psi,m),2,max)
  arr <- function() array(NA_real_,c(B,length(designs),2),dimnames=list(NULL,designs,methods))
  # Independent fitted Gaussian calibration. Both methods see the same draws.
  if (cached) {
    draw <- readRDS(file.path(reference, raw_name))
    cal <- draw$cal; actual <- draw$actual; ref <- draw$reference
  } else {
  set.seed(840000L + stream)
  cal <- t(replicate(C, {
    e <- rnorm(n); maximum(H*fitted_reference(e))
  }))
  # Evaluate the clean and shifted reference on common Gaussian paths. The
  # training bridge is unchanged for the uniformly dispersed mean profile.
  # Subcritical finite-dose drift uses its actual n_rep/sqrt(m); its limit is 0.
  set.seed(850000L + stream)
  ref <- arr()
  for(i in seq_len(B)) {
    e <- rnorm(n)
    base <- fitted_reference(e)
    for(j in seq_along(designs)) {
      z <- base
      z[(m+1L):n] <- z[(m+1L):n] - theta*dose[j]/m*v[(m+1L):n]
      ref[i,j,] <- maximum(H*z)
    }
  }
  # Gaussian factor f and idiosyncratic errors are independent across dates.
  # A replacement changes the loading on that date to loading + delta*direction.
  # Sharing f and the errors gives the prescribed replacement covariance C1
  # and increases precision for clean-minus-contaminated comparisons.
  set.seed(860000L + stream)
  actual <- arr()
  for(i in seq_len(B)) {
    f <- rnorm(n)
    X0 <- (tcrossprod(f,loading) + matrix(rnorm(n*N),n,N))*sqrt(v)
    for(j in seq_along(designs)) {
      X <- X0
      t <- dates[[j]]
      amplitude <- if(designs[j] == "diffuse") delta*dose[j]/m else delta
      if(length(t)) X[t,] <- X[t,,drop=FALSE] + amplitude*tcrossprod(f[t]*sqrt(v[t]), direction)
      pc <- eigen(crossprod(X[seq_len(m),,drop=FALSE])/m,symmetric=TRUE)
      u <- pc$vectors[,1]
      xp <- as.numeric(X%*%u)
      score <- xp/sqrt(pc$values[1])*(rowMeans(X)-mean(u)*xp)
      actual[i,j,] <- maximum(H*score)
    }
  }
  saveRDS(list(m=m,N=N,B=B,C=C,delta=delta,theta=theta,cal=cal,actual=actual,
               reference=ref,profile=if(front) "front-cycle" else "flat",
               seeds=stream+c(840000,850000,860000,870000)),file.path(output,raw_name))
  }
  critical <- apply(cal, 2, quantile, .95, type=8)
  # Resample calibration draws separately; preserve pairing of the methods
  # and all three training conditions within actual/reference evaluation draws.
  rates <- function(arr,idx,crit) apply(sweep(arr[idx,,,drop=FALSE],3,crit,">"),2:3,mean)
  measure <- function(ic,ia,ir) {
    cv <- apply(cal[ic,,drop=FALSE],2,quantile,.95,type=8)
    a <- rates(actual,ia,cv); g <- rates(ref,ir,cv)
    ad <- sweep(a,2,a[1,],"-"); gd <- sweep(g,2,g[1,],"-")
    c(as.vector(a),as.vector(g),as.vector(a-g),as.vector(ad),as.vector(gd),as.vector(ad-gd),
      as.vector(sweep(a,2,a[3,],"-")))
  }
  est <- measure(seq_len(C),seq_len(B),seq_len(B))
  set.seed(870000L + stream)
  bs <- replicate(999L,measure(sample.int(C,C,TRUE),sample.int(B,B,TRUE),sample.int(B,B,TRUE)))
  ci <- t(apply(bs,1,quantile,probs=c(.025,.975),type=7))
  rows <- expand.grid(design=designs,method=methods,
                      measure=c("actual","reference","gap","actual_change","reference_change","change_gap","actual_vs_boundary"),
                      KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  rows$m <- m; rows$replaced <- count[match(rows$design,designs)]; rows$estimate <- est
  rows$mc_se <- apply(bs,1,sd); rows$lo <- ci[,1]; rows$hi <- ci[,2]
  rows$critical <- unname(critical[match(rows$method,methods)])
  rows$profile <- if(front) "front-cycle" else "flat"
  write.csv(rows,file.path(output,paste0("training_",tag,".csv")),row.names=FALSE)
  rows
}

cores <- as.integer(Sys.getenv("TRAINING_CORES","3"))
studies <- if("--front-only" %in% args) 4L else if("--flat-only" %in% args) 1:3 else seq_len(4L)
one <- function(i) run_training(c(grid,1600L)[i],front=(i==4L))
out <- if(.Platform$OS.type == "windows" || cores == 1L) lapply(studies,one) else
  parallel::mclapply(studies,one,mc.cores=min(cores,length(studies)),mc.set.seed=FALSE)
stopifnot(all(vapply(out,is.data.frame,logical(1))))

if (length(studies) == 4L) {
  all <- do.call(rbind,out)
  write.csv(all,file.path(output,"training_summary.csv"),row.names=FALSE)
  # The two estimates in each cell share a calibration sample; independent
  # paired panel and Gaussian evaluations produce the outside/inside numbers.
  lines <- c("\\begin{table}[!ht]\\centering\\color{revisionred}",
      "\\caption{Training contamination: panel and Gaussian-reference alarm probabilities (percent).}\\label{tab:contamination}",
      "\\setlength{\\tabcolsep}{3pt}\\begin{tabular}{llrrrr}\\toprule",
      "Profile, $m$ & Method & Clean & Sparse $m^{1/4}$ & Sparse $m^{1/2}$ & Diffuse \\\\ \\midrule")
  for(i in seq_len(4L)) for(k in seq_along(methods)) {
    cell <- all[all$profile==c(rep("flat",3),"front-cycle")[i] &
                       all$m==c(grid,1600L)[i] & all$method==methods[k],]
    vals <- vapply(c(designs,"diffuse"),function(j) {
      a <- cell$estimate[cell$design==j & cell$measure=="actual"]
      g <- cell$estimate[cell$design==j & cell$measure=="reference"]
      if(!length(a)) return("--")
      sprintf("%.2f (%.2f)",100*a,100*g)
    },character(1))
    profile_label <- if(i==4L) "Front" else "Flat"
    lab <- if(k==1L) paste0(profile_label,", ",c(grid,1600L)[i]) else ""
    lines <- c(lines,paste0(paste(c(lab,methods[k],vals),collapse=" & ")," \\\\"))
  }
  mc <- 100*max(all$mc_se[all$measure %in% c("actual","reference")])
  lines <- c(lines,"\\bottomrule\\end{tabular}",
    paste0("\\par\\smallskip\\begin{minipage}{0.97\\linewidth}\\footnotesize ",
    "Panel alarm frequencies are followed by the fitted Gaussian predictions in parentheses. ",
    "Sparse counts are rounded down. Diffuse and sparse $m^{1/2}$ changes have the same projected training dose. ",
    "Each profile and training length uses 4,000 independent calibration draws and 2,000 draws per evaluation sample. ",
    sprintf("Bootstrap Monte Carlo standard errors, including calibration uncertainty, are at most %.2f percentage points.",mc),
    "\\end{minipage}\\end{table}"))
  writeLines(lines,file.path(output,"contamination_table.tex"))
}
