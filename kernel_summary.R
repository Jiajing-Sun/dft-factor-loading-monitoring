# Compare rules at separately calibrated null cutoffs on paired evaluation paths.
source("setup.R")
x <- readRDS(file.path(OUTPUT,"kernel_front_rho05.rds"))
methods <- dimnames(x$stop)[[2]]
rules <- c("KS","CvM")
rows <- list()
for (rule in rules) for (method in methods) {
    null <- x$stop[,method,rule,"null"]
    alt <- x$stop[,method,rule,"oscillation"]
    detected <- is.finite(alt) & alt>60
    rows[[length(rows)+1]] <- data.frame(rule,method,size=100*mean(is.finite(null)),
        power=100*mean(detected),ADD=mean(alt[detected]-60),
        restricted=mean(pmin(alt[alt>60]-60,181)))
}
tab <- do.call(rbind,rows)
write.csv(tab,file.path(OUTPUT,"kernel_comparison.csv"),row.names=FALSE)
lines <- c("\\begin{tabular}{llrrrrrrrr}\\toprule",
    " & & \\multicolumn{4}{c}{KS} & \\multicolumn{4}{c}{CvM} \\\\",
    "Method & $h$ & Size & Power & ADD & RMD & Size & Power & ADD & RMD \\\\\\midrule")
for (method in methods) {
    d <- tab[tab$method==method,]
    h <- sub("^[A-Za-z]+","",method)
    label <- sub("[0-9]+$","",method)
    values <- unlist(lapply(1:2,function(i) c(sprintf("%.2f",unlist(d[i,c("size","power")])),
        sprintf("%.1f",unlist(d[i,c("ADD","restricted")])))))
    lines <- c(lines,paste(paste(c(label,if(nzchar(h)) h else "---",values),collapse=" & "),"\\\\"))
}
writeLines(c(lines,"\\bottomrule\\end{tabular}"),file.path(OUTPUT,"kernel_table.tex"))

# Resample shared calibration indices and, independently, shared evaluation
# indices. Each bootstrap draw re-estimates all method- and rule-specific cutoffs.
# Report all contrasts, including a simultaneous band across both rules.
left <- c(rep("RSMS",8),methods[3:6],"Bartlett2","Parzen2")
right <- c(methods[3:10],methods[7:10],"Bartlett16","Parzen16")
effects <- function(ii,cuts) {
    out <- matrix(NA_real_,length(left),2,dimnames=list(NULL,rules))
    for (r in rules) {
        p <- vapply(methods,function(j) mean(x$pre[ii,j,r,"oscillation"]<=cuts[j,r] &
                            x$post[ii,j,r,"oscillation"]>cuts[j,r]),numeric(1))
        out[,r] <- 100*(p[left]-p[right])
    }
    as.vector(out)
}
n <- dim(x$stop)[1]; nc <- dim(x$cal)[3]
estimate <- effects(seq_len(n),x$cut)
set.seed(2611017)
boot <- replicate(2000,{
    cut <- apply(x$cal[,,sample.int(nc,nc,TRUE),drop=FALSE],1:2,quantile,.95,type=8)
    effects(sample.int(n,n,TRUE),cut)
})
width <- unname(quantile(apply(abs(boot-estimate),2,max),.95))
interval <- t(apply(boot,1,quantile,c(.025,.975)))
write.csv(data.frame(rule=rep(rules,each=length(left)),left,right,estimate,
    lower=interval[,1],upper=interval[,2],joint_lower=estimate-width,joint_upper=estimate+width),
    file.path(OUTPUT,"kernel_paired_power.csv"),row.names=FALSE)
