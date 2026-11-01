# Main Figure 2 uses the same unrounded, size-adjusted powers as Table S31.
source("setup.R")
tab <- read.csv(file.path(OUTPUT, "kernel_comparison.csv"))
bandwidth <- c(2, 4, 8, 16)
kernels <- c("Bartlett", "Parzen")
rules <- c("KS", "CvM")

review <- "--review" %in% commandArgs(TRUE)
filename <- if (review) "front_cycle_bandwidth_review.pdf" else "front_cycle_bandwidth.pdf"
pdf(file.path(OUTPUT, filename), width=7.2, height=3.65,
    useDingbats=FALSE, title="Size-adjusted power by HAC bandwidth")
layout(matrix(c(1,2,3,3), 2, byrow=TRUE), heights=c(8.4,1))
par(mar=c(3.2,3.6,1.9,0.65), mgp=c(2.1,.55,0), tcl=-.25)
for (rule in rules) {
    d <- tab[tab$rule == rule, ]
    power <- setNames(d$power, d$method)
    stopifnot(!anyDuplicated(d$method))
    # Powers of two are evenly spaced on the log-bandwidth axis in both panels.
    plot(bandwidth, rep(NA_real_,4), type="n", log="x", xaxt="n", yaxt="n",
         ylim=c(34,66), xlab="", ylab="Size-adjusted power (%)",
         main=paste(if (rule=="KS") "(a)" else "(b)",rule), bty="l", las=1)
    title(xlab="HAC bandwidth h (in lags)", col.lab=if (review) "#8B0000" else "black")
    axis(1, at=bandwidth, labels=bandwidth)
    axis(2, at=seq(35,65,5), las=1)
    level <- power[["RSMS"]]
    stopifnot(is.finite(level))
    abline(h=level, lty=3, lwd=1.5)
    text(16,level+1.6,sprintf("%.1f%%",level),adj=1,cex=.82)
    for (j in seq_along(kernels)) {
        y <- power[paste0(kernels[j], bandwidth)]
        stopifnot(all(is.finite(y)))
        lines(bandwidth, y, lty=j, lwd=1.4)
    }
    # Open triangles remain visible around the small solid circles at close values.
    points(bandwidth, power[paste0("Parzen",bandwidth)], pch=2, cex=1.18, lwd=1.3)
    points(bandwidth, power[paste0("Bartlett",bandwidth)], pch=16, cex=.64)
}
par(mar=rep(0,4))
plot.new()
legend("center",c("Bartlett HAC","Parzen HAC","RSMS"),lty=1:3,pch=c(16,2,NA),
       pt.cex=c(.68,1.15,1),lwd=1.4,horiz=TRUE,bty="n",cex=.9)
dev.off()
