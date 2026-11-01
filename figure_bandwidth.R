# Main Figure 2: compare separately size-adjusted monitors on the same panels.
source("setup.R")
tab <- read.csv(file.path(OUTPUT,"all_results.csv"))
cell <- subset(tab, profile=="front" & rho==0.5 & shape=="oscillation")
cell <- cell[match(c("RSMS","SSMS","HAC2","HAC4","HAC8","HAC16"),cell$method),]
stopifnot(nrow(cell)==6L, !anyNA(cell$method))
pdf(file.path(OUTPUT,"front_cycle_bandwidth.pdf"),width=7.2,height=3.8,useDingbats=FALSE)
layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(9,1))
par(mar=c(3.7,3.8,2,0.7),oma=rep(0,4),mgp=c(2.2,.6,0))
for (v in c("power","restricted_delay")){
  y <- cell[[v]]
  detection <- v=="power"
  plot(1:4,y[3:6],type="n",xlab="Bartlett bandwidth",xaxt="n",las=1,
       ylim=if(detection)c(30,70) else c(115,155),
       ylab=if(detection)"Detection (%)" else "Delay (observations)",
       main=if(detection)"(a) Size-adjusted power" else "(b) Restricted mean delay",cex.main=.95)
  axis(1,1:4,c("2","4","8","16"))
  abline(h=y[1],lty=2,lwd=1.7); abline(h=y[2],lty=3,lwd=1.7)
  lines(1:4,y[3:6],type="b",pch=16,lty=1,lwd=1.4)
}
par(mar=rep(0,4))
plot.new()
legend("center",c("HAC","RSMS","SSMS"),lty=1:3,pch=c(16,NA,NA),lwd=1.5,horiz=TRUE,bty="n",cex=.85,inset=0)
dev.off()
