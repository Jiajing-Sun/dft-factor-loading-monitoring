source("setup.R")
setwd(OUTPUT)
# Each resample preserves pairing across methods and alternatives; calibration
# draws and evaluation draws are resampled independently. Bands cover the 45
# pre-specified RSMS-minus-default-HAC power contrasts jointly.
profiles<-c("flat","front","back","rising","falling")
grid<-expand.grid(profile=profiles,rho=c(0,.5,.85),stringsAsFactors=FALSE)
methods<-c("RSMS","SSMS","HAC2","HAC4","HAC8","HAC16")
shapes<-c("oscillation","pulse","reversal")
RNGkind("Mersenne-Twister","Inversion","Rejection")
bootN<-999L
one<-function(id){
 x<-readRDS(sprintf("cell_%02d.rds",id));B<-dim(x$stop)[1];C<-nrow(x$cal)
 base<-readRDS(sprintf("cell_%02d.rds",match(grid$profile[id],profiles)))$cut
 tab<-list()
 for(j in 1:3)for(k in 1:6){
   t<-x$stop[,k,j+1];nul<-pmax(x$pre[,k,1],x$post[,k,1]);hit<-is.finite(t)&t>60
   tab[[length(tab)+1L]]<-data.frame(grid[id,],shape=shapes[j],method=methods[k],cut=x$cut[k],
      size=100*mean(nul>x$cut[k]),iid_cut_size=100*mean(nul>base[k]),power=100*mean(hit),
      early=100*mean(t<=60),median_delay=median(t[hit]-60),restricted_delay=mean(pmin(t[t>60]-60,181)))
 }
 rates<-function(ii,cv) {
  p<-vapply(1:6,function(k)colMeans(x$pre[ii,k,2:4]<=cv[k]&x$post[ii,k,2:4]>cv[k]),numeric(3))
  c(p[,1]-p[,4],p[,1]-p[,2],p[,3]-p[,6])*100
 }
 est<-rates(1:B,x$cut)
 set.seed(1230010+id)
 bs<-replicate(bootN,{
   cv<-apply(x$cal[sample.int(C,C,TRUE),,drop=FALSE],2,quantile,.95,type=8)
   rates(sample.int(B,B,TRUE),cv)
 })
 list(tab=do.call(rbind,tab),est=est,bs=bs)
}
res<-parallel::mclapply(1:15,one,mc.cores=if(.Platform$OS.type=="windows") 1L else as.integer(Sys.getenv("DFT_CORES","5")),mc.set.seed=FALSE)
all<-do.call(rbind,lapply(res,function(z)z$tab))
write.csv(all,"all_results.csv",row.names=FALSE)
con<-list()
for(k in 1:3){
  idx<-(k-1)*3+1:3
  dev<-do.call(rbind,lapply(res,function(z)sweep(z$bs[idx,,drop=FALSE],1,z$est[idx])))
  width<-unname(quantile(apply(abs(dev),2,max),.95))
  for(i in 1:15)for(j in 1:3){
    item<-res[[i]];at<-idx[j]; ci<-quantile(item$bs[at,],c(.025,.975))
    con[[length(con)+1]]<-data.frame(grid[i,],shape=shapes[j],contrast=c("RSMS-HAC4","RSMS-SSMS","HAC2-HAC16")[k],
       gain=item$est[at],low=ci[1],high=ci[2],joint_low=item$est[at]-width,joint_high=item$est[at]+width)
  }
}
write.csv(do.call(rbind,con),"paired_power_intervals.csv",row.names=FALSE)
# One monochrome page per variance profile, with shared axes across bandwidths.
pdf("bandwidth_power.pdf",width=8,height=7)
for(profile in profiles){
 par(mfrow=c(3,3),mar=c(2.6,2.8,2,0.3),mgp=c(1.6,.5,0),oma=c(1.7,1.5,2,0))
 for(rho in c(0,.5,.85))for(shape in shapes){
  z<-all[all$profile==profile & all$rho==rho & all$shape==shape,]
  plot(1:4,z$power[3:6],type="b",pch=16,xaxt="n",ylim=c(0,100),xlab="",ylab="",main=paste(shape," / ",rho),cex.main=.9)
  axis(1,1:4,c("2","4","8","16"));abline(h=z$power[1],lty=2,lwd=1.5);abline(h=z$power[2],lty=3)
 }
 mtext(paste(profile," | points: HAC, dashed: RSMS, dotted: SSMS"),outer=TRUE,side=3,line=.5)
 mtext("Bartlett bandwidth h",outer=TRUE,side=1);mtext("Size-adjusted post-break detection (%)",outer=TRUE,side=2)
}
dev.off()
