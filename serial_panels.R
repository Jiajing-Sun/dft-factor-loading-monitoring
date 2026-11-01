# Fixed design grid: the original 4-dimensional Gaussian panel and three normalized changes.
# All methods share each fitted PCA. Calibration and evaluation have independent RNG streams.
source("setup.R")
source(file.path(ROOT,"monitoring_core.R"))
setwd(OUTPUT)
RNGkind("Mersenne-Twister","Inversion","Rejection")
m<-120L; n<-360L; N<-40L; k0<-60L; C<-2000L; B<-1500L
hs<-c(2L,4L,8L,16L)
methods<-c("RSMS","SSMS",paste0("HAC",hs))
rr<-(1:n)/m; H<-cbind(cos(2*pi*rr),sin(2*pi*rr),cos(4*pi*rr),sin(4*pi*rr))
L<-sin(2*pi*(1:N)/N)+.5; U<-L/sqrt(sum(L^2)); a<-1-U*sum(U); a<-a/sqrt(mean(a^2))
d<-1:180; g<-cbind(oscillation=cos(2*pi*(m+k0+d)/m),pulse=as.numeric(d<=45),reversal=ifelse(d<=90,1,-1))
g<-rbind(matrix(0,m+k0,3),sweep(g,2,sqrt(colMeans(g^2)),"/"))
profiles<-c("flat","front","back","rising","falling")
grid<-expand.grid(profile=profiles,rho=c(0,.5,.85),stringsAsFactors=FALSE)
run<-function(id,kernel_study=FALSE){
 profile<-grid$profile[id];rho<-grid$rho[id]
 v<-switch(profile,flat=rep(1,n),front=1+.65*sin(2*pi*rr),back=1-.65*sin(2*pi*rr),rising=exp(.35*(rr-.5)),falling=exp(-.35*(rr-.5)))
 v<-v/mean(v[1:m])
 draw<-function(alternatives,check=FALSE){
   z<-matrix(rnorm(n*(N+1)),n,N+1)
   # Stationary unit-variance Gaussian AR(1), independent across factor/error streams.
   if(rho!=0){z[-1,]<-sqrt(1-rho^2)*z[-1,];z<-matrix(as.numeric(stats::filter(z,rho,method="recursive")),n,N+1L)}
   f<-z[,1];X<-(tcrossprod(f,L)+z[,-1])*sqrt(v)
   pc<-eigen(crossprod(X[1:m,])/m,symmetric=TRUE);u<-pc$vectors[,1];ev<-pc$values[1]
   xp<-drop(X%*%u);xb<-rowMeans(X)
   scores<-function(amp){p<-xp+amp*sum(a*u);H*(p/sqrt(ev)*(xb+amp*mean(a)-p*mean(u)))}
   psi<-scores(0); tr<-psi[1:m,]; mu<-colMeans(tr)
# The training geometry is fitted once per panel and reused for every change shape.
   G<-aersn::aersn_path(tr)$G
   R<-solve(qr.R(qr(G)));GG<-G%*%R
   V<-GG[unique(as.vector(geometry::convhulln(GG))),,drop=FALSE]
   ij<-expand.grid(i=seq_len(nrow(V)),j=seq_len(nrow(V)))
   facets<-geometry::convhulln(unique(V[ij$i,]-V[ij$j,]),options="Qt",output.options="n")$normals
   A<-R%*%t(facets[,1:4]/(-facets[,5]))
   W<-sweep(tr,2,mu);base<-crossprod(W)/m
   # Bartlett bandwidth h gives positive weights to lags 1 through h-1.
   covs<-lapply(1:15,function(l){x<-crossprod(W[(l+1):m,],W[1:(m-l),])/m;x+t(x)})
   inv<-c(list(solve(crossprod(G)/m)),lapply(hs,function(h){
     S<-base;for(l in seq_len(h-1L))S<-S+(1-l/h)*covs[[l]]
     solve(S+diag(1e-10,4))
   }))
   if (kernel_study) {
     # Match support at each bandwidth: both kernels use lags 1 through h-1.
     parzen <- lapply(hs,function(h) {
       S <- base
       weights <- sandwich::kweights(seq_len(h-1L)/h,kernel="Parzen")
       for (lag in seq_along(weights)) S <- S + weights[lag]*covs[[lag]]
       solve(S+diag(1e-10,4))
     })
     inv <- c(inv,parzen)
   }
   path<-function(s){
     Z<-apply(sweep(s[-(1:m),],2,mu),2,cumsum)/sqrt(m)/(1+(1:(n-m))/m)
     cbind(apply(Z%*%A,1,max)^2,vapply(inv,function(V)rowSums((Z%*%V)*Z),numeric(n-m)))
   }
   null<-path(psi)
   if(check)stopifnot(max(abs(null[,1:2]-seqstats(psi,m,TRUE)))<1e-5)
   if (kernel_study) {
     if (!alternatives) return(cbind(KS=apply(null,2,max),CvM=colSums(null)/m))
     return(list(null=null,oscillation=path(scores(.12*sqrt(v)*f*g[,"oscillation"]))))
   }
   if(!alternatives)return(apply(null,2,max))
   out<-array(0,c(n-m,6,4));out[,,1]<-null
   for(j in 1:3)out[,,j+1]<-path(scores(.12*sqrt(v)*f*g[,j]))
   out
 }
 # Calibrate on null panels before drawing the independent evaluation sample.
 if (kernel_study) {
   source(file.path(ROOT,"kernel_draws.R"),local=TRUE)
   return(id)
 }
 set.seed(813000+id);cal<-t(replicate(C,draw(FALSE)));cuts<-apply(cal,2,quantile,.95,type=8)
 stops<-array(Inf,c(B,6,4));pre<-post<-array(0,c(B,6,4))
 set.seed(913000+id)
 for(i in 1:B){s<-draw(TRUE,i==1)
   pre[i,,]<-apply(s[1:k0,,,drop=FALSE],2:3,max)
   post[i,,]<-apply(s[-(1:k0),,,drop=FALSE],2:3,max)
   for(j in 1:4)for(k in 1:6){h<-which(s[,k,j]>cuts[k]);if(length(h))stops[i,k,j]<-h[1]}
 }
 saveRDS(list(cal=cal,stop=stops,pre=pre,post=post,cut=cuts),sprintf("cell_%02d.rds",id))
 id
}
{
if ("--kernel" %in% commandArgs(TRUE)) {
  run(which(grid$profile=="front" & grid$rho==0.5),kernel_study=TRUE)
  quit(status=0)
}
ids<-as.integer(commandArgs(TRUE))
if(!length(ids)) ids<-seq_len(nrow(grid))
jobs<-parallel::mclapply(ids,run,mc.cores=if(.Platform$OS.type=="windows") 1L else as.integer(Sys.getenv("DFT_CORES","7")),mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(all(vapply(jobs,is.integer,logical(1))))
}
