
library(Iso)

toxeffprob <- read.csv("./true prob/True probbaility.csv")
toxprob    <- toxeffprob[(toxeffprob$Item=="Tox")&(toxeffprob$Grade!="DLT")&(toxeffprob$Grade!="ETS")&(toxeffprob$Grade!="nETS"),]
effprob    <- toxeffprob[(toxeffprob$Item=="Eff")&(toxeffprob$Grade!="ORR")&(toxeffprob$Grade!="EES")&(toxeffprob$Grade!="nEES"),]

ls1   <- length(unique(toxprob$Scenario))
lds   <- 6
dosen <- 1:lds
dose  <- paste("Dose",1:lds,sep="")

ncoh <- 3
nmax <- 36
nesc <- nmax/ncoh

lambdaeta <- read.csv("mylambdaeta.csv")
phi       <- lambdaeta$phi
phi1      <- phi*0.1
phi2      <- phi*1.4
delta     <- lambdaeta$delta
delta1    <- delta*0.6
lambda1   <- lambdaeta$lambda1
lambda2   <- lambdaeta$lambda2
eta1      <- lambdaeta$eta1

weight.T <- c(0,0.5,1,1.5)
weight.E <- c(0,0.25,1,3)

tau.T <- 30
tau.E <- c(45,60,90)
ls2   <- length(tau.E)

accrual <- c(5,10)
ls3   <- length(accrual)

pr.alpha <- 1
pr.beta  <- 1

alpha.T1 <- 0.5
alpha.T2 <- 0.5
alpha.E1 <- 0.5
alpha.E2 <- 0.5

lmm <- 2
lss <- 10000

data.obs.n   <- array(0,dim=c(lmm,lss,ls1,ls2,ls3,lds))
data.obs.tox <- array(0,dim=c(lmm,lss,ls1,ls2,ls3,lds))
data.obs.eff <- array(0,dim=c(lmm,lss,ls1,ls2,ls3,lds))
data.dur     <- array(0,dim=c(lmm,lss,ls1,ls2,ls3))

for(mm in 1:lmm){
for(ss in 1:lss){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 efftoxp <- list(effp=effprob[effprob$Scenario==s1,dose],toxp=toxprob[toxprob$Scenario==s1,dose])

 obs.n   <- numeric(lds)
 obs.tox <- numeric(lds)
 obs.eff <- numeric(lds)
 pe      <- numeric(lds)
 pt      <- numeric(lds)

 t.enter    <- NULL
 t.decision <- 0

 curdose <- 1

 for(i in 1:nesc){

  dlab <- paste("Dose",curdose,sep="")
  obs.n[curdose] <- obs.n[curdose] + ncoh

  for(j in 1:ncoh){
   if(j==1){
    t.enter <- c(t.enter,t.decision)
   }else{
    t.enter <- c(t.enter,t.enter[length(t.enter)]+runif(1,0,2*accrual[s3]))
  }}
  t.decision <- t.enter[length(t.enter)]+max(tau.T,tau.E[s2])

  psi.T    <- sum(efftoxp$toxp[2:4,dlab])
  zetta.T1 <- log(log(1-psi.T)/log(1-psi.T+alpha.T1*psi.T))/log(1/(1-alpha.T2))
  zetta.T2 <- tau.T/(-log(1-psi.T))^(1/zetta.T1)
  time.tox <- rweibull(ncoh,shape=zetta.T1,scale=zetta.T2)
  event.T  <- as.numeric(time.tox<=tau.T)
  grade    <- event.T*((1:3)%*%rmultinom(3,1,efftoxp$toxp[2:4,dlab]))+1
  DLT      <- as.numeric(grade>=3)
  ETS      <- apply(grade,2,function(x){return(weight.T[x])})
  nETS     <- ETS/max(weight.T)

  psi.E    <- sum(efftoxp$effp[2:4,dlab])
  zetta.E1 <- log(log(1-psi.E)/log(1-psi.E+alpha.E1*psi.E))/log(1/(1-alpha.E2))
  zetta.E2 <- tau.E[s2]/(-log(1-psi.E))^(1/zetta.E1)
  time.eff <- rweibull(ncoh,shape=zetta.E1,scale=zetta.E2)
  event.E  <- as.numeric(time.eff<=tau.E[s2])
  response <- event.E*((1:3)%*%rmultinom(3,1,efftoxp$effp[2:4,dlab]))+1
  ORR      <- as.numeric(response>=3)
  EES      <- apply(response,2,function(x){return(weight.E[x])})
  nEES     <- EES/max(weight.E)

  if(mm==1){ obs.tox[curdose] <- obs.tox[curdose]+sum(DLT)  }
  if(mm==2){ obs.tox[curdose] <- obs.tox[curdose]+sum(nETS) }
  pt[curdose] <- obs.tox[curdose]/obs.n[curdose]

  if(mm==1){ obs.eff[curdose] <- obs.eff[curdose]+sum(ORR)  }
  if(mm==2){ obs.eff[curdose] <- obs.eff[curdose]+sum(nEES) }
  pe[curdose] <- obs.eff[curdose]/obs.n[curdose]

  if((pt[curdose]<=lambda1[mm])&(pe[curdose]<=eta1[mm])){
   nxtdose <- curdose+1
  }else if((pt[curdose]<lambda2[mm])&(pe[curdose]>eta1[mm])){
   nxtdose <- curdose
  }else if(pt[curdose]>=lambda2[mm]){
   nxtdose <- curdose-1

  }else if((pt[curdose]>lambda1[mm])&(pt[curdose]<lambda2[mm])&(pe[curdose]<=eta1[mm])){

   if(curdose==lds){
    three   <- c(curdose-1,curdose)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)

   }else if(obs.n[curdose+1]==0){
    nxtdose <- curdose+1

   }else if(curdose==1){
    three   <- c(curdose,curdose+1)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)

   }else{
    three   <- c(curdose-1,curdose,curdose+1)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)
   }

  }

  tpo.shape1 <- pr.alpha + obs.tox
  tpo.shape2 <- pr.beta  + (obs.n-obs.tox)
  tterm      <- pbeta(phi[mm],tpo.shape1,tpo.shape2)

  epo.shape1 <- pr.alpha + obs.eff
  epo.shape2 <- pr.beta  + (obs.n-obs.eff)
  eterm      <- 1-pbeta(delta1[mm],epo.shape1,epo.shape2)

  admflg  <- !((eterm<0.05)|(tterm<0.05))
  admdose <- dosen[admflg]

  if(sum(admflg)==0){
   break
  }else{

   if(nxtdose==0){
    if(admflg[1]){
     curdose <- 1
    }else{
     break
    }

   }else if(nxtdose==(lds+1)){
    curdose <- lds

   }else if(is.element(nxtdose,admdose)){
    curdose <- nxtdose

   }else if(curdose<nxtdose){
    if(sum(admdose>=nxtdose)!=0){
     curdose <- min(admdose[admdose>=nxtdose])
    }
   }else if(curdose>=nxtdose){
    if(sum(admdose<=nxtdose)!=0){
     curdose <- max(admdose[admdose<=nxtdose])
    }else{
     break
    }
   }

  }

 }

 data.obs.n[mm,ss,s1,s2,s3,]   <- obs.n
 data.obs.tox[mm,ss,s1,s2,s3,] <- obs.tox
 data.obs.eff[mm,ss,s1,s2,s3,] <- obs.eff
 data.dur[mm,ss,s1,s2,s3]      <- t.decision

}}}}}

#####

futility3 <- function(pe,pt,psi00,psi11)
{
 psi.e0t1 <- 0
 psi.e0t0 <- psi00
 psi.e1t1 <- psi11
 psi.e1t0 <- 100

 ut = (  psi.e0t1*(1-pe)*pt
       + psi.e0t0*(1-pe)*(1-pt)
       + psi.e1t1*pe    *pt
       + psi.e1t0*pe    *(1-pt))

 return(ut)
}

mydose1 <- function(pevec,ptvec,ph,ph2,probe,probt)
{
 candose <- which((probe>=0.05)&(probt>=0.05))
 adddose <- which((probe>=0.05)&(probt>=0.05)&(ptvec<=ph2))

 fu31 <- futility3(pe = pevec,
                   pt = ptvec,
                   psi00=40,psi11=60)
 fu31.max <- which(fu31==max(fu31[candose]))

 re311 <- min(intersect(fu31.max,candose))

 if(length(adddose)==0){
  re312 <- 0
 }else{
  re312.max <- which(fu31==max(fu31[adddose]))
  re312 <- min(intersect(re312.max,adddose))
 }

 mdif1 <- min(abs(ptvec[candose]-ph))
 mtd1  <- max(which(abs(ptvec-ph)==mdif1))
 meff1 <- max(pevec[intersect(candose,1:mtd1)])
 deff1 <- which(pevec==meff1)
 re41  <- min(intersect(candose,deff1))

 if(length(adddose)==0){
  re42 <- 0
 }else{
  mdif2 <- min(abs(ptvec[adddose]-ph))
  mtd2  <- max(which(abs(ptvec-ph)==mdif2))
  meff2 <- max(pevec[intersect(adddose,1:mtd2)])
  deff2 <- which(pevec==meff2)
  re42  <- min(intersect(adddose,deff2))
 }

 rval <- c(re311,re312,re41,re42)

 return(rval)
}

bi.pr.alpha <- 1
bi.pr.beta  <- 1

res <- array(0,dim=c(lmm,lss,ls1,ls2,ls3,4))

for(mm in 1:lmm){
for(ss in 1:lss){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 bi.prob.pe <- numeric(lds)
 bi.prob.pt <- numeric(lds)

 obsn <- data.obs.n[mm,ss,s1,s2,s3,]
 obst <- data.obs.tox[mm,ss,s1,s2,s3,]
 obse <- data.obs.eff[mm,ss,s1,s2,s3,]
 evadose <- dosen[data.obs.n[mm,ss,s1,s2,s3,]!=0]

 obspe <- obse[evadose]/obsn[evadose]
 obspt <- obst[evadose]/obsn[evadose]

 for(i in evadose){
  po.shape1 <- bi.pr.alpha + obse[i]
  po.shape2 <- bi.pr.beta  + (obsn[i]-obse[i])
  bi.prob.pe[i] <- 1-pbeta(delta1[mm],po.shape1,po.shape2)

  po.shape1 <- bi.pr.alpha + obst[i]
  po.shape2 <- bi.pr.beta  + (obsn[i]-obst[i])
  bi.prob.pt[i] <- pbeta(phi[mm],po.shape1,po.shape2)
 }

 if(sum(obsn)<nmax){

   res[mm,ss,s1,s2,s3,] <- 0

 }else if(length(evadose)==1){

  if((bi.prob.pe[evadose]>=0.05)&(bi.prob.pt[evadose]>=0.05)){
   res[mm,ss,s1,s2,s3,] <- evadose
  }

 }else if(sum((bi.prob.pe[evadose]>=0.05)&(bi.prob.pt[evadose]>=0.05))>=1){

  iso.obspt <- pava(obspt[evadose])
  res[mm,ss,s1,s2,s3,] <- mydose1(obspe,iso.obspt,phi[mm],phi2[mm],bi.prob.pe[evadose],bi.prob.pt[evadose])

 }

}}}}}

sprop <- array(0,dim=c(lmm,ls1,ls2,ls3,4,7))

for(mm in 1:lmm){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 sprop[mm,s1,s2,s3,,1] <- apply(res[mm,,s1,s2,s3,]==0,2,mean)*100
 sprop[mm,s1,s2,s3,,2] <- apply(res[mm,,s1,s2,s3,]==1,2,mean)*100
 sprop[mm,s1,s2,s3,,3] <- apply(res[mm,,s1,s2,s3,]==2,2,mean)*100
 sprop[mm,s1,s2,s3,,4] <- apply(res[mm,,s1,s2,s3,]==3,2,mean)*100
 sprop[mm,s1,s2,s3,,5] <- apply(res[mm,,s1,s2,s3,]==4,2,mean)*100
 sprop[mm,s1,s2,s3,,6] <- apply(res[mm,,s1,s2,s3,]==5,2,mean)*100
 sprop[mm,s1,s2,s3,,7] <- apply(res[mm,,s1,s2,s3,]==6,2,mean)*100

}}}}

#####

method <- c("BOIN-ET","gBOIN-ET")
nobd   <- as.numeric(apply(toxeffprob[toxeffprob$Grade=="Utility",dose],1,which.max))

obdm <- c(1,3)
ls4  <- length(obdm)

pcs.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))
pca.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))
pos.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))
poa.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))
pet.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))
dur.s1 <- array(0,dim=c(lmm,ls1,ls2,ls3,ls4))

for(mm in 1:lmm){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){

 pcs.s1[mm,s1,s2,s3,s4] <- mean(res[mm,,s1,s2,s3,obdm[s4]]==nobd[s1])*100
 pca.s1[mm,s1,s2,s3,s4] <- sum(data.obs.n[mm,,s1,s2,s3,nobd[s1]])/sum(data.obs.n[mm,,s1,s2,s3,])*100
 if(nobd[s1]<max(dosen)){
  pos.s1[mm,s1,s2,s3,s4] <- mean((res[mm,,s1,s2,s3,obdm[s4]]!=0)&(res[mm,,s1,s2,s3,obdm[s4]]>nobd[s1]))*100
  poa.s1[mm,s1,s2,s3,s4] <- sum(data.obs.n[mm,,s1,s2,s3,(nobd[s1]+1):lds])/sum(data.obs.n[mm,,s1,s2,s3,])*100
 }
 pet.s1[mm,s1,s2,s3,s4] <- mean(res[mm,,s1,s2,s3,obdm[s4]]==0)*100
 dur.s1[mm,s1,s2,s3,s4] <- mean(data.dur[mm,,s1,s2,s3])

}}}}}

obdml <- c("Utility","Highest Efficacy")

out.df <- NULL
for(mm in 1:lmm){
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){

 out    <- cbind(pcs.s1[mm,,s2,s3,s4],pca.s1[mm,,s2,s3,s4],pos.s1[mm,,s2,s3,s4],poa.s1[mm,,s2,s3,s4],pet.s1[mm,,s2,s3,s4],dur.s1[mm,,s2,s3,s4])
 out.df <- rbind(out.df,data.frame(M  = method[mm],
                                   Z1 = tau.E[s2],
                                   Z2 = accrual[s3],
                                   Z3 = obdml[s4],
                                   Z4 = 1:ls1,
                                   data.frame(out)))
}}}}

write.csv(out.df,".\\output\\Table1_BOIN-ET_gBOIN-ET.csv")

out.df <- NULL
for(mm in 1:lmm){
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){
for(s1 in 1:ls1){

 out <- rbind(sprop[mm,s1,s2,s3,obdm[s4],],c(NA,apply(data.obs.n[mm,,s1,s2,s3,],2,mean)))
 out.df <- rbind(out.df,data.frame(M  = method[mm],
                                   Z1 = tau.E[s2],
                                   Z2 = accrual[s3],
                                   Z3 = obdml[s4],
                                   Z4 = s1,
                                   Z5 = c("Rec","Pats"),
                                   data.frame(out)))
}}}}}

write.csv(out.df,".\\output\\TableS2_BOIN-ET_gBOIN-ET.csv")



