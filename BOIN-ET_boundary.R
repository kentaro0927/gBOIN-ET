
# phi1 < lambda1 < phi =< lambda2 < phi2
# delta1 < eta1 < delta

mypinc <- function(pi,phi,phi1,phi2,delta,delta1,lambda1,lambda2,eta1)
{
 nsub <- 100

 l1p  <- pbinom(nsub*lambda1,nsub,phi)
 l1p1 <- pbinom(nsub*lambda1,nsub,phi1)
 l1p2 <- pbinom(nsub*lambda1,nsub,phi2)

 l2p  <- pbinom(nsub*lambda2-1,nsub,phi)
 l2p1 <- pbinom(nsub*lambda2-1,nsub,phi1)
 l2p2 <- pbinom(nsub*lambda2-1,nsub,phi2)

 e1d  <- pbinom(nsub*eta1,nsub,delta)
 e1d1 <- pbinom(nsub*eta1,nsub,delta1)

 pinc <- (pi[1]*(l1p1*(1-e1d1)             + 2/3*(l2p1-l1p1)*e1d1 + (l2p1-l1p1)*(1-e1d1) + (1-l2p1))
        + pi[2]*(l1p1*e1d                  + 2/3*(l2p1-l1p1)*e1d  + (1-l2p1))
        + pi[4]*(l1p*e1d                   + 2/3*(l2p-l1p1)*e1d   + (1-l2p))
        + pi[5]*(l1p2*e1d1 + l1p2*(1-e1d1) + 2/3*(l2p2-l1p2)*e1d1 + (l2p2-l1p2)*(1-e1d1))
        + pi[6]*(l1p2*e1d  + l1p2*(1-e1d)  + 2/3*(l2p2-l1p2)*e1d  + (l2p2-l1p2)*(1-e1d)))

 return(pinc)
}

s.phi   <- c(0.313,0.33)
lsph    <- length(s.phi)
s.delta <- c(0.583,0.70)
lsde    <- length(s.delta)

mylameta <- NULL

for(sph in 1:lsph){
for(sde in 1:lsde){

phi   <- s.phi[sph]
delta <- s.delta[sde]

l1  <- round(seq(phi*0.1,phi,by=0.01),digits=2)
l1s <- length(l1)
l2  <- round(seq(phi,phi*1.4,by=0.01),digits=2)
l2s <- length(l2)
e1  <- round(seq(delta*0.6,delta,by=0.01),digits=2)
e1s <- length(e1)

pincval <- array(0,dim=c(l1s,l2s,e1s))

for(s1 in 1:l1s){
for(s2 in 1:l2s){
for(s3 in 1:e1s){

 pincval[s1,s2,s3] <- mypinc(
                       pi      = rep(1/6,6),
                       phi     = phi,
                       phi1    = phi*0.1,
                       phi2    = phi*1.4,
                       delta   = delta,
                       delta1  = delta*0.6,
                       lambda1 = l1[s1],
                       lambda2 = l2[s2],
                       eta1    = e1[s3])

}}}

pnum <- which(pincval==min(pincval),arr.ind=TRUE)
ledf <- data.frame(phi=phi,delta=delta,lambda1=l1[pnum[1]],lambda2=l2[pnum[2]],eta1=e1[pnum[3]])

mylameta <- rbind(mylameta,ledf)

}}

mylameta
write.csv(mylameta,"lambdaeta.csv")


# opt = optim(
#         par     = c(l1[pnum[1]],l2[pnum[2]],e1[pnum[3]]),
#         fn      = mypinc,
#         method  = "L-BFGS-B",
#         lower   = c(phi1,phi,delta1),
#         upper   = c(phi,phi2,delta))


