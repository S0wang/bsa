########################################################
##   example code for BSA for unmeasured confouding   ##
##   developed by: Shuang Wang                        ##
##   last updated: 12 Dec 2021                        ##
########################################################

##----------BSA-------------##
library(mvtnorm)
library(mcmcse)
library(survival)
library(GIGrvg)
library(heavy)
library(statmod)
library(pracma)
library(SuppDists)
library(dplyr)


############################
####   simulated data   ####
############################

## setup bias parameters and covariate coeff. in both treatment and outcome models
set.seed(1228)
n=1000 # 15000
lambda.a.list = c(-0.75,-0.5,-0.1,0,0.1,0.5,0.75)
lambda.s.list = c(-0.75,-0.5,-0.1,0,0.1,0.5,0.75)
sim.list=list()
totdata=0
data.num = c(0,0)

true.beta.A = c(-1,0.4,0.4,0.4,0.4)
true.beta = true.beta=c(3,-1,rep(-0.5,4))
true.sigma2 = 4


for( t in 1:length(seeds)){
  sim.list.all = simdata.bsa(lambda.a.list,lambda.s.list,true.beta.A,true.beta,true.sigma2,seeds=seeds[t])
  for( i in  1: sim.list.all$totdata){
    true.lambda = sim.list.all$data.num[i,]
    simdata = sim.list.all$sim.list[[i]]
    n =dim(simdata)[1]
    d = simdata$status
    #y = simdata$eventtime
    colnames(simdata)
    U =   as.matrix(simdata %>% dplyr::select(U) )
    X.s = as.matrix(simdata %>% dplyr::select(V1,A,C1,C2,C3,C4))
    C = as.matrix(simdata %>% dplyr::select(C1,C2,C3,C4))
    X.a = as.matrix(simdata %>% dplyr::select(V1,C1,C2,C3,C4)) 
    A = as.matrix(simdata %>% dplyr::select(A)) 
    id1 = which(A==1);nid1=length(id1)
    id0 = which(A==0);nid0=length(id0)
    X.a1 = X.a[id1,]
    X.a0 = X.a[id0,]
    
    head(X.s)
    p = dim(X.s)[2]
    p.a = dim(X.a)[2]
    
    censored.id <- which(d == 0)
    n.censored  <- length(censored.id)  # number of censored observations
    X.censored = cbind(X.s, U)[censored.id,]
    y.s <- logtime <- log(simdata$eventtime)   # for coding convenience, since the whole code is written with y
    #y <- scale(y)
    y.s.censored <- logtime[censored.id]
    
    #S <- 25000 ; burn = 5000; thin =100
    S <- 10000 ; burn = 5000; thin =25
    outlength = (S - burn)/thin
    ateout.reg          <- rep(0,outlength)
    ateout.ipw          <- rep(0,outlength)
    drout           <- rep(0,outlength)
    betasout        <- matrix(0,nrow=p,ncol=outlength)
    
    #---priors
    # prior for sigma2 - invGamma/ improper prior 1/\sigma2
    sigma2.a0    <- 0.1    #  
    sigma2.b0    <- 0.1    # 
    
    # prior for beta.s - MVN
    beta0     <- rep(0,p)
    Sigma20   <- diag(c(5^2, rep(2^2,(p-1))))
    
    # prior for beta.a - MVN
    beta0.a     <- rep(0,p.a)
    Sigma20.a   <- diag(c(5^2, rep(2^2,(p.a-1))))
    
    #---initials
    beta.s = rep(0,p)
    beta.a = rep(0,p.a)
    sigma2 = 3^2
    #lambda.s = 0
    #lambda.a = 0
    pi.0 = 0.5
    ranks = match(A, sort(unique(A)));uranks = sort(unique(ranks))
    KK =length(uranks)
    z = qnorm(rank(A, ties.method = "random")/(n+1))
    
    beta.s = true.beta[1:p]
    sigma2 = true.sigma2
    lambda.s = -0.5
    lambda.a = -0.5
    #lambda.s = true.lambda[2]
    #lambda.a = true.lambda[1]
    
    #set.seed(1258)
    # Gibbs sampling
    for(s in 1:S){
      y.s       <- sample.ycens.impute(n.censored, X.censored,y.s,sigma2,beta.s,lambda.s,p)
      U         <- sample.U.impute(n,y.s,A,C,beta.s,beta.a, lambda.s, lambda.a, pi.0, sigma2)
      beta.s    <- sample.beta.s(X.s,y.s,U,sigma2,lambda.s,p)
      lambda.s  <- sample.lambda.s(X.s,y.s,U,sigma2,beta.s,a=-3,b=3)
      sigma2    <- sample.sigma2(X.s,y.s,U,lambda.s,beta.s,p)
      beta.a    <- sample.beta.a(X.a,z,U,lambda.a,p.a)
      lambda.a  <- sample.lambda.a(X.a,z,U,beta.a,a=-3,b=3)
      z         <- sample.z(X.a,z,U,beta.a,lambda.a)
      ate.list  <- sample.DR(y.s, beta.s,lambda.s, beta.a,lambda.a,U,1,0)
      ate.reg   <- ate.list$ate.reg
      ate.ipw   <- ate.list$ate.ipw
      dr        <- ate.list$dr
      
      if (s>burn & s%%thin == 0)
      {
        betasout[,(s-burn)/thin]       <- beta.s
        #lambdasout[i]      <- lambda.s
        #betaaout[,i]       <- beta.a
        #lambdaaout[i]      <- lambda.a
        #sigma2out[i]       <- sigma2
        #zout[,i]           <- z
        #Uout[,i]           <- U
        ateout.reg[(s-burn)/thin]      <- ate.reg
        ateout.ipw[(s-burn)/thin]      <- ate.ipw
        drout[(s-burn)/thin]           <- dr
        #print(s)
      }
    }
    #print(i)
    ate.REG[,i] = c(mean(ateout.reg), median(ateout.reg), quantile(ateout.reg,probs = c(0.025,0.975)))
    ate.IPW[,i] = c(mean(ateout.ipw), median(ateout.ipw), quantile(ateout.ipw,probs = c(0.025,0.975)))
    ate.DR[,i]  = c(mean(drout), median(drout), quantile(drout,probs = c(0.025,0.975)))
  }
  print(t)
  ATE.REG[[t]] = ate.REG
  ATE.IPW[[t]] = ate.IPW
  ATE.DR[[t]]  = ate.DR
}
Sys.time()-start.time


## calculate coverage rate
CR.reg.bsa = CR.ipw.bsa = CR.dr.bsa = matrix(0,50,49)
CR.reg.bsa.point = CR.ipw.bsa.point = CR.dr.bsa.point = array(0,dim=c(2,50,49))
CR.reg.bsa.mean = CR.ipw.bsa.mean = CR.dr.bsa.mean = matrix(0,50,49)
for( i in  1: 50){
  for (iter in 1:49){
    CR.reg.bsa[i,iter] =  ifelse(ATE.REG[[i]][4,iter] >= -1 & -1 >= ATE.REG[[i]][3,iter],1,0)
    CR.reg.bsa.point[,i,iter] = ATE.REG[[i]][1:2,iter] 
    CR.reg.bsa.mean[i,iter] = ATE.REG[[i]][1,iter] 
  }
  
}

for( i in  1: 50){
  for (iter in 1:49){
    CR.ipw.bsa[i,iter] =  ifelse(ATE.IPW[[i]][4,iter] >= -1 & -1 >= ATE.IPW[[i]][3,iter],1,0)
    CR.ipw.bsa.point[,i,iter] = ATE.IPW[[i]][1:2,iter] 
    CR.ipw.bsa.mean[i,iter] = ATE.REG[[i]][1,iter] 
  }
}

for( i in  1: 50){
  for (iter in 1:49){
    CR.dr.bsa[i,iter] =  ifelse(ATE.DR[[i]][4,iter] >= -1 & -1 >= ATE.DR[[i]][3,iter],1,0) 
    CR.dr.bsa.point[,i,iter] = ATE.DR[[i]][1:2,iter] 
    CR.dr.bsa.mean[i,iter] = ATE.REG[[i]][1,iter] 
  }
}


CRpost.reg.bsa = apply(CR.reg.bsa, 2, sum)/50
CRpost.ipw.bsa = apply(CR.ipw.bsa, 2, sum)/50
CRpost.dr.bsa = apply(CR.dr.bsa, 2, sum)/50

CR.reg.bsa.repeat=apply(CR.reg.bsa.point,c(1,3),quantile,probs = c(0.025,0.975))
CR.ipw.bsa.repeat=apply(CR.ipw.bsa.point,c(1,3),quantile,probs = c(0.025,0.975))
CR.dr.bsa.repeat=apply(CR.dr.bsa.point,c(1,3),quantile,probs = c(0.025,0.975))





