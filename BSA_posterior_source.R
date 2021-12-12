
########################################################
##   example code for BSA for unmeasured confouding   ##
##   developed by: Shuang Wang                        ##
##   last updated: 12 Dec 2021                        ##
########################################################


################################
#### Gibbs Sampler functions ###
################################
## Update \pi(U=1)
post.pi.U = function(n,y.s, A, C, U, beta.s, lambda.s, beta.a, lambda.a,pi.0,sigma2){
  exp.y = cbind(rep(1,n), A, C, U) %*% as.matrix(c(beta.s, lambda.s),col=1)
  ll.y  = dnorm(y.s, exp.y, sqrt(sigma2))
  porm.a = pnorm(cbind(rep(1,n),  C, U) %*% as.matrix(c(beta.a, lambda.a),col=1))
  ll.a = porm.a^A*(1-porm.a)^(1-A)
  post.pi = pi.0*ll.a*ll.y
  return(post.pi)
}

## Impute survival time for censored
sample.ycens.impute = function(n.censored, X.censored,y.s,sigma2,beta.s,lambda.s,p){
  beta.imput    =  as.matrix(c(beta.s, lambda.s) ,col=1)  # used for imputation
  mean.impute   = X.censored %*% beta.imput
  sd.impute     = sqrt(sigma2) #sqrt(sigmaSq.imput)
  ## update censored data ##
  time.censored = msm::rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = y.s.censored) # truncated at log(time) for censored data
  y.s[censored.id] = time.censored
  return(y.s)
}

## Sample U from its posterior
sample.U.impute = function(n,y.s,A,C, beta.s, beta.a, lambda.s, lambda.a, pi.0, sigma2){
  U1 = rep(1,n); U0 = rep(0,n)
  probU1 = post.pi.U(n,y.s, A, C, U1, beta.s, lambda.s, beta.a, lambda.a,pi.0,sigma2)
  probU0 = post.pi.U(n,y.s, A, C, U0, beta.s, lambda.s, beta.a, lambda.a,pi.0,sigma2)
  post.pi = probU1/(probU1+probU0)
  U = rbinom(n,1,post.pi)
  return(U)
}

## Update beta coeff in outcome model - aft with lognormal distribution
sample.beta.s = function(X.s,y.s,U,sigma2,lambda.s,p){
  beta.s = rep(NA,p)
  #V.beta = n0/(n0+1)*solve(t(X)%*%X)
  V.beta = solve(solve(Sigma20)+t(X.s)%*%X.s/sigma2)
  #E.beta = V.beta%*%t(X)%*%z
  E.beta = V.beta%*%(solve(Sigma20)%*%beta0 + t(X.s)%*%(y.s - lambda.s*U))/sigma2
  beta.s = mvtnorm::rmvnorm(1,E.beta,V.beta)
  return(beta.s)
}

## Update bias parameter in outcome model
sample.lambda.s = function(X.s,y.s,U,sigma2,beta.s,a,b){
  lambda.s = 0 
  V.lambda.s = solve(t(U)%*%U/sigma2)
  E.lambda.s = V.lambda.s%*%(t(U)%*%(y.s -X.s%*%as.matrix(c(beta.s), col=1)) )/sigma2
  lambda.s   = msm::rtnorm(1, mean = E.lambda.s, sd = sqrt(V.lambda.s), lower =a,upper = b )
  return(lambda.s)
}

## Update sigma^2 in outcome model
sample.sigma2 = function(X.s,y.s,U,lambda.s,beta.s,p){
  beta.all = as.matrix(c(beta.s,lambda.s),col=1)
  nu = (n + p + 1) 
  ssn = sum((y.s - cbind(X.s,U)%*%beta.all)^2) 
  sigma2 = 1/rgamma(1, nu/2 + sigma2.a0, ssn/2 +sigma2.b0)
  return(sigma2)
}

## Update beta coeff in treatment model - probit model 
sample.beta.a = function(X.a,z,U,lambda.a,p.a){
  beta.a = rep(NA,p.a)
  #X.a = cbind(rep(1,n), C)
  V.beta.a = solve(solve(Sigma20.a)+t(X.a)%*%X.a)
  E.beta.a = V.beta.a%*%(solve(Sigma20.a)%*%beta0.a + t(X.a)%*%(z - lambda.a*U))
  beta.a = mvtnorm::rmvnorm(1,E.beta.a,V.beta.a)
  return(beta.a)
}

## Update beta coeff in treatment model - probit model - g-prior
sample.beta.a2 = function(X.a,z,U,lambda.a,p.a){
  beta.a = rep(NA,p.a)
  X.a = cbind(rep(1,n), C1)
  iXX<-solve(t(X.a)%*%X.a)  ; V<-iXX*(n/(n+1)) ; cholV<-chol(V)
  E.beta = V%*%( t(X.a)%*%(z - lambda.a*U) )
  beta.a = cholV%*%rnorm(p.a) + E.beta
  return(beta.a)
}

## Update bias parameter in treatment model
sample.lambda.a = function(X.a,z,U,beta.a,a,b){
  lambda.a = 0
  V.lambda.a = solve(t(U)%*%U)
  E.lambda.a = V.lambda.a%*%(t(U)%*%(z -X.a%*%as.matrix(c(beta.a),col=1)) ) 
  lambda.a   = msm::rtnorm(1, mean = E.lambda.a, sd = sqrt(V.lambda.a), lower = a, upper  = b)
  return(lambda.a)
}

## Update auxiliary variable in treatment model - apply data aug-memtation approach followed by Tanner and Wong (1987)
sample.z = function(X.a,z,U,beta.a,lambda.a){
  z=rep(0,n)
  mean.impute.b1 <- X.a1 %*% as.matrix(c(beta.a),col=1) + lambda.a*U[id1]
  mean.impute.b0 <- X.a0 %*% as.matrix(c(beta.a),col=1) + lambda.a*U[id0]
  z[id1]       <- msm::rtnorm(nid1, mean = mean.impute.b1, lower = 0)  # Equation (6) of Albert and Chib (1993)
  z[id0]       <- msm::rtnorm(nid0, mean = mean.impute.b0, upper = 0)
  return(z)
}


#################################
#### Causal effect estimators ###
#################################

## Doubly Robust Estimaters
sample.DR = function(y.s, beta.s,lambda.s, beta.a,lambda.a,U,a1,a0){
  pA = pnorm( X.a %*% as.matrix(c(beta.a),col=1) + lambda.a*U)
  E.exp = cbind(rep(1,n),rep(a1,n),C,U)%*%as.matrix(c(beta.s,lambda.s),col=1)
  E.unexp = cbind(rep(1,n),rep(a0,n),C,U)%*%as.matrix(c(beta.s,lambda.s),col=1)
  dr = mean( (y.s*A - (A-pA)*E.exp)/pA -  (y.s*(1-A) + (A-pA)*E.unexp)/(1-pA))
  ate.reg = mean(E.exp - E.unexp) 
  ate.ipw = mean(y.s*A/pA -y.s*(1-A)/(1-pA) )
  return(list(dr=dr,ate.reg=ate.reg,ate.ipw=ate.ipw))
} 

sample.DR2 = function(y.s, beta.s,lambda.s, beta.a,lambda.a,U,a1,a0){
  pA = pnorm( X.a %*% as.matrix(c(beta.a),col=1) + lambda.a*U)
  E.exp = cbind(rep(1,n),rep(a1,n),C,U)%*%as.matrix(c(beta.s,lambda.s),col=1)
  E.unexp = cbind(rep(1,n),rep(a0,n),C,U)%*%as.matrix(c(beta.s,lambda.s),col=1)
  dr = mean( y.s*A/pA - y.s*(1-A)/(1-pA) - (A-pA)*E.exp/pA - (A-pA)*E.unexp/(1-pA) )
  ate.reg = mean(E.exp - E.unexp) 
  ate.ipw = mean(y.s*A/pA -y.s*(1-A)/(1-pA) )
  return(list(dr=dr,ate.reg=ate.reg,ate.ipw=ate.ipw))
} 

sample.DR.naive = function(y.s,fit.s.coeff ,lambda.s,fit.a.coeff,lambda.a,U,a1=1,a0=0,data){
  d = data$status
  y.s = log(data$eventtime)
  U =   as.matrix(data %>% dplyr::select(U) )
  X.s = as.matrix(data %>% dplyr::select(V1,A,C1,C2,C3,C4))
  C = as.matrix(data %>% dplyr::select(C1,C2,C3,C4))
  X.a = as.matrix(data %>% dplyr::select(V1,C1,C2,C3,C4)) 
  A = as.matrix(data %>% dplyr::select(A)) 
  
  dr = dr1 =dr0= ate.reg = ate.ipw = NULL
  #new1 =as.data.frame(cbind(rep(1,n), rep(a1 ,n),C)) ; new0 = as.data.frame(cbind(rep(1,n), rep(a0,n),C))
  new1 =as.matrix(cbind(rep(1,n), rep(a1 ,n),C)) ; new0 = as.matrix(cbind(rep(1,n), rep(a0,n),C))
  #pA = pnorm(X.a %*% as.matrix(c(fit.a.coeff),col=1) + lambda.a*U)
  pA = pnorm(X.a %*% as.matrix(c(fit.a.coeff),col=1))
  E.exp = new1%*%as.matrix(c(fit.s.coeff),col=1)  #stats::predict(fit.s.coeff, new1, type="linear") 
  E.unexp = new0%*%as.matrix(c(fit.s.coeff),col=1)  #stats::predict(fit.s.coeff, new0, type="linear")   
  dr = mean((y.s*A - (A-pA)*E.exp)/pA -  (y.s*(1-A) + (A-pA)*E.unexp)/(1-pA))
  ate.reg =  mean(E.exp - E.unexp )
  ate.ipw =  mean(y.s*A/pA -y.s*(1-A)/(1-pA)  )
  return(list(dr=dr,ate.reg=ate.reg,ate.ipw=ate.ipw))
} 


#################################
#### data simulation function ###
#################################
simdata.bsa = function(lambda.a.list,lambda.s.list,beta.A,beta.Y,true.sigma2,seeds){
  for (s in 1:length(lambda.a.list)){
    for (j in 1:length(lambda.s.list)){
      data.num = rbind(data.num,c(lambda.a.list[s],lambda.s.list[j]))
      totdata = totdata+1
      # covariates
      set.seed(seeds)
      inter = rep(1,n)
      U = rbinom(n, 1L, 0.5)
      C1 = rnorm(n, 0, 1)
      C2 = rnorm(n, 0, 1)
      C3 = rbinom(n, 1, 0.5)
      C4 = rbinom(n, 1, 0.5)
      C = cbind(C1,C2,C3,C4)
      
      # trt (binary): A
      beta.A.all = c(beta.A,lambda.a.list[s])
      A  = rbinom(n, 1, pnorm(cbind(inter,C,U)%*%beta.A.all)) 
      #A = rbinom(n, 1L, 0.2)
      mean(A)
      
      # outcome (survival weibull dist.): Y
      #W = cbind( W1,W2)
      Xs<- cbind(rep(1,n), A, C, U)
      
      set.seed(1258) #9568
      beta.Y.all=c(beta.Y,lambda.s.list[j])
      #true.mu = t(rmnorm(1, X%*%true.beta, 0.1*diag(n)))
      true.mu= rep(0,n)
      for (i in 1:n){
        true.mu[i] = rnorm(1,Xs[i,]%*%beta.Y.all, sqrt(true.sigma2))
      }
      #mean(true.mu)
      #true.lambda = X.all%*%true.beta
      
      survt = exp(true.mu)
      #--add administrative censoring  
      tmax = quantile(survt,prob=0.3)
      d = tmax < survt # censoring indicator
      table(d==1)
      survt[d==1] = tmax 
      status = rep(1,n)
      status[d==1] = 0
      table(status) 
      s2 = data.frame(status, eventtime=survt, Xs)
      mean(s2$status)
      sim.list[[totdata]] = data.frame(s2)
      colnames(sim.list[[totdata]]) = colnames(s2)
    }
  }
 colnames(data.num) = c("lambda.a","lambda.s")
 data.num = data.num[-1,]
 return(list(sim.list=sim.list,data.num=data.num,totdata=totdata) )
}






