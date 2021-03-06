########################################################
##   example code for BSA for unmeasured confouding   ##
##   developed by: Shuang Wang                        ##
##   last updated: 12 Dec 2021                        ##
########################################################

##########################
####  **** NOTE ****  ####
##########################

# assuming in the absense of unmeasured confouding
# Regression and bootstap

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

set.seed(1258)
#seeds = sample(seq(1,500),100)
#seeds =  seeds[1:50]
seeds = sample(seq(1,5000),500)
seeds = seeds[1:50]
ate.REG.Boot = ate.IPW.Boot = ate.DR.Boot = matrix(0,nrow=6,ncol=49)
ATE.REG.Boot = ATE.IPW.Boot = ATE.DR.Boot = list()

for( t in 1:length(seeds)){
  sim.list.all = simdata.bsa(lambda.a.list,lambda.s.list,true.beta.A,true.beta,true.sigma2,seeds=seeds[t])
  for( iter in  1: sim.list.all$totdata){
    true.lambda = sim.list.all$data.num[iter,]
    simdata = sim.list.all$sim.list[[iter]]
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
    y.s.censored <- logtime[censored.id]
    
    fit.surv = survreg(data = simdata, Surv(eventtime,d) ~ A+C1+C2+C3+C4 , dist = "lognormal" )
    fit.prob = glm(data = simdata, A ~ C1+C2+C3+C4 , family = binomial(link = "probit"))
    
    DR.naive.list = sample.DR.naive2(y.s, fit.surv$coefficients, 0,fit.prob$coefficients, 0,U,1,0,simdata)
    
    #fit.a.coeff = fit.prob$coefficients
    #pnorm(X.a %*% as.matrix(c(fit.a.coeff),col=1))[1:10]
    #predict(fit.prob, type = "response")[1:10]
    #REG.naive  = mean(DR.naive.list$ate.reg)
    #IPW.naive  = mean(DR.naive.list$ate.ipw)
    #DR.naive   =  mean(DR.naive.list$dr) 
    
    B <- 200
    DR.Boot = REG.Boot = IPW.Boot = rep(0,B)
    for (b in 1:B) {
      index <- sample(n, size = n, replace = T)
      newdf <- simdata[index, ]
      fit.surv.boot = survreg(data = newdf, Surv(eventtime,status) ~ A+C1+C2+C3+C4 , dist = "lognormal" )
      fit.prob.boot = glm(data = newdf, A ~ C1+C2+C3+C4 , family = binomial(link = "probit"))
      DR.boot.list = sample.DR.naive(y.s, fit.surv.boot$coefficients, 0, fit.prob.boot$coefficients, 0,U,1,0, newdf)
      DR.Boot[b]  = mean(DR.boot.list$dr)
      REG.Boot[b] = mean(DR.boot.list$ate.reg)
      IPW.Boot[b] = mean(DR.boot.list$ate.ipw)
    }
    ate.REG.Boot[,iter] = c(mean(DR.naive.list$ate.reg), median(DR.naive.list$ate.reg), mean(REG.Boot),median(REG.Boot), quantile(REG.Boot,probs = c(0.025,0.975)))
    ate.IPW.Boot[,iter] = c(mean(DR.naive.list$ate.ipw), median(DR.naive.list$ate.ipw), mean(IPW.Boot),median(IPW.Boot), quantile(IPW.Boot,probs = c(0.025,0.975)))
    ate.DR.Boot[,iter]  = c(mean(DR.naive.list$dr), median(DR.naive.list$dr), mean(DR.Boot),median(DR.Boot), quantile(DR.Boot,probs = c(0.025,0.975)))
  }
  print(t)
  ATE.REG.Boot[[t]] = ate.REG.Boot
  ATE.IPW.Boot[[t]] = ate.IPW.Boot
  ATE.DR.Boot[[t]]  = ate.DR.Boot
}


#--------calculate coverage rate
CR.reg.boot = CR.ipw.boot = CR.dr.boot = matrix(0,50,49)
CR.reg.boot.point = CR.ipw.boot.point = CR.dr.boot.point = array(0,dim=c(4,50,49))
CR.reg.boot.mean = CR.ipw.boot.mean = CR.dr.boot.mean = matrix(0,50,49)
for( i in  1: 50){
  for (iter in 1:49){
    CR.reg.boot[i,iter] =  ifelse(ATE.REG.Boot[[i]][6,iter] >= -1 & -1 >= ATE.REG.Boot[[i]][5,iter],1,0) 
    CR.reg.boot.point[,i,iter] = ATE.REG.Boot[[i]][1:4,iter] 
    CR.reg.boot.mean[i,iter] =ATE.REG.Boot[[i]][1,iter] 
  }
  
}

for( i in  1: 50){
  for (iter in 1:49){
    CR.ipw.boot[i,iter] =  ifelse(ATE.IPW.Boot[[i]][6,iter] >= -1 & -1 >= ATE.IPW.Boot[[i]][5,iter],1,0) 
    CR.ipw.boot.point[,i,iter] = ATE.IPW.Boot[[i]][1:4,iter] 
    CR.ipw.boot.mean[i,iter] =ATE.IPW.Boot[[i]][1 ,iter] 
  }

}

for( i in  1: 50){
  for (iter in 1:49){
    CR.dr.boot[i,iter] =  ifelse(ATE.DR.Boot[[i]][6,iter] >= -1 & -1 >= ATE.DR.Boot[[i]][5,iter],1,0) 
    CR.dr.boot.point[,i,iter] = ATE.DR.Boot[[i]][1:4,iter] 
    CR.dr.boot.mean[i,iter] =ATE.DR.Boot[[i]][1,iter] 
  }
  
}


CRpost.reg.boot = apply(CR.reg.boot, 2, sum)/50
CRpost.ipw.boot = apply(CR.ipw.boot, 2, sum)/50
CRpost.dr.boot  = apply(CR.dr.boot, 2, sum)/50
CR.reg.boot.repeat=apply(CR.reg.boot.point,c(1,3),quantile,probs = c(0.025,0.975))
CR.ipw.boot.repeat=apply(CR.ipw.boot.point,c(1,3),quantile,probs = c(0.025,0.975))
CR.dr.boot.repeat=apply(CR.dr.boot.point,c(1,3),quantile,probs = c(0.025,0.975))

