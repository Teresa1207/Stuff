
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(mvtnorm)
library(boot)
library(MCMCpack)
library(coda)


########################################################################################
########################################################################################

args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################

#install.packages("mvtnorm"), use for generating Bivariate Normal rv.
#install.packages("boot"), use for inverse logit function
#install.packages("emdbook"), use for multivariate normal density
#install.packages("MCMCpack"), MCMC package 

library(mvtnorm)
library(boot)
#library(emdbook)
library(MCMCpack)
library(coda)

#setwd("~/GitHub/Stuff/HW1/BayesLogit/data")

#load data
data<-read.csv(file=paste("data/blr_data_",sim_num,".csv",sep=""))
beta<-read.csv(file=paste("data/blr_pars_",sim_num,".csv",sep=""))

#Set Initial Matrices

##Generate Proposal Beta values
##Bivariate(u=Beta at time t, v^2)

##Define data, variables

mu<-c(0,0) #mean of the prior distribution of Beta 
s<-diag(2) #covariance matrix of prior on Beta
y<-as.matrix(data$y)
m<-as.matrix(data$n)
X<-as.matrix(cbind(data$X1,data$X2)) #design matrix
beta<-as.matrix(beta) #true beta values



##functions##
log.f.y<-function(B){dbinom(y,m,inv.logit(X%*%B),log=TRUE)} #this is f(y|B), log scale
log.f.B<-function(B){dmvnorm(B,mu,s,log=TRUE)} #this is p(B), prior dist. of B, log scale
log.target.density<-function(B){log.f.B(B)+sum(log.f.y(t(B)))} #this is the posterior distr. (proportionate)
alpha<-function(A,B){log.target.density(A)-log.target.density(B)}


#Generate Proposed Betas
#Proposed Density is bivariate normal, mean(beta(t),covariance matrix)

#I am thinking a good place to start for the covariance matrix is the covariance matrix
#of beta under linear regression. But also, another option is the covariance matrix of
#coefficients under logistic regression (thanks Chris Aden). I computed both.
#Compute LSE of Covariance Matrix of B
y.lin<-lm(y~data$X2)
MSE<-sum((y.lin$residuals)^2)/y.lin$df.residual


#Compute the covariance variance matrix of the Logistic Regression estimates

data["y.c"]<-NA
data$y.c<-data$n-data$y


y.log<-glm(formula=cbind(y,y.c)~X2, family="binomial", data=data)
log.beta<-matrix(c(y.log$coef[1],y.log$coef[2]))

data["pi"]<-NA
data$pi<-inv.logit(X%*%log.beta)

data["vi"]<-NA
data$vi<-data$n*data$pi*(1-data$pi)

vi<-as.vector(data$vi)
V<-diag(vi)

#sigma values:

#sigma<-MSE*solve(t(X)%*%X)
sigma<-solve(t(X)%*%V%*%X)
#sigma<-matrix(c(.5,.5,.5,.5),2,2)



#Generate proposed values of Beta
# USE THIS IN ALGORITHM AT EACH step t
#prop.beta<-rmvnorm(1,mean = beta.samples[t,],sigma=sigma,method="svd")

###METROPOLIS-ALGORTHIM!##


#Create Sample Beta Matrix 
burnin<-1000
niter<-15000
b1=0 #intial B1 at t=0
b2=0 #intial B2 at t=0
beta.samples<-matrix(NA,niter+burnin,2)
n.accept<-0
n.accept.burn<-0
beta.current<-matrix(c(b1,b2),1,2)
log.U<-log(runif(burnin+niter))
##Tuning Parameter for Sigma
v.tune<-1.75

#MH

bayes.log.reg<-for (t in 1:(niter+burnin)){
  prop.beta<-rmvnorm(1,mean = beta.current,sigma=sigma,method="svd") #propose a new beta
  log.alpha<-alpha(prop.beta,beta.current) #calculate alpha (accept or reject ratio)
  U<-log.U[t]
  if(U < log.alpha){beta.current<-prop.beta}
  if(t>burnin & U < log.alpha)
  {n.accept<-n.accept+1} 
  beta.samples[t,]<-beta.current
  if(t<(burnin+1) & U < log.alpha)
  {n.accept.burn<-n.accept.burn+1}
  if(t==(.25)*burnin|t==(.5)*burnin|t==(.75)*burnin|t==(1.0)*burnin){
    if(n.accept.burn/t> .40){sigma<-v.tune*sigma}
    else{if(n.accept.burn/t<.20)
    {sigma<-sigma/v.tune}}
    print(cat(paste("Burnin Acceptance rate is  ",100*round(n.accept.burn/t,2),"%",sep="")))
    flush.console()}
  if(t>burnin & t%%1000==0){print(cat(paste("Acceptance rate is ",100*round(n.accept/(t-burnin),2),"%",sep="")))
                            flush.console()
  }
  
}

#Discard the burnin beta.samples

beta.samples<-beta.samples[(burnin +1):(burnin+niter),]

#Get information about your run

cat(paste("Acceptance rate was ",100*round(n.accept/niter,2),"%\n",sep="")) 
acf(beta.samples)
ess <- effectiveSize(beta.samples)
cat("Effective sample size:\n") ; print(ess)

beta.samples<-mcmc(beta.samples)
pdf(file=paste("results/Plots",sim_num,".pdf",sep=""))
plot(beta.samples)
dev.off()
pv <- c(1:99)/100
beta.q <- apply(beta.samples,2,quantile,probs=pv)
colnames(beta.q)<-c("B1","B2")
write.table(beta.q,file = paste("results/beta.q.",sim_num,".csv",sep=""))






