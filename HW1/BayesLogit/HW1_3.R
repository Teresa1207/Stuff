#Homework 3
#setwd("~/GitHub/Stuff/HW1/BayesLogit")

library(mvtnorm)
library(coda)
library(boot)
library(MCMCpack)
library(xtable)


##Load and Format the Data

data<-read.table("breast_cancer.txt", stringsAsFactors=FALSE)
#col.names<-c("area","compactness","concavepts","concavity","fracdim", "perimeter", "radius", "smoothness", "symmetry", "texture", "diagnosis")
col.names<-c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","diagnosis")
colnames(data)<-col.names
data<-data[2:570,]
data["D"]<-NA
data$D<-ifelse(data$diagnosis=="M",1,0) # 1 if M, 0 if B

      
y<-matrix(data$D,569,1)


## I realize I could have done this more efficient! Time time time!
F1<-scale(as.numeric(as.character(data$F1)))
F2<-scale(as.numeric(as.character(data$F2)))
F3<-scale(as.numeric(as.character(data$F3)))
F4<-scale(as.numeric(as.character(data$F4)))
F5<-scale(as.numeric(as.character(data$F5)))
F6<-scale(as.numeric(as.character(data$F6)))
F7<-scale(as.numeric(as.character(data$F7)))
F8<-scale(as.numeric(as.character(data$F8)))
F9<-scale(as.numeric(as.character(data$F9)))
F10<-scale(as.numeric(data$F10))
scaled.X<-cbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,data$D)
colnames(scaled.X)<-col.names<-c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","diagnosis")
scaled.X<-data.frame(scaled.X)
C<-c(rep(1,569))
X<-as.matrix(cbind(C,scaled.X[,1:10])) #standardize covariates

#Setup prior distribution of Beta#
mu<-c(rep(0,11))
s.var<-c(rep(1000,11))
s<-diag(s.var)

##functions##
log.f.y<-function(B){dbinom(y,1,inv.logit(X%*%B),log=TRUE)} #this is f(y|B), log scale
log.f.B<-function(B){dmvnorm(B,mu,s,log=TRUE)} #this is p(B), prior dist. of B, log scale
log.target.density<-function(B){log.f.B(B)+sum(log.f.y(t(B)))} #this is the posterior distr. (proportionate)
alpha<-function(A,B){log.target.density(A)-log.target.density(B)}


#Generate Proposed Betas
#Proposed Density is multivariate normal, mean(beta(t),covariance matrix)

#A good place to start for the covariance matrix is the covariance matrix
#of the coefficients under logistic regression (thanks Chris Aden).

#Compute the covariance variance matrix of the Logistic Regression estimates

y.log<-glm(formula=diagnosis~F1+F2+F3+F4+F5+F6+F7+F8+F9+F10,family="binomial", data=scaled.X, control = list(maxit = 50))

log.beta<-as.matrix(y.log$coef)

scaled.X["pi"]<-NA
scaled.X$pi<-inv.logit(X%*%log.beta)

scaled.X["vi"]<-NA
scaled.X$vi<-scaled.X$pi*(1-scaled.X$pi)

vi<-as.vector(scaled.X$vi)
V<-diag(vi)

#sigma values:
sigma<-solve(t(X)%*%V%*%X)



#Create Sample Beta Matrix 
burnin<-1000
niter<-50000
beta.samples<-matrix(NA,niter+burnin,11)
n.accept<-0
n.accept.burn<-0
beta.current<-matrix(c(rep(0,11)),1,11)
log.U<-log(runif(burnin+niter))
##Tuning Parameter for Sigma
v.tune<-1.75

#MH algorithm

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
beta_check<-beta.samples
#Get information about your run

cat(paste("Acceptance rate was ",100*round(n.accept/niter,2),"%\n",sep="")) 
beta.samples<-mcmc(beta.samples)
ess <- effectiveSize(beta.samples)
cat("Effective sample size:\n") ; print(ess)

beta.samples<-mcmc(beta.samples)
pdf("Plots_BC.pdf")
plot(beta.samples)
dev.off()

##Generate quantiles for beta.samples
pv1 <- c(1,5,25,50,75,95,99,100)/100
beta.q <- apply(beta.samples,2,quantile,probs=pv1)
beta.q<-round(beta.q,2)

#colnames
col.names1<-c("Intercept","area","compactness","concavepts","concavity","fracdim", "perimeter", "radius", "smoothness", "symmetry", "texture")
colnames(beta.q)<-col.names1

#Decide if covariate is related or not
sig.beta<-matrix(0,1,11)
for (i in 1:11){if(beta.q[1,i]<0 && beta.q[7,i]>0){sig.beta[,i]="NR"}else{sig.beta[,i]="R"}}
beta.q<-rbind(beta.q,sig.beta)


##posterior predictive check##
y_check<-matrix(0,569,niter)
for(i in 1:niter){
b<-as.matrix(beta_check[i,])
y.curr<-as.matrix(rbinom(569,1,inv.logit(X%*%b)))
y_check[,i]<-y.curr
}

png("sample.png")
sample_means<-apply(y_check,2,mean)
true_mean<-mean(y)
hist(sample_means, main="Histogram of Sample Means")
abline(v=true_mean,col="red")
dev.off()


###auto lag 1

lag1 = sapply(1:11, function (i) acf(beta.samples[,i],plot=F)$acf[2])

t = t(as.table(lag1))
colnamest<-c("\beta1","\beta2","\beta3","\beta4","\beta5","\beta6","\beta7","\beta8","\beta9","\beta10","\beta11")
colnames(t)<-colnamest

