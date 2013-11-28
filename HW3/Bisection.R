#==========Homework 3==========#

#Problem 1
#Bisection Algorithm

#The Bisection Algorithm function takes in the following arguments:
#1. Function to find the root (f)
#2. Lower bound, (lower)
#3. Upper bound, (upper)
#4. Level of acceptance (Tolerance)
#5. Maximum number of iterations (max.iterations)
# 
#setwd("~/Davis2013-2014/Fall2013/S250")

#Bisection Algorithm

Bisection<-function(f,lower,upper,tolerance,max.iterations,debug){
l=lower
u=upper
count=1
c=(l+u)/2
if ((debug=TRUE & is.finite(f(lower)*f(upper))&f(lower)*f(upper)>0)){
  print("Poor choice of Upper and Lower bound, f(upper)*f(lower)>0, no root between. Choose new bounds")
  }else{
   repeat{
  if ((debug=TRUE & is.infinite(f(c)))){
    print("Diverged or a function is attempting to divide by zero. Try different Upper/Lower bounds")
    break}else{
  if(abs(f(c))<tolerance || count == max.iterations){
    print(paste("MLE of ",expression(lambda),"=",c))
    print(paste("Number of Iterations = ",count))
    print(paste("Derivative of Likelihood evaluated f'(lambda)=",f(c)))
    break
  }else{
  count=count+1
  l<-ifelse(f(c)*f(l)<0,l,c)
  u<-ifelse(f(c)*f(l)<0,c,u)
  c=(l+u)/2
  
  }}}}}


#Implement For the linkage problemf.x =derivative of loglikelihood function
l.x = function(x){(2+x)^125*(1-x)^38*x^34}
f.x = function(x){
  125/(2+x) - 38/(1-x) + 34/x
}


lambda<-function(x){1-2*x+x^2}
png("link.png")
par(mfrow=c(2,1))
plot(l.x,main="Classic Linkage Likelihood Function", ylab="Likelihood",xlab=expression(lambda),xlim=c(0,1))
plot(f.x,main="Classic Linkage Derivative of log-Likelihood Function", ylab="Derivate of log-Likelihood",xlab=expression(lambda),xlim=c(-4,4),ylim=c(-1200,1200))
abline(h=0,v=c(-2,1,0),col=c("red","blue","blue","blue"),lty=c(1,2,2,2))
dev.off()


Bisection(f.x,0,4,0.000001,10000,TRUE)


