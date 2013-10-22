#Bootstrapping Excersize

x0<-c(rep(1,100))
x1<-rnorm(100)
x2<-rnorm(100)
B<-c(1.2,0.3,-0.9)
X<-cbind(x0,x1,x2)
y<-X%*%B
Boot<-matrix(0,100,1000)
for(i in 1:1000){Boot[,i]=sample(100,size=100,replace=TRUE,prob=NULL)}
Boot_y<-matrix(0,100,1000)
for (i in 1:1000){for (j in 1:100){Boot_y[j,i]=y[Boot[j,i]]}}

B.hat<-matrix(0,3,1000)

for(i in 1:1000){B.hat[,i]=solve(t(X)%*%X)%*%t(X)%*%Boot_y[,i]}
S<-matrix(0,3,1)
for (i in 1:3){S[i]=sd(B.hat[i,])}