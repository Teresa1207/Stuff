#Newton-Raphson Algorithm

#The Newton Rhapson function takes in the following arguments:
#1. Function to find the root (f)
#2. First Derivative of func, (deriv)
#3. Initial starting point, (start)
#4. Level of acceptance (Tolerance)
#5. Maximum number of iterations (max.iterations)
# If either the func or deriv take in a value that forces the function to divide by zero,
# we need the function to not evaluate and move on. 

#Newton-Raphson Algorithm

NewtRhap<-function(f,deriv,start,tolerance,max.iterations,debug){
  x=start
  count=1
  repeat{
    if ((debug=TRUE & is.infinite(f(x))) | (debug=TRUE & is.infinite(deriv(x)))){
      print("Diverged or a function is attempting to divide by zero. Try a new start value")
      break}else{
    x=x-f(x)/deriv(x)
    count=count+1
    if(abs(f(x))<tolerance || count==max.iterations){
      print(paste("MLE of ",expression(lambda),"=",x))
      print(paste("Number of Iterations = ",count))
      print(paste("Derivative of Likelihood evaluated f'(lambda)=",f(x)))
        break
      }
    }
  }
}




#Implement For the linkage problem

f.y<-function(r){125/(2+r) - 38/(1-r) + 34/r} #Derivative of the loglikelihood function
d.y<-function(r){-125/(2+r)^2 - 38/(1-r)^2 - 34/r^2} #Second Derivative of the loglikelihood function

NewtRhap(f.y,d.y,0,0.000001,100,TRUE)