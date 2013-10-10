#Homework 00
#Teresa Filshtein
# Q1 Write a program that prints the numbers from 1 to 100. 
# But for multiples of three print "Fizz" instead of the number and 
# for the multiples of five print "Buzz". For numbers which are multiples 
# of both three and five print "FizzBuzz". 
A<-c(1:100)
A3<-A/3
A5<-A/5
A15<-A/15
FB<-vector()
for(i in 1:100)
  + {if(is.wholenumber(A15[i])){FB[i]="FizzBuzz"} 
     +  else{if(is.wholenumber(A5[i])){FB[i]="Buzz"}
             +  else{if(is.wholenumber(A3[i])){FB[i]="Fizz"} 
                     +  else{FB[i]=A[i]}}}}
#Q2
#Write a program that generates 10,000 uniform random numbers 
#between 0 and equation(x), and 10,000 
#uniform random numbers between 0 and 1 (y)
#You will then have 10,000 pairs of random numbers.
#Transform equation to equation where: equation, and, equation.
#Make a 2D scatterplot of the 10,000 (u,v) pairs. 
#What is the distribution of r=(sqrt(u^2+v^2))
x<-runif(10000,0,2*pi)
y<-runif(10000,0,1)
u<-y*cos(x)
v<-y*sin(x)
plot(u,v,main="2D Scatterplot of u vs v", xlab="u=y*cosx", ylab="v=y*sinx")

#Q3
#Consider the following snippet:
#"Hello, my name is Bob. I am a statistician. I like statistics very much."
#Part A
#Write a program to spit out every character in the snippet to a separate file

filename<-vector()
for(i in 1:nchar(sentence,type="chars")){filename[i]=paste("ch", i, ".txt", sep = "")}
for(i in 1:72){write(substring(sentence,i,i),filename[i])}
#Part B
#Write a program to combine everything back together!



