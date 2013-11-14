##Bag of Little BootStraps
# Fit a simple linear regression with 1,000 covariates (no intercept term)
# Goal is to find SE(beta.hat.1),...,SE(beta.hat.1000)


# Clear out everything in memory for a nice clean run and easier debugging
rm(list=ls())

s = 5           # s = number of distinct subsamples 
r = 50          # r = number of bootstrap samples to run for each subsample

###===Setup for running on Gauss===###

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  
  
  
  
  # Get the job number from the argument 1:(s*r) (i.e. 1:5*50 = 1:250)
  job = as.numeric(args[1])
  
  ##Need to get the s and r index
  # Get the s_index by using mod s.  
  # Also, if job mod s == 0, then it is subsample dataset s (i.e. 5)
  s_index = job %% s
  if (s_index == 0){
    s_index = s
  }
  
  # Get r_index 
  # Also, if job mod r == 0, then it is bootstrap sample r (i.e. 50) within subsample dataset s
  r_index = ceiling(job / s)
  if (r_index == 0){
    r_index = r
  }
  
  # The seed must be a function of the s_index to ensure that the subsample is the same
  # for same values of s (Thanks MBissel for your help with this)
  sim_seed <- (762*(s_index) + 121231)
  set.seed(sim_seed)
}


cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))





##++++++++++++ Bag of Little BootStraps R-coding ++++++++++++##


library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

##setwd("~/GitHub/Stuff/HW2")

#Set path names and file names for use on Gauss
datapath <- "/home/pdbaines/data/"
outpath <- "output/"
rootfilename <- "blb_lin_reg_data"
datafile = paste0(datapath, rootfilename, ".desc")

# Attach data and set parameters for BOLBS

data<-attach.big.matrix(datafile)
n=dim(data)[1]
g=0.7
b=ceiling(n^(g))
d=dim(data)[2]-1


# Choose a random subset of size b from the original data
# First pick the rows, then create the sample.
rows.subsample = sample(1:n, size=b, replace=FALSE)
X.subsample = data[rows.subsample,1:d]
Y.subsample = data[rows.subsample,d+1]

# Need to reset the simulation seed, now we will set the seed for the r bootstrap samples
# The seed for the bootstrap sample must be different for each bootstrap sample r within
# each subsample s. This is a function of s_index and r_index
sim_seed <- (762*(s_index) + 121231 + r_index)
set.seed(sim_seed)

# Sampling from the subsample. Sample n=10,0000 from your subsample.
# Use the multinomial dist. to assign a number of repeats to each of the 
# rows in the subsample.
boot.weights<-rmultinom(1, size = n, prob=rep(1/b, b))

# Fit the linear regression and get the coefficients. Apply the boot.weights.
# Make sure to account for no intercept term.
model = lm(Y.subsample ~ 0 + X.subsample, weights=boot.weights)
beta.hat = model$coefficients

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:
write.table(x=beta.hat,file=outfile, sep=",", col.names=TRUE, row.names=FALSE)

cat("done. :)\n")

q("no")

