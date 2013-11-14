#======= Hive =======+
setwd("~/GitHub/Stuff/HW2")
#Load data from 'final' folder

data<-read.table('BigSummaryStats.txt', header=FALSE,sep="\001")

png('ScatterPlot.png')
plot(data$V2,data$V3,main='ScatterPlot: Group Means vs Group Variances',
     ylab='Within Group Variances',xlab='Within Group Means',col='blue')
dev.off()