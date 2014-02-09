#Problem Set 3
#Jae Hee Jung


library(plyr)
library(foreach)
library(doMC)
library(multicore)

####SAMPLING DISTRIBUTIONS AND P-VALUES####
##1
random.data <- rnorm(10000) 
#A random data of size 10000 drawn from normal distribution
array.3d <- array(random.data,c(20,5,1000))
#Rearranges the data to an array of 20x5x1000

##2
Beta <- matrix(c(1,2,0,4,0),ncol=1)
#Beta is information from Jacob
Y.generator <- function(X){
	(X%*%Beta) + rnorm(20)
}
#Y.generator is a function that creates Y values 
Ys <- apply(X=array.3d, MARGIN=3, FUN=Y.generator)
#Applies the Y.generator function over the third dimensions of data array.3d, creating a 20x1000 array

##3
repeated.regression <- function(i,Y,X){
	lm(Y[,i]~X[,,i])[[1]]
}
#repeated.regression() gives the coefficients of the linear regression model of Y[,i] regressed on X[,,i]
coefficients <- laply(.data=1:1000, .fun=repeated.regression, Ys, array.3d)
#coefficients is a 1000x6 array, obtained by running repeated.regression() 1000 times with the corresponding arguments

##4
density.plots <- function(i){
	plot(density(coefficients[,i]))
}
#density.plots() draws a density plot of the coefficient in the ith column of the array coefficients
par(mfrow=c(3,2))
#This divides the plotting screen into 3 by 2
density.plots(1)
density.plots(2)
density.plots(3)
density.plots(4)
density.plots(5)
density.plots(6)
#These plots represent the sampling distributions of the parameters. We can see that the means of the coefficients approximate the covariates that Jacob provided us: Beta <- matrix(c(1,2,0,4,0),ncol=1)

##5
t.statistics <- function(i,Y,X){
	as.list(coef(summary(lm(Y[,i]~X[,,i])))[,"t value"])
}
#t.statistics() extracts t values from the coefficients table included in the summary of the linear regression model. I turn the t values into a list in order to use the laply function below.
all.t.statistics <- laply(.data=1:1000, .fun=t.statistics, Ys, array.3d)
#all.t.statistics is a 1000 by 6 array that contains the t statistics for all six coefficients in the 1000 regressions.

##6
critical.t <- qt(1-(0.05/2),df=14)
#This gives the value of the t statistic with 14 degrees of freedom at 0.05 significance level. The value is about 2.145. 14 degrees of freedom comes from 20-5-1.
significant.ts <- function(i){
	length(all.t.statistics[,i][all.t.statistics[,i] > critical.t | all.t.statistics[,i] < -(critical.t)])
}
#significant.ts() calculates the length of the ith column of all.t.statistics whose absolute values are larger than critical.t
laply(.data=1:6,.fun=significant.ts)
#This laply function applies significant.ts() from the first to the sixth column of all.t.statistics, generating a vector of the number of statistically significant t-statistics

##7
system.time(laply(.data=1:6,.fun=significant.ts))
#Measures system time without parallel. The results differ every time this function is run, but the changes are minimal. In general, the results are approximately user:0.007, system:0.001, elapsed:0.007.
registerDoMC(cores=4)
#Registers parallel backend
system.time(laply(.data=1:6,.fun=significant.ts,.parallel=TRUE))
#In parallel, system time is significantly larger for all three time measures. An instance is user:0.041, system:0.080, elapsed:0.064. I tried running system.time with different numbers of cores, but system time only gets larger as the number of cores increases. In essence, no time is saved.

####CALCULATING FIT STATISTICS####





