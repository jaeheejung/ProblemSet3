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
##1
OutofStep <- read.table("http://pages.wustl.edu/montgomery/incumbents.txt",header=TRUE,row.names="x")
#Reads in data from Jacob's website with the first row of the dataset as the names of the variables and the variable "x" as the numbers of observations
index <- sample(1:nrow(OutofStep),nrow(OutofStep)/2)
#Creates a vector of randomly sampled row numbers from OutofStep dataset (half)
training.data <- OutofStep[index,]
#This partition of OutofStep is my training set.
testing.data <- OutofStep[-index,]
#This partition of OutofStep is my testing set.

model.1 <- lm(voteshare~chalspend+incspend,data=training.data)
#model.1 has challenger spending and incumbent spending as explanatory variables
model.2 <- lm(voteshare~presvote+inparty+midterm,data=training.data)
#model.2 has presidential candidate vote share, the party of the incumbent as the president's party, and midterm election as explanatory variables
model.3 <- lm(voteshare~seniority+urban,data=training.data)
#model.3 has seniority of the incumbent and urban population as explanatory variables

predict.1 <- predict(model.1,newdata=testing.data)
#Predicts values of voteshare in testing.data using model.1
predict.2 <- predict(model.2,newdata=testing.data)
#Predicts values of voteshare in testing.data using model.2
predict.3 <- predict(model.3,newdata=testing.data)
#Predicts values of voteshare in testing.data using model.3

##2
P <- matrix(c(predict(model.1,newdata=testing.data), predict(model.2,newdata=testing.data), predict(model.3,newdata=testing.data)),ncol=3)
#P is a 3344x3 matrix of predictions about testing.data that I calculated with the three models I constructed above.

y <- testing.data$voteshare
#y is a vector of "true" observed voteshares separated out from testing.data

library(DataCombine)
#This package is needed to construct lagged variables.
data.with.lagged.var <- slide(testing.data,Var="voteshare",slideBy=-1)
#slide() creates lagged observations of voteshare.
r <- data.with.lagged.var$"voteshare-1"
#r is naive forecast, a vector of lagged observations of voteshare

fit.stats.matrix <- function(y,P,r){
#This function creates a matrix of six fit statistics for the three models I constructed before.
rmse <- function(P){
	sqrt(sum((P-y)^2,na.rm=TRUE)/length(y))
}
#rmse() calculates RMSE
rmse.vec <- aaply(.data=P,.margins=2,.fun=rmse)
#A vector of the values of RMSE for the three models is obtained by using aaply() across the columns of the prediction matrix.
mad <- function(P){
	median(abs(P-y),na.rm=TRUE)
}
#mad() calculates MAD
mad.vec <- aaply(.data=P,.margins=2,.fun=mad)
#A vector of the values of MAD for the three models is obtained by using aaply() across the columns of the prediction matrix.
rmsle <- function(P){
	sqrt(sum((log(P+1)-log(y+1))^2,na.rm=TRUE)/length(y))
}
#rmsle() calculates RMSLE
rmsle.vec <- aaply(.data=P,.margins=2,.fun=rmsle)
#A vector of the values of RMSLE for the three models is obtained by using aaply() across the columns of the prediction matrix.
mape <- function(P){
	sum((abs(P-y)/abs(y))*100,na.rm=TRUE)/length(y)
}
#mape() calculates MAPE
mape.vec <- aaply(.data=P,.margins=2,.fun=mape)
#A vector of the values of MAPE for the three models is obtained by using aaply() across the columns of the prediction matrix.
meape <- function(P){
	median((abs(P-y)/abs(y))*100,na.rm=TRUE)
}
#meape() calculates MEAPE
meape.vec <- aaply(.data=P,.margins=2,.fun=meape)
#A vector of the values of MEAPE for the three models is obtained by using aaply() across the columns of the prediction matrix.
mrae <- function(P){
	median(abs(P-y)/abs(r-y),na.rm=TRUE)
}
#mrae() calculates MRAE
mrae.vec <-aaply(.data=P,.margins=2,.fun=mrae)
#A vector of the values of MRAE for the three models is obtained by using aaply() across the columns of the prediction matrix.
return(cbind(rmse.vec,mad.vec,rmsle.vec,mape.vec,meape.vec,mrae.vec))
#The returning matrix is created by column binding the vectors of the six fit statistics.
}

#3





