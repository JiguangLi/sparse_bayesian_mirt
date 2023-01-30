# ***** Factor regression simulated example ****** #


source("FACTOR_CODE_update.R")

library(mvtnorm)
library(partitions)
library(nloptr)
library(glmnet)

n<-100	# Number of observations
G<-1956	# Number of responses 	
p<-10
dd<-5   # True number of factors

# TRUE LOADING MATRIX


B<-matrix(0,G,dd)
length=500
offset=135
end<-1
for(i in (1:dd)){
	start<-end-(i!=1)*offset
	end<-start+length-1
	B[start:end,i]<-1
}


myImagePlot(B,T)

# TRUE COVARIANCE MATRIX

Sigma<-B%*%t(B)+diag(G)

# GENERATE DATA

Y<-rmvnorm(n,numeric(G),Sigma)

# INITIALIZATIONS

K<-20

startB<-matrix(rnorm(G*K),G,K)

alpha<-1/G

lambda1<-0.001

epsilon<-0.05

myImagePlot(abs(startB),F)

title("Initialization")


# PXEM: Dynamic Posterior Exploration (Approximate M-step)

start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))

lambda0<-5
result_5<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-10
result_10<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-20
result_20<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-30
result_30<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100,TRUE)




# PXEM: Dynamic Posterior Exploration (Exact M-step): a bit slower

start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))

lambda0<-5
result_5_2<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,FALSE,100,TRUE)

lambda0<-10
result_10_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5_2,K,epsilon,alpha,TRUE,FALSE,100,TRUE)

lambda0<-20
result_20_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10_2,K,epsilon,alpha,TRUE,FALSE,100,TRUE)

lambda0<-30
result_30_2<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20_2,K,epsilon,alpha,TRUE,FALSE,100,TRUE)

