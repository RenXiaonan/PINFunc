###############################################################
### Sample code for generating and estimating               ###
### simulated datasets.                                     ###
###############################################################

##################### Call packages and functions #################

library(MASS)
source("PINFunc.R")

######################### Homogeneity model ########################
### Generate simulated datasets:            

# There are THREE studies, set values for input variables

n1=200
n2=175
n3=150

p=50 # Each study has the same number of variables
corval=0.5 # Covariance of the covariate are AR(0.5) structure

b1=matrix(c(4,4,-8,-8,2,rep(0,45)),ncol = 1)
b2=matrix(c(-2,-2,-4,-4,1,rep(0,45)),ncol=1)
b3=matrix(c(-1.5,-1.5,3,-3,-0.75,rep(0,45)),ncol=1)

# Generate simulated datasets under the setting of Case II: Non-linear under homogeneity model.

# tuning parameters
m=3 # number of layers
eta=0.3 # learning rate
lambda=0.02 # tuning parameter for variable selection in the objective function
Lambda=0.0001 # tuning parameter for preventing overfitting in the objective function

# Each time 80% of the samples are randomly selected to fit the network and obtain important variables, repeated 10 times.
eva=NULL

s=1
for (s in 1:10) {

  set.seed(s)
  data=DataHoNonLinear(n1, n2, n3, b1, b2, b3, p, corval)
  trainrate=0.8
  dataDiv=DataDiv(data, trainrate)
  dataClean=DataC(dataDiv)
  
  x=dataClean[[1]]
  y=dataClean[[2]]
  
  set.seed(1)
  ho_pin=hopin(x,y,m,varepsilon=1e-3,maxiter=500,alpha=0.1,eta,lambda,Lambda,a=0.04,threshold=0.15,node=round(sqrt(p)))
  index=ho_pin$index
  eva=rbind(eva,Eva(index,b1))
  
  s=s+1
  
}

ho_result = matrix(apply(eva,2,mean),nrow = 1)
colnames(ho_result) = c("SEN","SPE","GM","CCR")
ho_result

######################### Heterogeneity model ########################
### Generate simulated datasets:            

# There are THREE studies, set values for input variables
n1=200
n2=175
n3=150

p=50 # Each study has the same number of variables
corval=0.5 # Covariance of the covariate are AR(0.5) structure

b1=matrix(c(4,4,-8,-8,2,-2,0,0,0,0,rep(0,40)),ncol = 1)
b2=matrix(c(-2,-2,4,-4,1,0,1,0,0,0,rep(0,40)),ncol=1)
b3=matrix(c(-1.5,-1.5,3,-3,-0.75,0,0,0.75,0,0,rep(0,40)),ncol=1)

# Generate simulated datasets under the setting of Case IV: Non-linear under heterogeneity model.

# tuning parameters
m=3 # number of layers
eta=0.5 # learning rate
lambda=0.03 # tuning parameter for variable selection in the objective function
beta=0.05 # tuning parameter for variable selection in the objective function
Lambda=0.001 # tuning parameter for preventing overfitting in the objective function

# Each time 80% of the samples are randomly selected to fit the network and obtain important variables, repeated 10 times.
eva1=NULL
eva2=NULL
eva3=NULL

s=1
for (s in 1:10) {
  
  set.seed(s)
  data=DataHeNonLinear(n1, n2, n3, b1, b2, b3, p, corval)
  trainrate=0.8
  dataDiv=DataDiv(data, trainrate)
  dataClean=DataC(dataDiv)
  
  x=dataClean[[1]]
  y=dataClean[[2]]
  
  set.seed(2)
  he_pin=hepin(x,y,m,varepsilon=1e-3,maxiter=500,alpha=0.1,eta,lambda,beta,Lambda,a=0.04,sig=0.04,threshold=0.05,node=round(sqrt(p)))
  index1=he_pin$index[[1]]
  index2=he_pin$index[[2]]
  index3=he_pin$index[[3]]
  
  eva1=rbind(eva1,Eva(index1,b1))
  eva2=rbind(eva2,Eva(index2,b2))
  eva3=rbind(eva3,Eva(index3,b3))
  
  s=s+1
  
}

he_1 = matrix(apply(eva1,2,mean),nrow = 1)
he_2 = matrix(apply(eva2,2,mean),nrow = 1)
he_3 = matrix(apply(eva3,2,mean),nrow = 1)
he_result = rbind(he_1,he_2,he_3)
he_result = rbind(he_result,matrix(apply(he_result,2,mean),nrow = 1))
colnames(he_result) = c("SEN","SPE","GM","CCR")
rownames(he_result) = c("data1","data2","data3","mean")
he_result

