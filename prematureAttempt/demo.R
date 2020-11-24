setwd("C:/Users/Peter/Desktop/phd/R/demo")
samples <- 100

#'
#' # Simulation case study
#'
library(VineCopula)
set.seed(2019)
#' let U1, U2 be uniformly distributed random variables created from V1 and V2 by probability integral transformation. We will define bivariate distribution of U1 and U2 by bi-copula from student copula family with 6 degrees of freedom and strong tail dependency
copula_obj <- BiCop(family = 2, par = 0.8, par2 = 6)
copula_obj

copula_samples <- BiCopSim(samples, copula_obj)
plot(copula_samples)
U1 <- copula_samples[,1]
U2 <- copula_samples[,2]
#' for each value of sample for U1 calculate V1 by inverse normal transformation from normal marginal distribution N(0,1)
V1 <- qnorm(U1)
#' we see that V1 is normally distributed
hist(V1)
var(V1)
#' V1 is defined as sum of two random variables X1 and X2. we will make simplyfing assumption that X1 and X2 are normally distribtued. We will generate samples of X1 and X2 such as: (1) first generate samples of X1
par1 <- 0.7
par2 <- 0.3
X1 <- vector(length = samples)
X2 <- vector(length = samples)
for(i in 1:samples){
  X1[i] <- par1 * V1[i]
  X2[i] <- par2 * V1[i]
}
#' Let's look at the situation of second cluster. Here we assume that value of both random variables is driven by value of aggregation variable. For simplicity we will aslo assume that aggregation variable is normally distributed N(0,1)
V2 <- qnorm(U2)
hist(V2)
var(V2)
#' for purpose of simulation we will generate Y1 and Y2 to satisfy the constraint as before.
delta <- 0.01
Y1 <- vector(length = samples)
Y2 <- vector(length = samples)
for(i in 1:samples){
  Y1[i] <- par1 * V2[i]
  Y2[i] <- par2 * V2[i]
}

plot(X1, type = "l")
lines(X2, col = "red")
lines(Y1, col = "green")
lines(Y2, col = "blue")


#' we will define two agent in this cluster by two functions.
agent1 <- function(V, other_input){
  #code goes here
  return(output)
}

agent2 <- function(V, other_input){
  #code goes here
  return(output)
}
#' The idea is that both agents will describe data Y1 and Y2. For purpuse of this simulation, since we know what proces generates the data, we can define agents as this proces. This gives us both qualitative and quantitative exact describtion of data. We will include idiosynchratic shock, but will be kept on low value for now
agent1 <- function(V2){
  epsilon <- rnorm(1, mean = 0, sd = 1)
  output <- par1 * V2 + delta * epsilon
  return(output)
}

agent2 <- function(V2){
  epsilon <- rnorm(1, mean = 0, sd = 1)
  output <- par2 * V2 - delta * epsilon
  return(output)
}

#'
#' # joint density over the whole network
#'

#' Having copula over aggregation variables, margins of aggregation variables and distribution of both clusters, we can combine this to one joint density over the network. One could oblige that there is autocorrelation in agents, but we claim this autocorreltaion to be so small that it will not have significant impact on results. In practice we would need to explicity separate autocorrelation. For general agent we would have to rely on statistical methods like GARCH model and estimate the copula on residuals. This autocorrelation also produces dependency among agents, so we have to estimate copula between margins of each agent. In our example this dependece is very small so we will consider that not only clusters are conditionally independend give aggregation variables, but also variables in both clusters are conditionally independent given aggregation variable of each cluster. This simplifies the joint density formula, since  copula over particular cluster is independence copula for both clusters.
library(lmom)
dat <- as.data.frame(cbind(X1,X2,X1,X2,V1,V2))
sdX1 <- sd(X1)
sdX2 <- sd(X2)
sdY1 <- sd(Y1)
sdY2 <- sd(Y2)
sdV1 <- sd(V1)
sdV2 <- sd(V2)

density <- function(v1,v2,x1,x2,y1,y2){
  cluster1 <- dnorm(x1,mean = 0, sd = sdX1)*dnorm(x2,mean = 0, sd = sdX2)
  cluster2 <- dnorm(y1,mean = 0, sd = sdY1)*dnorm(y2,mean = 0, sd = sdY2)
#  agg_margins <- dnorm(v1,mean = 0, sd = sdV1)*dnorm(v2,mean = 0, sd = sdV2) #this is not needed because it cancels out
  output <- cluster1*cluster2*BiCopPDF(pnorm(v1,mean = 0, sd = sdV1),pnorm(v2,mean = 0, sd = sdV2), copula_obj)
  return(output)
}


#' 
#' # Causlality
#' 

#' In principle when dealing with graphical models there cannot be any causal reasoning except when we have some information on causal strucutre. Since we have causal ABM on one of the clusters, we can reason about possible response of second cluster to behaviour of agents in first cluster. 

#' (1) set state of system on ABM cluster
ag1 <- agent1(V2 = -2*sdV2)
ag2 <- agent2(V2 = -2*sdV2)
aggV <- ag1 + ag2
#' we have choasen state far from the mean
aggV

#' (2) response of V1 to value of V2 expressed as conditional probability
tmpseq <- seq(-4,2,0.001)
tmplen <- length(tmpseq)
tmp <- vector(length = tmplen)
for(i in 1:tmplen){
#  tmp[i] <- BiCopHfunc2(pnorm(tmpseq[i],mean = 0, sd = sdV1), pnorm(aggV,mean = 0, sd = sdV2), copula_obj)
  tmp[i] <- BiCopPDF(pnorm(tmpseq[i],mean = 0, sd = sdV1),pnorm(aggV,mean = 0, sd = sdV2),copula_obj)*dnorm(tmpseq[i],mean = 0, sd = sdV1)
}
plot(x=tmpseq, y=tmp, type = "l", main = "conditional density f(v1|v2)")
max(tmp) #skewed normal distiribution with heavier right tail
m<-tmpseq[which.max(tmp)]

#' (3) state of V1 cluster. Given the value of V1=m we can calculate conditional densities for X1|X1+X2=m and X2|X1+X2=m

cond_x1 <- Vectorize(function(x1){
  dnorm(x1,mean = 0, sd = sdX1)*dnorm(m-x1,mean = 0, sd = sdX2)/dnorm(m,mean = 0, sd = sdV1)
}, vectorize.args = "x1")
cond_x2 <- Vectorize(function(x2){
  dnorm(x2,mean = 0, sd = sdX2)*dnorm(m-x2,mean = 0, sd = sdX1)/dnorm(m,mean = 0, sd = sdV1)
}, vectorize.args = "x2")

tmpseq <- seq(-2.5,0.5,0.01)
plot(x=tmpseq,y=cond_x1(tmpseq), type = "l", col = "blue")
lines(x=tmpseq, y=cond_x2(tmpseq), col="green")
tmpseq[which.max(cond_x1(tmpseq))]
tmpseq[which.max(cond_x2(tmpseq))]


#'
#' # perturbation of equilibrium
#'


agent1 <- function(V2){
  if(V2 < sdV2){
    d <- 0.01
  }else{
    d <- 2
  }
  epsilon <- rnorm(1, mean = 0, sd = 1)
  output <- par1 * V2 + d * epsilon
  return(output)
}

agent2 <- function(V2){
  if(V2 < sdV2){
    d <- 0.01
  }else{
    d <- 2
  }
  epsilon <- rnorm(1, mean = 0, sd = 1)
  output <- par2 * V2 - d * epsilon
  return(output)
}

ag1series <- vector(length = samples)
ag2series <- vector(length = samples)
for(i in 1:samples){
  ag1series[i] <- agent1(V2[i])
  ag2series[i] <- agent2(V2[i])
}
plot(ag1series, type = "l")
lines(X1, col = "red")

plot(ag2series, type = "l")
lines(X2, col = "red")


