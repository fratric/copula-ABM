setwd("C:/Users/Peter/Desktop/phd/R/demo")
samples <- 800
set.seed(2019)

#'
#' # Simulation case study in equilibrium state
#'


#' define potential. Each agent will reside in local minima of his potential function.
potential <- function(x){
  return(((x-0.25)^2)*(x-1)*(x-1.5)+1)
}

tmpx <- seq(0,1.6,0.01)
tmpy <- vector(length = length(tmpx))
for(i in 1:length(tmpx)){
  tmpy[i] <- potential(tmpx[i])
}
plot(x=tmpx, y=tmpy, type = "l")


#' define agents in ABM cluster. Agents oscilate around local minima of their potential
Y1_locmin <- 0.25
Y2_locmin <- 1+sqrt(3/2)/4

agent1 <- function(){
  epsilon <- rnorm(1, mean = 0, sd = 1)
  delta <- 0.2
  output <- Y1_locmin + delta * epsilon
  return(output)
}

agent2 <- function(){
  epsilon <- rnorm(1, mean = 0, sd = 1)
  delta <- 0.2
  output <- Y2_locmin + delta * epsilon
  return(output)
}

#' aggregation variable V2 will represent potential energy sum of both agents
V2 <- rep(0, samples)
Y1 <- vector(length = samples)
Y2 <- vector(length = samples)
for(i in 1:samples){
  Y1[i] <- agent1()
  Y2[i] <- agent2()
  V2[i] <- potential(Y1[i]) + potential(Y2[i])
}
plot(Y1, type="l")
lines(rep(0.25,samples),col="green")
hist(Y1, breaks = 20)
pacf(Y1)

plot(Y2, type="l")
lines(rep(1+sqrt(3/2)/4,samples),col="red")
hist(Y2)
pacf(Y2)

plot(V2, type = "l")
lines(2*rep(potential(0.25),samples),col="green")
lines(rep(potential(0.25)+potential(1+sqrt(3/2)/4),samples),col="red")
hist(V2, breaks = 30)
pacf(V2)

YDF <- as.data.frame(cbind(Y1,Y2,V2))


#plot(density(V2))
#df <- approxfun(density(V2))
#tmpx <- seq(1.90,2.3,0.0001)
#Mode <- tmpx[which.max(df(tmpx))]

Ecdf <- ecdf(V2)
sorted_V2 <- sort(V2)
e_cdf <- 1:length(sorted_V2) / length(sorted_V2)

plot(sorted_V2, e_cdf, type = "s")
abline(h = 0.5, lty = 3)
Mode<-sorted_V2[which(e_cdf >= 0.5)[1]]
abline(h = 0.99, lty = 3)
quantile99 <- sorted_V2[which(e_cdf >= 0.99)[1]]

Ecdf(quantile99)
Ecdf(Mode)

#' Let V1 be aggregation variable on second cluster having normal distribution. V1 will aggregate normally distributed variables in cluster linked by multivariabe vine copula
library(VineCopula)

# define 5-dimensional R-vine tree structure matrix
Matrix <- c(5, 2, 3, 1, 4,0, 2, 3, 4, 1,0, 0, 3, 4, 1,0, 0, 0, 4, 1,0, 0, 0, 0, 1)
Matrix <- matrix(Matrix, 5, 5)

# define R-vine pair-copula familymatrix
family <- c(0, 1, 3, 4, 4,0, 0, 3, 4, 1,0, 0, 0, 4, 1,0, 0, 0, 0, 3,0, 0, 0, 0, 0)
family <- matrix(family, 5, 5)

# define R-vine pair-copula parameter matrix
par <- c(0, 0.2, 0.9, 1.5, 3.9,0, 0, 1.1, 1.6, 0.9,0, 0, 0, 1.9, 0.5,0, 0, 0, 0, 4.8,0, 0, 0, 0, 0)
par <- matrix(par, 5, 5)

# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 5, 5)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,par = par, par2 = par2,names = c("X1", "X2", "X3", "X4", "X5"))
RVM

#' Once we have our vine copula defined. We can draw samples from it and define aggregation function h1 and observe distribution V1=h1(X1,..,X5). We will use 10000 samples to get better idea of how V1 looks like

simdata <- RVineSim(50000, RVM)
 
#' we will set marginal distributions for some specific distribution function. To demonstrate the fact that we can use any type of distribution we will mix different distributions or different parameters of distributions for each variable

x1 <- qnorm(simdata[,1]) # normal with mean=0, sd=1
x2 <- qlnorm(simdata[,2]) #lognormal with meanlog = 0, sdlog = 1
x3 <- qt(simdata[,3], df=5, ncp=2) #t with 5 degrees of freedom and noncentraility parameter 2 
x4 <- qt(simdata[,3], df=4, ncp=1) #t with 4 degrees of freedom and noncentraility parameter 1
x5 <- qnorm(simdata[,5], mean = 0, sd = 0.5) # normal with mean=0, sd=1/2

v1 <- x1+x2+x3+x4+x5
v1min <- which.max(v1)
hist(v1, breaks = 50)


library(lmom)
pV1 <- function(x) cdfgev(x, pelgev(samlmu(v1)))
tmpx <- seq(min(v1), max(v1), 1)  
plot(v1, copula::pobs(v1))
lines(x=tmpx,y=pV1(tmpx), col = "red")


qV1 <- function(x) quagev(x, pelgev(samlmu(v1)))
xdf <- as.data.frame(cbind(x1,x2,x3,x4,x5,v1))

#' let U1 and U2 be uniformly distributed random variables created from V1 and V2 by probability integral transformation. We will define bivariate distribution of U1 and U2 by bi-copula from student copula family with 6 degrees of freedom and strong tail dependency
copula_obj <- BiCop(family = 2, par = 0.8, par2 = 6)
copula_obj


#' Observe causal relationshipe Y1,Y2 -> V2 and assume causal relationship V2 -> V1 -> X1,X2,X3,X3,X4. Further on we will investigate causal propagation on this network. For simplicity we will assume that V1,X1,X2,X3,X4,X5 are normally distributed N(0,1)

conditionalDensity <- function(x1, x2, x3, x4, x5, y1, y2){
  #conditional density over whole network conditioned on Y=y
  v1 <- x1+x2+x3+x4+x5
  v2 <- potential(y1) + potential(y2)
  
  
  marginsVine <- dnorm(x1)*dlnorm(x2)*dt(x3, df=5, ncp=2)*dt(x4, df=4, ncp=1)*dnorm(x5, mean = 0, sd = 0.5)
  
  vinePDF <- RVinePDF(c(pnorm(x1), plnorm(x2), pt(x3, df=5, ncp=2), pt(x4, df=4, ncp=1), pnorm(x5,mean = 0, sd = 0.5)), RVM)
  
  output <- BiCopPDF(pV1(v1), Ecdf(v2), copula_obj)*vinePDF*marginsVine 
  return(output)
}

#' we can see in expression above that copula acts as a scaling parameter for distribution on V1 cluster
tmpx <- seq(0.01,0.99,0.01)
tmpy <- vector(length = length(tmpx))
for(i in 1:length(tmpx)){
  tmpy[i] <- BiCopPDF(tmpx[i], Ecdf(Mode), copula_obj)
}
plot(x=tmpx, y=tmpy, type = "l")


#' now let us generate samples of V1 from conditional bivariate copula conditioned on V2

V1 <- BiCopCondSim(samples, Ecdf(V2), cond.var = 2, copula_obj)

X1 <- vector(length = samples)
X2 <- vector(length = samples)
X3 <- vector(length = samples)
X4 <- vector(length = samples)
X5 <- vector(length = samples)

#we use fitted distribution on generated and aggregated data from vine copula to transfrom simulated data from bivariate copula to [0,1]
V1 <- qV1(V1)

min(V1)
min(v1)

max(V1)
max(v1)



cou<-0
for(i in 1:samples){
  #for given tolerance check the situation
  if(length(v1[abs(v1-V1[i])<0.01])>0){
    cou<-cou+1
  }
}
cou

for(i in 1:samples){
  # for each value find the closest element
  index <- which.min(abs(v1 - V1[i]))
  
  X1[i] <- x1[index]
  X2[i] <- x2[index]
  X3[i] <- x3[index]
  X4[i] <- x4[index]
  X5[i] <- x5[index]
}

XDF <- as.data.frame(cbind(X1,X2,X3,X4,X5,V1))

plot(V1[order(V1)])
lines((X1+X2+X3+X4+X5)[order(X1+X2+X3+X4+X5)], col = "red", type="b")

DF <- as.data.frame(cbind(Y1,Y2,V2,X1,X2,X3,X4,X5,V1))

#' since we have data generated we can plot a correlation graph
library(qgraph)
#tu sa pohraj s correlacnym koefictinetom (tau) a so zobrazovanim a clustrovanim
qgraph(cor(DF, method = "kendall"), layout="circle", label.cex=0.9, label.scale=FALSE)


#' now let's go back to our problem ofdetermining the most probable set of states given known state on ABM cluster

condDen <- function(vectorX){
  return(conditionalDensity(vectorX[1], vectorX[2], vectorX[3], vectorX[4], vectorX[5], Y1_locmin, Y2_locmin))
}

colMeans(XDF)[-6]

Optima <- optim(colMeans(XDF)[-6], condDen, method = "SANN", control = list(fnscale = -1))
Optima

pdfC <- function(x){
#  return(RVinePDF(x, RVM))
  return(RVinePDF(c(pnorm(x[1]), pnorm(x[2]), pt(x[3], df=5, ncp=2), pt(x[4], df=4, ncp=1), pnorm(x[5],mean = 0, sd = 0.5)), RVM))
}

# FUCK FUCK FUCK

pdfC(c(0.5,0.5,0.5,0.5,0.5))

optim(colMeans(XDF)[-6],pdfC, control = list(fnscale = -1))


condDen(unname(Optima$par))
conditionalDensity(-3.9,0.08,2.29,-7.2,-1, Y1_locmin, Y2_locmin)

