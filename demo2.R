setwd("/home/peter/Desktop/phd/R/demo")
samples <- 502 #generated series from ABM will be of length 300, there are two initial values
#set.seed(2020)

#'
#' #ABM
#'
Nfund <- 300
Ntech <- 700
N <- Ntech + Nfund
Tt <- samples
sim <- 100

Ft <- as.vector(read.table("/home/peter/Desktop/phd/R/demo/fundamentalPrice.txt")$x)
#plot(Ft, type = "l", ylim = c(8.5,11.75), xlab = "time", ylab="fundamental price")

#' market maker parameters
a <- 0.2 * 10^(-2) #a>0
sigmadelta <- 0.005 #* 10^(-10)


#trader parameters
sigmaalfa <- 1 #* 10^(-10)
sigmabeta <- 1 #* 10^(-10)
meanf <- 0.6
sdf <- 0.1
meancc <- 0.8
sdcc <- 0.4

Msim <- matrix(nrow = sim, ncol = Tt)
Dsim <- matrix(nrow = sim, ncol = Tt)
f_ind <- list()
cc_ind <- list()
defaultCount <- 0 #number of defaults
for(s in 1:sim){
  pt <- vector(length = Tt) #price adjustment
  deltap <- vector(length = Tt)
  f <- matrix(nrow = Nfund, ncol = Tt) #intensity of fundamental trader
  cc <- matrix(nrow = Ntech, ncol = Tt) #intensity of technikal trader
  Dt <- vector(length = Tt) #number of buy
  St <- vector(length = Tt) #number of sell
  Efpt <- matrix(0, nrow = Nfund, ncol = Tt) #Expectation of fundamental trader
  Ecpt <- matrix(0, nrow = Ntech, ncol = Tt) #Expectation of technikal trader
  
  #initialization
  pt[1] <- Ft[1]
    pt[2] <- Ft[2]
  
  for(t in 3:Tt){
    countBuy <- 0
    countSell <- 0
    
    f[,t] <- rnorm(Nfund, mean = meanf, sd = sdf)
    for(i in 1:Nfund){
      Efpt[i,t] <- pt[t-1] + f[i,t] * (Ft[t-1] - pt[t-1]) + rnorm(1, mean = 0, sd = sigmaalfa)
      if(Efpt[i,t] < Ft[t]){
        countBuy <- countBuy + 1
      }else if(Efpt[i,t] > Ft[t]){
        countSell <- countSell + 1
      }
    }
    
    cc[,t] <- rnorm(Ntech, mean = meancc, sd = sdcc)
    for(i in 1:Ntech){
      Ecpt[i,t] <- pt[t-1] + cc[i,t] * (pt[t-1] - pt[t-2]) + rnorm(1, mean = 0, sd = sigmabeta)
      if(Ecpt[i,t] > pt[t-1]){
        countBuy <- countBuy + 1
      }else if(Ecpt[i,t] < pt[t-1]){
        countSell <- countSell + 1
      }
    }
    Dt[t-1] <- countBuy
    St[t-1] <- countSell
    
    deltap[t] <- a*(1 + Dt[t-1] - St[t-1]) + rnorm(1, mean=0, sd = sigmadelta) 
    pt[t] <- pt[t-1] + deltap[t] #aggregation variable
    
    if(pt[t] <= 0){
      print(paste("default at time",toString(t)))
      defaultCount <- defaultCount + 1
      break;
    }
  }
  f_ind[[s]] <- f
  cc_ind[[s]] <- cc
  
  Msim[s,] <- pt
  Dsim[s,] <- deltap
}

MsimInd <- Msim
DsimInd <- Dsim

defaultCount
matplot(t(MsimInd), type="l", ylab = "value of asset", xlab = "time", ylim = c(7,12.75))
#plot(colMeans(MsimInd), ylab = "mean value of assets", xlab = "time", type = "l", ylim = c(8.5,11.75))


library(ComplexHeatmap)
library(circlize)
mInd <- MsimInd
names_m <- c("1",rep("",98),"100",rep("",99),"200",rep("",99),"300",rep("",99),"400",rep("",99),"500","","")
colnames(mInd) <- names_m
col_fun = colorRamp2(c(0, 8), c("white", "black"))
library(ggplot2)
col_fun = c("white", "black")
densityHeatmap(mInd, ylim = c(7,12.75), title = "independence", ylab = "value of asset", col = col_fun, )

#'
#' # data analysis
#'

incrInd <- as.vector(DsimInd)
hist(incrInd, breaks = 25, main = "", xlab = "increment")

upper <- quantile(incrInd, 0.90)
upper
lower <- quantile(incrInd, 0.1)
lower


#'
#' # joe copula
#'


library(copula)
joe.cop <- joeCopula(8, dim = N)


Msim <- matrix(nrow = sim, ncol = Tt)
Dsim <- matrix(nrow = sim, ncol = Tt)
f_joe <- list()
cc_joe <- list()
u01_joe <- list()
defaultCount <- 0 #number of defaults
for(s in 1:sim){
  pt <- vector(length = Tt) #price adjustment
  deltap <- vector(length = Tt)
  f <- matrix(nrow = Nfund, ncol = Tt) #intensity of fundamental trader
  cc <- matrix(nrow = Ntech, ncol = Tt) #intensity of technikal trader
  Dt <- vector(length = Tt) #number of buy
  St <- vector(length = Tt) #number of sell
  Efpt <- matrix(0, nrow = Nfund, ncol = Tt) #Expectation of fundamental trader
  Ecpt <- matrix(0, nrow = Ntech, ncol = Tt) #Expectation of technikal trader
  
  #initialization
  pt[1] <- Ft[1]
  pt[2] <- Ft[2]
  
  u01 <- rCopula(Tt, joe.cop)
  for(t in 3:Tt){
    countBuy <- 0
    countSell <- 0
    
    #    f[,t] <- rnorm(Nfund, mean = meanf, sd = sdf)
    for(i in 1:Nfund){
      f[i,t] <- qnorm(u01[t,i], mean = meanf, sd = sdf) #values u01[t,1:Nfund]; Q[N(meanf,sdf)] 
      #      f[i,t] <- 0.7
      #      f[i,t] <- rnorm(1, mean = meanf, sd=sdf)
      Efpt[i,t] <- pt[t-1] + f[i,t] * (Ft[t-1] - pt[t-1]) + rnorm(1, mean = 0, sd = sigmaalfa)
      if(Efpt[i,t] < Ft[t]){
        countBuy <- countBuy + 1
      }else if(Efpt[i,t] > Ft[t]){
        countSell <- countSell + 1
      }
    }
    
    for(i in 1:Ntech){ #values u01[t,Nfund+1:N]; Q[N(meanf,sdf)]
      cc[i,t] <- qnorm(u01[t,Nfund + i], mean = meancc, sd = sdcc) 
      #     cc[i,t] <- 1.2
      Ecpt[i,t] <- pt[t-1] + cc[i,t] * (pt[t-1] - pt[t-2]) + rnorm(1, mean = 0, sd = sigmabeta)
      if(Ecpt[i,t] > pt[t-1]){
        countBuy <- countBuy + 1
      }else if(Ecpt[i,t] < pt[t-1]){
        countSell <- countSell + 1
      }
    }
    Dt[t-1] <- countBuy
    St[t-1] <- countSell
    
    deltap[t] <- a*(1 + Dt[t-1] - St[t-1]) + rnorm(1, mean=0, sd = sigmadelta) 
    pt[t] <- pt[t-1] + deltap[t] #aggregation variable
    
    if(pt[t] <= 0){
      print(paste("default at time",toString(t)))
      defaultCount <- defaultCount + 1
      break;
    }
  }
  f_joe[[s]] <- f
  cc_joe[[s]] <- cc
  u01_joe[[s]] <- u01
  
  Msim[s,] <- pt
  Dsim[s,] <- deltap
}

defaultCount

orderM <- order(sqrt(diag(cov(t(Msim)))), decreasing = TRUE)[1:5]

Msim[Msim > 11.5, ]

orderM <- which(Msim > 12, arr.ind = TRUE)
indexes <- unique(orderM[,1])
orderM <- which(Msim < 8, arr.ind = TRUE)
indexes <- c(indexes ,unique(orderM[,1]))

matplot(t(Msim[sample(indexes, 5),]), type="l", ylab = "value of asset", xlab = "time", ylim = c(7,12.75))
#plot(colMeans(Msim), ylab = "mean value of assets", xlab = "time", type = "l", ylim = c(8.5,11.75))


mSim <- Msim
names_m <- c("1",rep("",98),"100",rep("",99),"200",rep("",99),"300",rep("",99),"400",rep("",99),"500","","")
colnames(mSim) <- names_m
col_fun = colorRamp2(c(0, 8), c("white", "black"))
library(ggplot2)
col_fun = c("white", "black")
densityHeatmap(mSim, ylim = c(7,12.75), title = "dependence", ylab = "value of asset", col = col_fun)


#'
#' # data analysis
#'

pacf(abs(DsimInd[1,]), lag.max = 15)
pacf(abs(Dsim[1,]), lag.max = 15)

  incr <- as.vector(Dsim)
hist(incr, breaks = 25, main = "", xlab = "increment")

upper <- quantile(incr, 0.90)
upper
lower <- quantile(incr, 0.1)
lower


#' qqplots

qqnorm(incrInd, main = "Independence")
qqline(incrInd, col = "steelblue", lwd = 2)
qqnorm(incr, main = "Dependence")
qqline(incr, col = "steelblue", lwd = 2)


shapiro.test(sample(incrInd,5000))
shapiro.test(sample(incr,5000))


