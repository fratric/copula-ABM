setwd("C:/Users/Peter/Desktop/phd/R/demo")
samples <- 3 #generated series from ABM will be of length 300, there are two initial values
#set.seed(2020)

#'
#' #ABM
#'
Nfund <- 300
Ntech <- 700
N <- Ntech + Nfund
Tt <- samples

Ft <- as.vector(read.table("C:/Users/Peter/Desktop/phd/R/demo/fundamentalPrice.txt")$x)
#plot(Ft, type = "l", ylim = c(8.5,11.75), xlab = "time", ylab="fundamental price")


#'
#' # joe copula
#'


library(copula)


simulation <- function(a, sigmadelta, sigmaalfa, sigmabeta, meanf, sdf, meancc, sdcc, pt1, pt2, joeParam){
  defaultCount <- 0 #number of defaults
 
  
  pt <- vector(length = Tt) #price adjustment
  deltap <- vector(length = Tt)
  f <- matrix(nrow = Nfund, ncol = Tt) #intensity of fundamental trader
  cc <- matrix(nrow = Ntech, ncol = Tt) #intensity of technikal trader
  Dt <- vector(length = Tt) #number of buy
  St <- vector(length = Tt) #number of sell
  Efpt <- matrix(0, nrow = Nfund, ncol = Tt) #Expectation of fundamental trader
  Ecpt <- matrix(0, nrow = Ntech, ncol = Tt) #Expectation of technikal trader
  
  #initialization
  pt[1] <- pt1
  pt[2] <- pt2
  
  joe.cop <- joeCopula(joeParam, dim = N)
  u01 <- rCopula(Tt, joe.cop)
  for(t in 3:Tt){
    countBuy <- 0
    countSell <- 0
    
    
    #herding behaviour description
#    if(t == T1){
#      joe.cop <- joeCopula(7, dim = N)
#      u01 <- rCopula(Tt, joe.cop)
#    }else if(t == T2){
#      joe.cop <- joeCopula(1.02, dim = N)
#      u01 <- rCopula(Tt, joe.cop)
#    }
    
    # brocken-homes paper
  
    for(i in 1:Nfund){
      f[i,t] <- qnorm(u01[t,i], mean = meanf, sd = sdf) #values u01[t,1:Nfund]; Q[N(meanf,sdf)] 

      Efpt[i,t] <- pt[t-1] + f[i,t] * (Ft[t-1] - pt[t-1]) + rnorm(1, mean = 0, sd = sigmaalfa)
      if(Efpt[i,t] < Ft[t]){
        countBuy <- countBuy + 1
      }else if(Efpt[i,t] > Ft[t]){
        countSell <- countSell + 1
      }
    }
    
    for(i in 1:Ntech){ #values u01[t,Nfund+1:N]; Q[N(meanf,sdf)]
      cc[i,t] <- qnorm(u01[t,Nfund + i], mean = meancc, sd = sdcc) 
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
  
  

  returnList <- list("pt" = pt, "deltap" = deltap)
  return(returnList)
}


#' #' market maker parameters
#' a <- 0.2 * 10^(-2) 
#' sigmadelta <- 0.005
#' 
#' 
#' #trader parameters
#' sigmaalfa <- 1 #* 10^(-10)
#' sigmabeta <- 1 #* 10^(-10)
#' meanf <- 0.6
#' sdf <- 0.1
#' meancc <- 0.8
#' sdcc <- 0.4

#' market maker parameters
a <- 0.5 * 10^(-2) 
sigmadelta <- 0.005


#trader parameters
sigmaalfa <- 1 #* 10^(-10)
sigmabeta <- 1 #* 10^(-10)
meanf <- 0.6
sdf <- 0.1
meancc <- 0.8
sdcc <- 0.4



sim <- 1000

MsimX <- matrix(nrow = sim, ncol = Tt)
deltaPsimX<- matrix(nrow = sim, ncol = Tt)

for(s in 1:sim){
  x <- simulation(a = a, sigmadelta = sigmadelta, sigmaalfa = sigmaalfa, sigmabeta = sigmabeta, meanf = meanf, sdf = sdf, meancc = meancc, sdcc = sdcc, pt1 = 10, pt2 = 13, joeParam = 1.02)
  MsimX[s,] <- x$pt
  deltaPsimX[s,] <- x$deltap
}
hist(deltaPsimX)
matplot(t(deltaPsimX[,2:3]), type="l", ylab = "delta p", xlab = "time", ylim = c(0,3), main = "Independence")



for(s in 1:sim){
  x <- simulation(a = a, sigmadelta = sigmadelta, sigmaalfa = sigmaalfa, sigmabeta = sigmabeta, meanf = meanf, sdf = sdf, meancc = meancc, sdcc = sdcc, pt1 = 10, pt2 = 13, joeParam = 7)
  MsimX[s,] <- x$pt
  deltaPsimX[s,] <- x$deltap
}
hist(deltaPsimX)
matplot(t(deltaPsimX[,2:3]), type="l", ylab = "delta p", xlab = "time", ylim = c(0,3), main = "Joe")



#'
#' # data analysis
#'
 
# incr <- as.vector(Dsim)
# hist(incr, breaks = 25, main = "", xlab = "increment")
# 
# upper <- quantile(incr, 0.90)
# upper
# lower <- quantile(incr, 0.1)
# lower


#' qqplots

# qqnorm(incrInd)
# qqline(incrInd, col = "steelblue", lwd = 2)
# qqnorm(incr)
# qqline(incr, col = "steelblue", lwd = 2)
# 
# 
# shapiro.test(sample(incrInd,5000))
# shapiro.test(sample(incr,5000))


#'
#' #multiple markets - independence bettween groups of predictor variables
#'

tau <- 0.5
theta <- iTau(gumbelCopula(), tau = tau)
d <- 2
gc <- gumbelCopula(theta, dim = d)
n <- 200
U. <- matrix(runif(n*d), ncol = d) # U(0,1)^d

## Transform to Gumbel sample via conditional distribution method
U <- cCopula(U., copula = gc, inverse = TRUE) 
# slow for ACs except Clayton
splom2(U) # scatter-plot matrix copula sample


## Rosenblatt transform back to U(0,1)^d (as a check)
U. <- cCopula(U, copula = gc)
splom2(U.) # U(0,1)^d again


#'
#' #estimation of aggregation copula from simulations agent-based models
#'



#'
#' #evaluation of conditional probability
#'

#'estimate marginal distribution function

library(lmom)
V1cdf_ind <- function(x) cdfnor(x, pelnor(samlmu(incrInd)))
plot(incrInd, copula::pobs(incrInd),  ylab="distribution function", xlab="quantile")
tmpx <- seq(min(incrInd), max(incrInd), 0.01)  
lines(tmpx, V1cdf_ind(tmpx), col="red")

library(MASS)
library(stable)
V1cdf_joe <- function(x) pstable(x, disp = 0.131)
plot(incr, copula::pobs(incr),  ylab="distribution function", xlab="quantile")
tmpx <- seq(min(incr), max(incr), 0.01)
lines(tmpx, V1cdf_joe(tmpx), col="red")

s <- 1 #simulation
#g01 <- list() #synthetic data
#g01[[s]] <- rCopula(Tt, t.cop)




# likelihood_ratio <- matrix(0, nrow = s, ncol = Tt)
# t <- 3
# sequenc1 <- seq(0.01,0.99,0.01)
# tmp <- vector(length = length(sequenc1))
# index <- 1
# 
# for(v2 in sequenc1){
#   
#   v1indep <- DsimInd[s,t]
#   FV1v1 <- V1cdf_ind(v1indep)
#   cind <- dCopula(c(FV1v1,v2), t.cop) #global copula - t, margin distributed normally
#   
#   v1joe <- Dsim[s,t]
#   FV1v1 <- V1cdf_joe(v1joe)
#   cjoe <- dCopula(c(FV1v1,v2), t.cop) #global copula - t, margin distributed t
#   
# 
#   clocaljoe <- dCopula(u01_joe[[s]][t,] , joe.cop) #local copula - joe
#   clocaljoe1 <- dCopula(u01_joe[[s]][t,] , joe.cop1)
#   
#   tmp[index] <- clocaljoe*cjoe/(cind*clocaljoe1)
#   index <- index + 1
# }
# plot(x = seq(0.01,0.99,0.01), y = tmp, type = "l")
# 
# dCopula(rep(0.5,200) , joe.cop)
# 
# library(copula)
# dCopula(rep(0.5,200) , joe.cop)

# f=N(0.2,0.1); cc = 0.1   -very certain
# f=N(0.2,0.1); cc = 1.1   -very volatile
# f=N(0.2,0.1); cc = 0.6   -uncertain around mean, but not volatile

#todo: 
# 0. rewrite deltap as function to see the form aggregation function
# 1. play with psychological parameters f, cc distributed normally (or differently) 
# 2. bayes theorem - given f,cc what is distribution of network. Are there similarities between bayes theorem and conditional desity over the network?
# 3. confounding (same variable influences two clusters)
# 4. vytried hodnoty deltap kde nastal extremny rast/pokles a najdi k nim prisluchajuce hodnoty ff alebo c, potom sa bavme o tom ako bivariate copula prenasa toto extremne spravanie
