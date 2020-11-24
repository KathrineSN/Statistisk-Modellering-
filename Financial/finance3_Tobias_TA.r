########################################################################################################################
rm(list=ls())
setwd("~/Documents/DTU/Matematik/Statistisk_modellering/Assignments/3")
library(MASS)
library(car)
library(numDeriv)
library(mvtnorm)

########################################################################################################################
## Import the data

## Read the finans_data.csv file containing the data
wr <- read.csv("finance_data.csv", header=TRUE, sep=";")

########################################################################################################################
## Overview of the data

## Dimension of HE (number of rows and columns)
dim(wr)
## Column names
names(wr)
## The first rows
head(wr)
## The last rows
tail(wr)
## Default summary
summary(wr)
## Another summary function also including the data type
str(wr)


########################################################################################################################
## 1: Descriptive statistics and initial analysis
########################################################################################################################

## a) ##################################################################################################################

## Presenting data 
time <- wr[, 1]
time <- as.Date(time)
SPY <- wr[, 2]
x <- seq(1, length(time))

par(mfrow=c(1,1))
plot(time, SPY, type = 'l', col = "blue", xlab = 'time', ylab = 'SPY', main = 'SPY as a function of time')

## Estimate paramters in normal model 

## Profile likelihoods 
lp.mu <- function(mu, x){
  n <- length(x)
  sigma2.hat <- sum((x - mu)^2)/n
  -n/2*log(sigma2.hat) - 1/(2*sigma2.hat)*sum((x - mu)^2)
}

lp.sigma2 <- function(sigma2, x){
  n <- length(x)
  sigma2.hat <- sum((x - mean(x))^2)/n
  -n/2*log(sigma2) - sigma2.hat * n/(2 * sigma2)
}

## Plot 
par(mfrow=c(2,1))

mu <- seq(-mean(wr$SPY),mean(wr$SPY)*3,len=length(wr$SPY))
plot(mu, exp(sapply(mu, lp.mu, wr$SPY) - lp.mu(mean(wr$SPY), wr$SPY)), type = "l")

sigma2 <- seq(var(wr$SPY)/1.25,var(wr$SPY)*1.25,len=length(wr$SPY))
plot(sigma2, exp(sapply(sigma2, lp.sigma2, wr$SPY) - lp.sigma2(var(wr$SPY), wr$SPY)), type = "l")

## Negative log lilelihood (numerical method)
theta = c(0,0.05)
l.norm = function(theta, x) {
  -sum(dnorm(x, theta[1], theta[2], log = TRUE))
}

fit.norm = optim(theta, l.norm, NULL, SPY)
fit.norm$par

## Check distribution
par(mfrow=c(1,1))
hist(SPY, prob=TRUE, nclass=60)
lines(seq(min(SPY), max(SPY), 0.001), 
      dnorm(seq(min(SPY), max(SPY), 0.001), fit.norm$par[1], fit.norm$par[2]), 
      lwd = 1, col="2")

qqPlot(SPY)

AIC.norm <- 2*l.norm(fit.norm$par, SPY) + 2*2


## b) ##################################################################################################################



########################################################################################################################
## 2: Mixture model
########################################################################################################################

## a) ##################################################################################################################

acf(SPY)

## The autocorrelation plot for weekly returns (SPY) shows that most of the spikes are not statistically significant. 
## This indicates that the returns are not highly correlated, as shown here.
## They mostly fluctuate around zero correlation.

########################################################################################################################

## Normal mixture: transform
## Natural to working parameters
norm.mix.pn2pw <- function(m, mu, sigma, delta) {
  if(sum(delta) >= 1) {
    print("sum(delta) should be < 1")
    return()
  }
  t.sigma <- log(sigma)
  t.delta <- log(delta/(1 - sum(delta)))
  return(list(mu = mu, t.sigma = t.sigma, t.delta = t.delta))
}

## Working to natural parameters
norm.mix.pw2pn <- function(m, mu, t.sigma, t.delta){
  if(m == 1){
    return(exp(t.sigma))
  }
  sigma <- exp(t.sigma)
  delta <- exp(t.delta)/(1 + sum(exp(t.delta)))
  delta <- c(1 - sum(delta),delta)
  return(list(mu = mu, sigma = sigma, delta = delta))
}

## Negative log likelihood
nll <- function(theta, m, x){
  if(m == 1) {
    return(-sum(pnorm(x, theta[1], exp(theta[2]), log=TRUE))) 
  }
  mu <- theta[1:m]
  t.sigma <- theta[(m+1):(2*m)]
  t.delta <- theta[(2*m+1):(3*m-1)]
  n.pars <- norm.mix.pw2pn(m, mu, t.sigma, t.delta)
  n <- length(x)
  nll <- 0
  for(i in 1:n) {
    nll <- nll - log(sum(n.pars$delta * dnorm(x[i], mu, n.pars$sigma)))
  }
  return(nll)
}

########################################################################################################################

## Estimation with 2 distributions
m <- 2; 

## Initial values
mu <- mean(SPY)*c(1/2,3/2)
sigma <-sd(SPY)*c(1/2,3/2)
delta <- c(1/2)

## Working parameters
wpars <- norm.mix.pn2pw(m, mu, sigma, delta)
theta <- c(wpars$mu, wpars$t.sigma, wpars$t.delta)

## MLE
opt2 <- nlminb(theta, nll, m = m, x = SPY)

## Natural parameters
npars2 <- norm.mix.pw2pn(m, opt2$par[1:m], opt2$par[(m+1):(2*m)], opt2$par[(2*m+1):(3*m-1)])

########################################################################################################################

## Estimation with 3 distributions
m <- 3;

## Initial values 
mu <- mean(SPY)*c(1/2,1,3/2)
sigma <- sd(SPY)*c(1/2,1,3/2);
delta <- c(1/3,1/3)

## Working parameters
wpars <- norm.mix.pn2pw(m, mu, sigma, delta)
theta <- c(wpars$mu, wpars$t.sigma, wpars$t.delta)

## MLE
opt3 <-nlminb(theta, nll, m = m, x = SPY)

## Natural parameters
npars3 <- norm.mix.pw2pn(m, opt3$par[1:m], opt3$par[(m+1):(2*m)], opt3$par[(2*m+1):(3*m-1)])

########################################################################################################################

mix.dist <- function(x ,npars){
  sum(npars$delta * dnorm(x, mean = npars$mu, sd = npars$sigma))
}

## Plot
par(mfrow=c(1,1))
hist(SPY, prob=TRUE, nclass=60)
lines(seq(min(SPY), max(SPY), 0.001), sapply(seq(min(SPY), max(SPY), 0.001), mix.dist, npars=npars2), col=2)
lines(seq(min(SPY), max(SPY), 0.001), 0.8928764*dnorm(seq(min(SPY), max(SPY), 0.001), 0.002607149, 0.01878823), col=3)
lines(seq(min(SPY), max(SPY), 0.001), 0.1071236*dnorm(seq(min(SPY), max(SPY), 0.001), -0.009033854, 0.05156611), col=3)
lines(seq(min(SPY), max(SPY), 0.001), sapply(seq(min(SPY), max(SPY), 0.001), mix.dist, npars=npars3), col=4)
legend("topleft",  c("2 components","3 components"), col=c("red", "blue"), lty=1, cex=0.7)


## Model check
AIC <- 2*c(opt2$objective, opt3$objective) + 2*c(length(opt2$par), length(opt3$par))

## Deviance 
1-pchisq(-2*(opt3$objective-opt2$objective),df=length(opt3$par)-length(opt2$par))


## b) ##################################################################################################################

## CI for working parameters (t.sigma and t.delta) 

## Using Wald while in unrestricted domain
H <- nlm(nll, theta, m = m, x = SPY, hessian = TRUE)$hessian
V.pars <- solve(H)
se <- sqrt(diag(V.pars))
tab <- round(cbind(opt3$par, se), digit=10)
rownames(tab) <- c("mu1","mu2","mu3",
                   "t.sigma1","t.sigma2","t.sigma2",
                   "t.delta1","t.delta2")
tab

########################################################################################################################

## CI for natural parameters by simulation from Wald distribution (normal)
k <- 100000
PARS <- rmvnorm(k, mean=opt3$par, sigma = V.pars)

## Simulated (mu1)
(CI.mu <- quantile(PARS[ ,1], probs=c(0.025,0.975)))
## Wald based
opt3$par[1]+qnorm(c(0.025,0.975))*se[1]

## Simulated (mu2)
CI.mu <- rbind(CI.mu, quantile(PARS[ ,2], probs=c(0.025,0.975)))
CI.mu[2, ]
## Wald based
opt3$par[2]+qnorm(c(0.025,0.975))*se[2]

## Simulated (mu3)
CI.mu <- rbind(CI.mu, quantile(PARS[ ,3], probs=c(0.025,0.975)))
CI.mu[3, ]
## Wald based
opt3$par[3]+qnorm(c(0.025,0.975))*se[3]


## Simulated (sigma1)
(CI.sigma <- quantile(exp(PARS[ ,4]), probs=c(0.025,0.975)))
## Wald based
exp(opt3$par[4]+qnorm(c(0.025,0.975))*se[4])

## Simulated (sigma2)
CI.sigma <- rbind(CI.sigma,quantile(exp(PARS[ ,5]), probs=c(0.025,0.975)))
CI.sigma[2, ]
## Wald based
exp(opt3$par[5]+qnorm(c(0.025,0.975))*se[5])

## Simulated (sigma3)
CI.sigma <- rbind(CI.sigma,quantile(exp(PARS[ ,6]), probs=c(0.025,0.975)))
CI.sigma[3, ]
## Wald based
exp(opt3$par[6]+qnorm(c(0.025,0.975))*se[6])


## Simulated (delta1, delta2 and delta3)
delta2 <- exp(PARS[ ,7])/(1 + rowSums(exp(PARS[ ,7:8])))
delta3 <- exp(PARS[ ,8])/(1 + rowSums(exp(PARS[ ,7:8])))
delta1 <- 1-delta2-delta3

## CI for delta
CI.delta <- c(npars3$delta[1], quantile(delta1, probs=c(0.025,0.5,0.975)))
CI.delta <- rbind(CI.delta, c(npars3$delta[2], quantile(delta2, probs=c(0.025,0.5,0.975))))
CI.delta <- rbind(CI.delta, c(npars3$delta[3], quantile(delta3, probs=c(0.025,0.5,0.975))))
CI.delta

## c) ##################################################################################################################

## Profile likelihood for sigma 1 given working parameters
lp.sigma1 <- function(sigma1, m, x, pars0){
  ## Fun for inner optim
  fun.tmp <- function(theta, sigma1, m, x){
    pars <- c(theta[1:m], log(sigma1), theta[-(1:m)])
    nll(pars, m, x)
  }
  nlminb(pars0, fun.tmp, sigma1 = sigma1, m = m, x = x)$objective    
}

## Estimation with 2 distributions
m <- 2; 

## Initial values
mu <- mean(SPY)*c(1/2,3/2)
sigma <-sd(SPY)*c(1/2,3/2)
delta <- c(0.1)

## Working parameters
wpars <- norm.mix.pn2pw(m, mu, sigma, delta)
theta0 <- c(wpars$mu, wpars$t.sigma, wpars$t.delta)
theta <- c(theta0[1:m],theta0[(m+2):(3*m-1)])

sigma1 <- seq(0.01, 0.1, length=100)

## profile likeihood
pnll <- sapply(sigma1, lp.sigma1, m = m, x = SPY, pars0 = theta)

## Plot the profile likelihood
plot(sigma1,exp(-(pnll-min(pnll))),type="l", ylim=c(0,1))
lines(range(sigma1),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),col=2,lty=2,lwd=2)
rug(npars2$sigma,col=2,lwd=2)

## d) ##################################################################################################################

## Reprametrizing mix-model 

## Natural to working parameters
re.norm.mix.pn2pw <- function(m, mu, sigma, delta) {
  if(sum(delta) >= 1) {
    print("sum(delta) should be < 1")
    return()
  }
  t.sigma <- log(c(sigma[1],sigma[-1]-sigma[-m]))
  t.delta <- log(delta/(1 - sum(delta)))
  return(list(mu = mu, t.sigma = t.sigma, t.delta = t.delta))
}

## Working to natural parameters
re.norm.mix.pw2pn <- function(m, mu, t.sigma, t.delta){
  if(m == 1){
    return(exp(t.sigma))
  }
  sigma <- cumsum(exp(t.sigma))
  delta <- exp(t.delta)/(1 + sum(exp(t.delta)))
  delta <- c(1 - sum(delta),delta)
  return(list(mu = mu, sigma = sigma, delta = delta))
}

## Negative log likelihood
re.nll <- function(theta, m, x){
  if(m == 1) {
    return(-sum(pnorm(x, theta[1], exp(theta[2]), log=TRUE))) 
  }
  mu <- theta[1:m]
  t.sigma <- theta[(m+1):(2*m)]
  t.delta <- theta[(2*m+1):(3*m-1)]
  n.pars <- re.norm.mix.pw2pn(m, mu, t.sigma, t.delta)
  n <- length(x)
  nll <- 0
  for(i in 1:n) {
    nll <- nll - log(sum(n.pars$delta * dnorm(x[i], mu, n.pars$sigma)))
  }
  return(nll)
}

#########################################################################################################################

## Profile likelihood for sigma 1 given working parameters
re.lp.sigma1 <- function(sigma1, m, x, pars0){
  ## Fun for inner optim
  fun.tmp <- function(theta, sigma1, m, x){
    pars <- c(theta[1:m], log(sigma1), theta[-(1:m)])
    re.nll(pars, m, x)
  }
  nlminb(pars0, fun.tmp, sigma1 = sigma1, m = m, x = x)$objective    
}

m <- 2; 

## Initial values
mu <- mean(SPY)*c(1/2,3/2)
sigma <-sd(SPY)*c(1/2,3/2)
delta <- c(0.1)

## Working parameters
wpars <- re.norm.mix.pn2pw(m, mu, sigma, delta)
theta0 <- c(wpars$mu, wpars$t.sigma, wpars$t.delta)
theta <- c(theta0[1:m],theta0[(m+2):(3*m-1)])

sigma1 <- seq(0.01, 0.05, length=100)

## profile likeihood
pnll <- sapply(sigma1, re.lp.sigma1, m = m, x = SPY, pars0 = theta)

## Plot the profile likelihood
plot(sigma1,exp(-(pnll-min(pnll))),type="l", ylim=c(0,1))
lines(range(sigma1),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),col=2,lty=2,lwd=2)
rug(npars2$sigma[1],col=2,lwd=2)


## e) ##################################################################################################################





########################################################################################################################
## 3: Hidden Markov Models
########################################################################################################################

## a) ##################################################################################################################

## Hidden Markov Model for Normal Distribution 
norm.HMM.pn2pw <- function(m, mu, sigma, gamma) {                                              
  t.sigma <- log(sigma)
  t.gamma <- NULL                              
  if (m > 1) {                                            
    foo     <- log(gamma / diag(gamma))           
    t.gamma <- as.vector(foo[!diag(m)])             
  }                                             
  parvect <- c(as.vector(mu), as.vector(t.sigma), 
               as.vector(t.gamma))
  parvect
}  

norm.HMM.pw2pn <- function(m, parvect) {
  mu    <- parvect[1:m]
  epar  <- exp(parvect[-(1:m)]) 
  sigma <- epar[1:m]
  gamma <- diag(m)                                    
  if (m > 1) {                                                  
    gamma[!gamma] <- epar[(m+1):(m*m)]
    gamma         <- gamma / apply(gamma, 1, sum)          
  }                                                   
  delta <- solve(t(diag(m) - gamma + 1), rep(1, m))          
  list(mu = mu, sigma = sigma, gamma = gamma, delta = delta)           
}  

norm.HMM.mllk <- function(parvect, x, m, ...) {
  n        <- length(x)                            
  pn       <- norm.HMM.pw2pn(m, parvect)            
  P        <- rep(NA, m)
  lscale   <- 0                                    
  foo      <- pn$delta  
  for (i in 1:n) {
    for (j in 1:m) {
      P[j] <- dnorm(x[i], mean = pn$mu[j], sd = pn$sigma[j])
    }
    foo    <- foo %*% pn$gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo    <- foo / sumfoo
  }
  mllk     <- -lscale
  mllk
}

norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0,...) {                                                      
  parvect0 <- norm.HMM.pn2pw(m, mu0, sigma0, gamma0)         
  mod      <- nlm(norm.HMM.mllk, parvect0, x = x, m = m) 
  pn       <- norm.HMM.pw2pn(m, mod$estimate)            
  mllk     <- mod$minimum                              
  np       <- length(parvect0)
  AIC      <- 2 * (mllk + np)                              
  n        <- sum(!is.na(x))                            
  BIC      <- 2 * mllk + np * log(n)                         
  list(mu = pn$mu, sigma = pn$sigma, gamma = pn$gamma, delta = pn$delta, 
       code = mod$code, mllk = mllk, AIC = AIC, BIC = BIC)   
}

########################################################################################################################

## 2 - state 
## Initial values
m <- 2
mu0 <- mean(SPY)*c(1/2,3/2)
sigma0 <- sd(SPY)*c(1/2,3/2)
gamma0 <- matrix(0.01,ncol = m, nrow = m)
diag(gamma0) <- 1-(m - 1)*gamma0[1,1]

## optimize
fit2 <- norm.HMM.mle(SPY, m, mu0, sigma0, gamma0)
fit2


## 3 - state 
## Initial values
m <- 3
mu0 <- mean(SPY)*c(1/2,1,3/2)
sigma0 <- sd(SPY)*c(1/2,1,3/2)
gamma0 <- matrix(0.01,ncol = m, nrow = m)
diag(gamma0) <- 1-(m - 1)*gamma0[1,1]

## optimize
fit3 <- norm.HMM.mle(SPY, m, mu0, sigma0, gamma0)
fit3

## b) ##################################################################################################################

## Wald statistics for 3 state model 
m <- 3

## Working parameters
parvect  <- norm.HMM.pn2pw(m,fit3$mu, fit3$sigma, fit3$gamma)

## Optimize to get Hessian 
mod <- nlm(norm.HMM.mllk, parvect, x = SPY, m = m, hessian = TRUE)  
mod

## Organize the result
parvect <- mod$estimate
names(parvect) <- c("mu1","mu2","mu3",
                    "t.sigma1","t.sigma2","t.sigma3",
                    "tau21","tau31","tau12","tau32","tau13","tau23")

se <- sqrt(diag(solve(mod$hessian)))

## Working pars + standard error
round(cbind(parvect,se),digits=10) ## note se of tau31
fit3$gamma

## c) ##################################################################################################################



## d) ##################################################################################################################

## Compare the following distributions (by plots)

## The long term distribution of the return.

# This is probably the trickiest part of our work with Markov Chains. The stationary distribution, which we notate by  s
# , can (generally, usually, under certain conditions) be thought of as the long-run distribution of our chain. 
# What does that mean? Well, if we have a Markov Chain, and we let our random particle bounce around for a very long time, then  s
# describes the distribution of the location of the particle. Even more specifically, you can think of  s1
# as the proportion of time that the particle is in State 1 over the long run, and  s
# codes the proportion of time that the particle is in each of the states (not just State 1). We will formalize this later.

## e) ##################################################################################################################



## f) ##################################################################################################################



