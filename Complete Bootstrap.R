####### Bootstrap for the ERF of VAR(2) model#######
# We look to simulate a VAR(2) with k = 4 equations
rm(list=ls())
set.seed(123) # Reset random number generator for reasons of reproducability
library(urca)
library(vars)

# Generate sample
t <- 50 # Number of time series observations
n <- t + 2
k <- 4 # Number of endogenous variables
p <- 2 # Number of lags

# Generate coefficient matrices
a <- -0.4 ; gamma <- 0.8 ; alpha <- t(t(c(a,0,0,0))); beta <- t(t(c(1,0,0,0)))
delta <- 0; # 0 , 0.1 , 0.2 , 0.3. 
A.1 <- alpha %*% t(beta) # Alpha matrix
A.2 <- diag(x = 1, k) # 4x4 identity matrix
Bmat <- matrix(c(gamma, delta, 0, 0, delta, gamma, 0, 0, 0, 0, gamma, 0, 0, 0, 0, gamma), k) # Gamma matrix
A <- A.1 + A.2 + Bmat

###### Monte Carlo Simulation to get Q_0,T and Q_1,T ######
# Number of simulations
nr.sim <- 10000
# Initialize a vector of 0s to store rejections
reject.0 <- rep(0, times = nr.sim)
reject.1 <- rep(0, times = nr.sim)
names <- c("V1", "V2", "V3", "V4") # Rename variables
Q <- matrix(data = NA,nrow= nr.sim, ncol = 4)
for (j in 1:nr.sim){
  ## Step 1: Simulate ##
  # Generate sample from VAR
  seriessim <- matrix(0, k, t + 2*p) # Raw series with zeros
  for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
    seriessim[, i] <- A%*%seriessim[, i-1] - Bmat%*%seriessim[, i-2] + rnorm(k, 0, 1)
  }
  X <- t(seriessim)
  colnames(X) <- names
  
  ##Step 2: Apply Trace test ##
  ca <- ca.jo(X, type = "trace", K = 2, ecdet = "none")
  S <- summary(ca)
  teststats <- rev(S@teststat)
  Q[j,] <- teststats
  ## Step 3: Evaluate ##
  # Check if null hypothesis is rejected
  if (teststats[1] > 71.1) {reject.0[j] <- 1}
  if (teststats[2] > 31.52) {reject.1[j] <- 1}
}

## Step 4: Summarize ##
# Empirical rejection frequency of rank = 0
ERF.0 <- mean(reject.0)
# Empirical rejection frequency of rank = 1
ERF.1 <- mean(reject.1)
# give the output on the screen
print(paste("Chance to reject 0: ", ERF.0))
print(paste("Chance to reject 1: ", ERF.1))

### Create all necessary things for the bootstrap ###
## Create Delta.Xt for OLS
series <- matrix(0, k, t + 2*p) # Raw series with zeros
Delta.Xt <- matrix(0, k, t + 2*p)
names <- names <- c("V1", "V2", "V3", "V4") # Rename variables
Xt.min1 <- matrix(0, k, t + 2*p)
Delta.Xt.min1 <- matrix(0, k, t + 2*p)
for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
  series[, i] <- A%*%series[, i-1] - Bmat%*%series[, i-2] + rnorm(k, 0, 1)
  Delta.Xt[,i] <- series[,i] - series[,i-1]
  Xt.min1[,i] <- series[,i-1]
  Delta.Xt.min1[,i] <- series[,i-1] - series[,i-2]
}
series <- t(series)
colnames(series) <- names
Delta.Xt <- t(Delta.Xt)
Xt.min1 <- t(Xt.min1)
Delta.Xt.min1 <- t(Delta.Xt.min1)

### OLS as step 1 in Hamilton
ols.lm1 <- lm(Delta.Xt~ Delta.Xt.min1) 
u <- ols.lm1$residuals
ols.lm2 <- lm(Xt.min1~Delta.Xt.min1)
v <- ols.lm2$residuals

### Step 2 in Hamilton
sigma.uu <- matrix(0, 4, 4)
sigma.vv <- matrix(0, 4, 4)
sigma.uv <- matrix(0, 4, 4)
sigma.vu <- matrix(0, 4, 4)
for(i in  1:length(v[,1])){
  sigma.uu <- sigma.uu + 1/length(v[,1])*(u[i,]%*%t(u[i,])) 
  sigma.vv <- sigma.vv + 1/length(v[,1])*(v[i,]%*%t(v[i,]))
  sigma.uv <- sigma.uv + 1/length(v[,1])*(u[i,]%*%t(v[i,]))
  sigma.vu <- sigma.vu + 1/length(v[,1])*(v[i,]%*%t(u[i,]))
}
Sigma <- solve(sigma.vv)%*%sigma.vu%*%solve(sigma.uu)%*%sigma.uv

### Step 3 in Hamilton
eigensigma <- eigen(Sigma)
eig.vect1 <- as.matrix(eigensigma$vectors[,1])

# create a normalized eigenvector for r = 1
norm.vec <- (1/(sqrt(t(eig.vect1)%*% sigma.vv %*% eig.vect1)))[1,1] * eig.vect1
# when sigma.vv is pre- and post-multiplied by this vector, the result is 1.
# We can then use this normalized vector to obtain the parameter estimates.

beta.hat.r0 <- t(cbind(0,0,0,0)) # rank of 0
beta.hat.r1 <- rev(eigensigma$vector[,1]) # rank of 1
nbeta.hat.r1 <- eig.vect1 # normalized eigenvector of rank 1
tbeta.hat.r1 <- t(beta.hat.r1)

## estimation of coefficients following the paper's notation for r = 0 ##
zeta.hat.0.r0 <- sigma.uv %*% beta.hat.r0 %*% t(beta.hat.r0) # alpha beta' of our paper
gamma.hat.1.r0 <- ols.lm1$coefficients[2:5,] - zeta.hat.0.r0 %*% ols.lm2$coefficients[2:5,] # pi.1 - zeta.hat.0 * chi.1
alpha.hat.r0 <- sigma.uv %*% beta.hat.r0 #alpha-hat as Sigma.UV A.hat (as zeta_0 in Ham corresponds to alpha beta' in the paper)

## estimation of coefficients following the paper's notation for r = 1 ##
zeta.hat.0.r1 <- sigma.uv %*% norm.vec %*% t(norm.vec) # alpha beta' of our paper
gamma.hat.1.r1 <- ols.lm1$coefficients[2:5,] - zeta.hat.0.r1 %*% ols.lm2$coefficients[2:5,] # pi.1 - zeta.hat.0 * chi.1
alpha.hat.r1 <- sigma.uv %*% nbeta.hat.r1 #alpha-hat as Sigma.UV A.hat (as zeta_0 in Ham corresponds to alpha beta' in the paper)

# Estimate VAR to get the recentered residuals
VARnew <- VAR(series, p = 2, type = "const")
res.VARnew <- residuals(VARnew)
sum <- summary(VARnew)
# mean matrix of residuals
mean.res <- cbind(mean(res.VARnew[,1]), mean(res.VARnew[,2]), mean(res.VARnew[,3]), mean(res.VARnew[,4])) 
mean.res1 <- mean(res.VARnew[,1]); mean.res2 <- mean(res.VARnew[,2]); mean.res3 <- mean(res.VARnew[,3]); mean.res4 <- mean(res.VARnew[,4])

# re-centered residuals
recenter.resids <- cbind(res.VARnew[,1] - mean.res1, res.VARnew[,2] - mean.res2, res.VARnew[,3] - mean.res3, res.VARnew[,4] - mean.res4) 

B = 99
Q.star1 <- matrix(data = NA,nrow= B, ncol = 4) 
reject.bstar.0 <- rep(0, times = B)
reject.bstar.1 <- rep(0, times = B)
for (b in 1:B) {
  est.VAR <- matrix(0, k, t + 2*p) # Raw series with zeros
  coef1 <- zeta.hat.0.r0 + gamma.hat.1.r0 + A.2
  J <- sample.int(n, size = n, replace = TRUE) # Draw J
  for (i in (p + 2):(t + 2*p)){ # Generate estimated series with recentered residuals
    est.VAR[, i] <- coef1%*%est.VAR[, i-1] - gamma.hat.1.r0%*%est.VAR[,i-2] + recenter.resids[J[i],] # formula 8 of paper
  }
  X.star <- t(est.VAR)
  colnames(X.star) <- names
  ca.star <- ca.jo(X.star, type = "trace", K = 2, ecdet = "const")
  S.star <- summary(ca.star)
  teststats.star <- rev(S.star@teststat) #stored as teststat
  Q.star1[b,] <- teststats.star
}
cv.star1 <- quantile(Q.star1[,1], probs=0.95) ## Crit value for r = 0
cv.star1

######### Bootstrap r = 0, to get Q.star_0,T#########
nr.sim <- 500; B <- 99;
n <- t + 2;
reject.star.0 <- rep(0, times = nr.sim)
names <- c("V1", "V2", "V3", "V4") # Rename variables
Q.star0 <- matrix(data = NA,nrow= B, ncol = 4) 
for (j in 1:nr.sim){
  ## Step 1: Simulate ##
  series <- matrix(0, k, t + 2*p) # Raw series with zeros
  for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
    series[, i] <- A%*%series[, i-1] - Bmat%*%series[, i-2] + rnorm(k, 0, 1)
  }
  X <- t(series)
  colnames(X) <- names
  ## Step 2: Apply ##
  ca <- ca.jo(X, type = "trace", K = 2, ecdet = "none")
  S <- summary(ca)
  teststats <- rev(S@teststat)
  for (b in 1:B) {
    est.VAR <- matrix(0, k, t + 2*p) # Raw series with zeros
    coef0 <- zeta.hat.0.r0 + gamma.hat.1.r0 + A.2
    J <- sample.int(n, size = n, replace = TRUE) # Draw J
    for (i in (p + 2):(t + 2*p)){ # Generate estimated series with recentered residuals
      est.VAR[, i] <- coef0%*%est.VAR[, i-1] - gamma.hat.1.r0%*%est.VAR[,i-2] + recenter.resids[J[i],] # formula 8 of paper
    }
    X.star <- t(est.VAR)
    colnames(X.star) <- names
    ca.star <- ca.jo(X.star, type = "trace", K = 2, ecdet = "const")
    S.star <- summary(ca.star)
    teststats.star <- rev(S.star@teststat) #stored as teststat
    Q.star0[b,] <- teststats.star
  }
  cv.star1 <- quantile(Q.star0[,1], probs=0.95)
  if (teststats[1] > cv.star1) {reject.star.0[j] <- 1}
}

## Step 4: Summarize ##
ERF.0 <- mean(reject.star.0)
print(paste("Rejection occurred in ", 100 *ERF.0, "% of the cases.")) 

######### Bootstrap r = 1, to get Q.star_1,T #########
nr.sim <- 500; B <- 99;
n <- t + 2;
reject.star.1 <- rep(0, times = nr.sim)
names <- c("V1", "V2", "V3", "V4") # Rename variables
Q.star1 <- matrix(data = NA,nrow= B, ncol = 4) 
for (j in 1:nr.sim){
  ## Step 1: Simulate ##
  series <- matrix(0, k, t + 2*p) # Raw series with zeros
  for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
    series[, i] <- A%*%series[, i-1] - Bmat%*%series[, i-2] + rnorm(k, 0, 1)
  }
  X <- t(series)
  colnames(X) <- names
  ## Step 2: Apply ##
  ca <- ca.jo(X, type = "trace", K = 2, ecdet = "none")
  S <- summary(ca)
  teststats <- rev(S@teststat)
  for (b in 1:B) {
    est.VAR <- matrix(0, k, t + 2*p) # Raw series with zeros
    coef1 <- zeta.hat.0.r1 + gamma.hat.1.r1 + A.2
    J <- sample.int(n, size = n, replace = TRUE) # Draw J
    for (i in (p + 2):(t + 2*p)){ # Generate estimated series with recentered residuals
      est.VAR[, i] <- coef1%*%est.VAR[, i-1] - gamma.hat.1.r1%*%est.VAR[,i-2] + recenter.resids[J[i],] # formula 8 of paper
    }
    X.star <- t(est.VAR)
    colnames(X.star) <- names
    ca.star <- ca.jo(X.star, type = "trace", K = 2, ecdet = "const")
    S.star <- summary(ca.star)
    teststats.star <- rev(S.star@teststat) #stored as teststat
    Q.star1[b,] <- teststats.star
  }
  cv.star2 <- quantile(Q.star1[,2], probs=0.95)
  if (teststats[2] > cv.star2) {reject.star.1[j] <- 1}
}

## Step 4: Summarize ##
ERF.1 <- mean(reject.star.1)
print(paste("Rejection occurred in ", 100 *ERF.1, "% of the cases.")) 

