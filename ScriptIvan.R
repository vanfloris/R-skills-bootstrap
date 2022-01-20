####### Initialise the variables #######
# We look to simulate a VAR(2) with k = 4 equations
rm(list=ls())
set.seed(123) # Reset random number generator for reasons of reproducability
library(urca)
library(vars)
library(matlib)
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

## Create DeltaXt
series <- matrix(0, k, t + 2*p) # Raw series with zeros
Delta.Xt <- matrix(0, k, t + 2*p)
Xt.min1 <- matrix(0, k, t + 2*p)
Delta.Xt.min1 <- matrix(0, k, t + 2*p)
Delta.Xt.min2 <- matrix(0, k, t + 2*p+1)
Delta.Xt.min3 <- matrix(0, k, t + 2*p+2)
for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
  series[, i] <- A%*%series[, i-1] - Bmat%*%series[, i-2] + rnorm(k, 0, 1)
  Delta.Xt[,i] <- series[,i] - series[,i-1]
  Xt.min1[,i] <- series[,i-1]
  Delta.Xt.min1[,i] <- series[,i-1] - series[,i-2]
  Delta.Xt.min2[,i+1] <- series[,i-1] - series[,i-2]
  Delta.Xt.min3[,i+2] <- series[,i-1] - series[,i-2]
}
Delta.Xt.min2 <- Delta.Xt.min2[,1:54]
Delta.Xt.min3 <- Delta.Xt.min3[,1:54]
Delta.Xt <- t(Delta.Xt)
Xt.min1 <- t(Xt.min1)
Delta.Xt.min1 <- t(Delta.Xt.min1)
Delta.Xt.min2 <- t(Delta.Xt.min2)
Delta.Xt.min3 <- t(Delta.Xt.min3)

#ols <- function(Y,X.ols){ # OLS function ourselves #
#  b<- solve(crossprod(X.ols), crossprod(X.ols,Y)) # coefficient estimates
#  y.hat <- X.ols%*%b # fitted values
#  out <- list(coef.estimates = b, fitted.values = y.hat)
#  return(out)
#}
#ols(Delta.Xt[2:52,], Xt.min1[1:51,]- Delta.Xt[1:51,])

ols.lm1 <- lm(Delta.Xt~Xt.min1 + Delta.Xt.min1)
u <- ols.lm1$residuals

ols.lm2 <- lm(Xt.min1~Delta.Xt.min1+Delta.Xt.min2+Delta.Xt.min3)
v <- ols.lm2$residuals


for(i in  1:length(u[,1])){
  sigma.uu <- 1/length(u[,1]) * u[i,]%*%t(u[i,])
  sigma.vv <- 1/length(v[,1]) * v[i,]%*%t(v[i,])
  sigma.uv <- 1/length(v[,1]) * u[i,]%*%t(v[i,])
  sigma.vu <- 1/length(v[,1]) * v[i,]%*%t(u[i,])
}
Sigma <- solve(sigma.vv)%*%sigma.vu%*%solve(sigma.uu)%*%sigma.uv


###### Monte Carlo Simulation ######
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
  series <- matrix(0, k, t + 2*p) # Raw series with zeros
  for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
    series[, i] <- A%*%series[, i-1] - Bmat%*%series[, i-2] + rnorm(k, 0, 1)
  }
  X <- t(series)
  colnames(X) <- names
  
  ##Step 2: Apply Trace test ##
  ca <- ca.jo(X, type = "trace", K = 2, ecdet = "none")
  S <- summary(ca)
  teststats <- rev(S@teststat)
  Q[j,] <- teststats
  ## Step 3: Evaluate ##
  # Check if null hypothesis is rejected
  if (teststats[1] > 48.28) {reject.0[j] <- 1}
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


############# Estimated VAR based on simul  ############# 
# lag selection
lag <- VARselect(X, lag.max = 3)
lag$selection

## OLS
#ols <- function(Y,X.ols){ # OLS function ourselves #
#  b<- solve(crossprod(X.ols), crossprod(X.ols,Y)) # coefficient estimates
#  y.hat <- X.ols%*%b # fitted values
#  out <- list(coef.estimates = b, fitted.values = y.hat)
#  return(out)
#}
#ols(X[3:53,], X[2:52,]- X[1:51,])

# Estimate VAR
VARnew <- VAR(X, p = 2, type = "const")
res.VARnew <- residuals(VARnew)
sum <- summary(VARnew)
# mean matrix of residuals
mean.res <- cbind(mean(res.VARnew[,1]), mean(res.VARnew[,2]), mean(res.VARnew[,3]), mean(res.VARnew[,4])) 
mean.res1 <- mean(res.VARnew[,1]); mean.res2 <- mean(res.VARnew[,2]); mean.res3 <- mean(res.VARnew[,3]); mean.res4 <- mean(res.VARnew[,4])

# re-centered residuals
recenter.resids <- cbind(res.VARnew[,1] - mean.res1, res.VARnew[,2] - mean.res2, res.VARnew[,3] - mean.res3, res.VARnew[,4] - mean.res4) 

# Create the estimated 1st and 2nd lag matrices of the model and the estimated constant
lag1coef <- rbind(t(as.matrix(sum$varresult$V1$coefficients[1:4,1])),t(as.matrix(sum$varresult$V2$coefficients[1:4,1])),
                  t(as.matrix(sum$varresult$V3$coefficients[1:4,1])),t(as.matrix(sum$varresult$V4$coefficients[1:4,1])))
lag2coef <- rbind(t(as.matrix(sum$varresult$V1$coefficients[5:8,1])),t(as.matrix(sum$varresult$V2$coefficients[5:8,1])),
                  t(as.matrix(sum$varresult$V3$coefficients[5:8,1])),t(as.matrix(sum$varresult$V4$coefficients[5:8,1])))
const <- rbind((as.matrix(sum$varresult$V1$coefficients[9,1])),(as.matrix(sum$varresult$V2$coefficients[9,1])),
               (as.matrix(sum$varresult$V3$coefficients[9,1])), (as.matrix(sum$varresult$V1$coefficients[9,1])))

# Test if we can resample from estimated serie 
estseries <- matrix(0, k, t + 2*p) # Raw series with zeros
J <- sample.int(n, size = n, replace = TRUE) # Draw J
for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
  estseries[, i] <- lag1coef%*%estseries[, i-1] + lag2coef%*%estseries[, i-2] + const + recenter.resids[J[i],]
}
X.star <- t(estseries)
ts.plot(X.star)

##################### THE BOOTSTRAP IN R #####################
# First draw indices of the bootstrap sample: draw n times with replacement
n = 52
#J <- ceiling(runif(n, min = 0, max = n))
# we do B bootstrap replications and store the quantities in a vector
B = 399
Q.star1 <- matrix(data = NA,nrow= B, ncol = 4) 
reject.bstar.0 <- rep(0, times = B)
reject.bstar.1 <- rep(0, times = B)
for (b in 1:B) {
  J <- sample.int(n, size = n, replace = TRUE) # Draw J
  estseries1 <- matrix(0, k, t + 2*p) # Raw series with zeros
  for (i in (p + 1):(t + 2*p)){ # Generate series with recentered residuals
    estseries1[, i] <- lag1coef%*%estseries1[, i-1] + lag2coef%*%estseries1[, i-2] + const + recenter.resids[J[i],]
  }
  X.star <- t(estseries1)
  colnames(X.star) <- names
  ca.star <- ca.jo(X.star, type = "trace", K = 2, ecdet = "const")
  S.star <- summary(ca.star)
  teststats.star <- rev(S.star@teststat) #stored as teststat
  Q.star1[b,] <- teststats.star
  #if (teststats.star[1] > 48.28) {reject.bstar.0[b] <- 1}
  #if (teststats.star[2] > 31.52) {reject.bstar.1[b] <- 1}
}

cv.star1 <- quantile(Q.star1[,1], probs=0.95) ## Crit value for r = 0
cv.star2 <- quantile(Q.star1[,2], probs=0.95) ## Crit value for r = 1


######### Putting everything together #########
nr.sim <- 100; B <- 199;
n <- t + 2;
reject.star.0 <- rep(0, times = nr.sim)
reject.star.1 <- rep(0, times = nr.sim)
names <- c("V1", "V2", "V3", "V4") # Rename variables
Q.star1 <- matrix(data = NA,nrow= B, ncol = 4) 
for (i in 1:nr.sim){
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
    J <- sample.int(n, size = n, replace = TRUE) # Draw J
    estseries1 <- matrix(0, k, t + 2*p) # Raw series with zeros
    for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
      estseries1[, i] <- lag1coef%*%estseries1[, i-1] + lag2coef%*%estseries1[, i-2] + const + recenter.resids[J[i],]
    }
    X.star <- t(estseries1)
    colnames(X.star) <- names
    ca.star <- ca.jo(X.star, type = "trace", K = 2, ecdet = "none")
    S.star <- summary(ca.star)
    teststats.star <- rev(S.star@teststat) #stored as teststat
    Q.star1[b,] <- teststats.star
  }
  cv.star1 <- quantile(Q.star1[,2], probs=0.95)
  cv.star2 <- quantile(Q.star1[,3], probs=0.95)
  if (teststats[1] > cv.star1) {reject.star.0[b] <- 1}
  if (teststats[2] > cv.star2) {reject.star.1[b] <- 1}
}

## Step 4: Summarize ##
ERF.0 <- mean(reject.star.0)
ERF.1 <- mean(reject.star.1)
print(paste("Rejection occurred in ", 100 *ERF.0, "% of the cases."))
print(paste("Rejection occurred in ", 100 *ERF.1, "% of the cases."))  


### NOTES ###

# Check eigenvalues
#eigen(A.1)
#eigen(A.2)
#eigen(A) # Unstable eigenvalues
#eigen(B) # 1 eigenvalue is equal to 1


# Generate series
#series <- matrix(0, k, t + 2*p) # Raw series with zeros
#for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,1)
#  series[, i] <- A%*%series[, i-1] - B%*%series[, i-2] + rnorm(k, 0, 1)
#}
#seriests <- ts(t(series[, -(1:p)])) # Convert to time series format
#ts.plot(seriests)
#transseries <- t(series)
#names <- c("V1", "V2", "V3", "V4") # Rename variables
#colnames(transseries) <- names

#ca <- ca.jo(transseries, type = "trace", K = 2, ecdet = "none")
#S <- summary(ca)
#cv <- c(48.28, 31.52, 17.95, 8.18)
#teststats <- rev(S@teststat)
#for(i in 1:4){
#  if(teststats[i]< cv[i]){
#    return(c.rank <- i-1) 
#  }
#}

