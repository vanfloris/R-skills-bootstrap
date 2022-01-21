#Merge test into master
# Merge test MIchiel
rm(list=ls())
?solve
# Ordinary Least Squares
X <- matrix(NA, 3, 3)
Y <- matrix(NA, 3, 3)
n <- c(NA)

beta.hat <- solve(crossprod(X),crossprod(X,Y)) # least squares estimator

e.hat <- Y - X%*%beta.hat # residuals
s2 <- crossprod(e.hat)/n # sigma hat squared

sigma.hat <- s2 * solve(crossprod(X)) # covariance matrix

t.stat <- (beta.hat[1] - 1)/sqrt(sigma.hat[1,1]) # test whether first element of beta is 1 over std. error (SD)

# Functions

ols <- function(y,x) {
  b <- solve(crossprod(X), crossprod(X,Y))
  return(b)
}
beta.hat <- ols(Y,X) # same OLS estimator, but can call back to it now

# can have multiple outputs for functions as a list
ols <- function(Y,X){
  b<- solve(crossprod(X), crossprod(X,Y)) # coefficient estimates
  y.hat <- x%*%b # fitted values
  out <- list(coef.estimates = b, fitted.values = y.hat)
  return(out)
}
# first value spits out the coefficient estimate, the second value spits out 
ols.out <- ols(Y,X)
beta.hat <- ols.out[[1]]



############ MONTE CARLO SIMULATION ############ 

# generate pseudo random variables (can be rnorm or runif, or...)
Z <- rnorm(100, mean = 10, sd = 5) # vector of 100 normally distributed variables with mean 10 and sd 5

# generate a correlated series, use Choleski decomposition where L'L is covariance matrix
Y <- matrix(rnorm(200), ncol = 2)
Sigma <- matrix(c(1,0.4,0.4,1), nrow = 2)
L <- chol(Sigma)
X <- Y%*%L # x is now correlated with the covariance matrix


# Generate an AR process
eps <- rnorm(100) # epsilon
X <- filter(x = eps, filter = c(phi.1, phi.2), method = "r", init = c(X.0, X.m1)) # with parameters, method is
#recursive, and initial values in reverse order
# to reduce impact of the initial observations, generate n+k variables and discard the first k. (burn in)


# Generate VAR process

phi <- matrix(c(0.5,0.1,0.2,0.4), nrow = 2)
eps <- matrix(rnorm(200), ncol = 2) # given n is 100
X <- matrix(nrow = 100, ncol = 2)
X[1,] <- eps[1,] # first starting point for "for" loop
for (t in 2:100) {
  X[t,] <- X[t-1,]%*%t(phi) + eps[t,]
}



####### Performing a Hypothesis Test #######

# set the seed for the random number generator
set.seed(515)
# Number of simulations
nr.sim <- 5000
# Sample Size
n <- 100
# Nominal level of the test
alpha <- 0.05
# set the true value of the mean: <=0 for size, >0 for power
mu <- 0
# Initialize a vector of 0s to store rejections
reject.n <- rep(0, times = nr.sim)
# Initialize a vector of 0s to store rejections
reject.t <- rep(0, times = nr.sim)

###### Start the Simulation ######
for (i in 1:nr.sim){
  ## Step 1: Simulate ##
  X <- rnorm(n, mean = mu) # Draw X
  ##Step 2: Apply ##
  X.bar <- mean(X) # Sample Mean
  St.Dev <- sd(X) # Sample Std. Dev
  Q.n <- sqrt(n) * X.bar / St.Dev # Test statistic
  cv.n <- qnorm(1-alpha) # normal critical value
  cv.t <- qt(1-alpha, n-1) # t critical value
  ## Step 3: Evaluate ##
  # Check if null hypothesis is rejected
  if (Q.n > cv.n) {reject.n[i] <- 1}
  if (Q.n > cv.t) {reject.t[i] <- 1}
}

## Step 4: Summarize ##
# Empirical rejection frequency (normal CV)
ERF.n <- mean(reject.n)
# Empirical rejection frequency (t CV)
ERF.t <- mean(reject.t)
# give the output on the screen
if (mu == 0) {
  print(paste("Size using normal CV: ", ERF.n))
  print(paste("Size using t CV: ", ERF.t))
} else if (mu>0) {
  print(paste("Power using normal CV: ", ERF.n))
  print(paste("Power using t CV: ", ERF.t))
}



##################### THE BOOTSTRAP IN R #####################

# First draw indices of the bootstrap sample: draw n times with replacement
n = 100
J <- ceiling(runif(n, min = 0, max = n))
X.star <- X[J] # put elements of X corresponding to the drawn indices in a vector X*
# we do B bootstrap replications and store the quantities in a vector
B = 499
Q.star <- rep(NA, times = B);
for (b in 1:B){
  J <- sample.int(n, size = n, replace = TRUE)
  X.star <- X[J]
  Q.star[b] <- Q.func(X.star, theta.star) # function and theta star need to be set prior (not here)
}

#bootstrap testing for the mean
set.seed(515)
b <- 499 # number of bootstrap replications
alpha <- 0.05 # Nominal level of the test
mu <- 0 # true value of the mean
X <- rnorm(n, mean = mu) # simulated data stored as X
n <- length(X) # sample size
X.bar <- mean (X) # sample mean
St.dev <- sd(X) # standard deviation
Q.n <- sqrt(n)*X.bar/St.dev # Test statistic

# now use bootstrap to find the critical value
# initialize vector for bootstrap statistics
Q.star <- rep(NA, times = B)
for (b in 1:B) {
  J <- sample.int(n, size = n, replace = TRUE) # Draw J
  X.star <- X[J] # Draw bootstrap sample
  X.bar.star <- mean(X.star) # bootstrap sample mean
  St.dev.star <- sd(X.star) # bootstrap standard deviation
  # bootstrap statistic must be stored
  Q.star[b] <- sqrt(n)*(X.bar.star-X.bar)/St.dev.star
}

# calculate the bootstrap critical value
cv <- quantile(Q.star, probs = 1-alpha)
p.val <- mean(Q.star>Q.n)
# p value of 0.06212 means we fail to reject the null hypothesis


######### Testing for the Mean #########
set.seed(515)
nr.sim <- 2000; B <- 499;
n <- 100; alpha <- 0.05;
mu <- 0 # true value of mean
reject <- rep(0, times = nr.sim) # vector that will store my rejections

for (i in 1:nr.sim){
  ## Step 1: Simulate ##
  X <- rnorm(n, mean = mu)
  ## Step 2: Apply ##
  X.bar <- mean(X)
  St.Dev <- sd(X)
  Q <- sqrt(n)*X.bar/St.Dev
  Q.s <- rep(NA, times = B)
  for (b in 1:B){ # nested bootstrap step
    J <- sample.int(n, size = n, replace = TRUE)
    X.star <- X[J]
    X.bar.star <- mean(X.star)
    St.Dev.star <- sd(X.star)
    Q.s[b] <- sqrt(n)*(X.bar.star - X.bar)/St.Dev.star
  }
  cv.star <- quantile(Q.s, probs = 1-alpha)
  ## Step 3: Evaluate ##
  if (Q>cv.star) {reject[i] <- 1}
}
## Step 4: Summarize ##
ERF <- mean(reject)
print(paste("Rejection occurred in ", 100 *ERF, "% of the cases."))

# rejection occurred in 6.15% of the cases. 


