#############################################################
# Generation of data 1 for Example 1 (beta0 = 0, p = 2, X is scaled)
#############################################################
n = 48 # sample size
beta = c(1, 0.5) # true coefficients vector
p = length(beta)
sigma = 0.4 # noise standard deviation

library(mnormt)
set.seed(234865) # set seed

# Generate matrix of covariates
X1 = c(rep(1, n/2), rep(-1, n/2))
X2 = c(rep(2, 6), rep(-2, 6), rep(0, n - 12))
# Verify has the right scaling
crossprod(X1)/n # should be 1
crossprod(X2)/n # should be 1
X = cbind(X1, X2)

# Generate response Y
Y = X %*% beta + sigma * rnorm(n) 

# Coordinate descent implementation for the case p = 2
#############################################################
# Helper functions source
source("LassoFunctions.R")

niter = 50 # fixed number of iterations (for simplicity)
lambda = 0 # tuning parameter 
beta_start = rep(0, 2) # starting point 

# Create storage vectors:
fvec <- vector(length = niter + 1)
fvec[1] <- lassoobj(X, Y, beta_start, lambda)

beta <- beta_start
# For loop
for(i in 1:niter){
  # update beta1
  x1 <- (1/n) * as.numeric(crossprod(X[,1], Y - X[,2]*beta[2]))
  beta[1] <- softthresh(x1, lambda)
  # update beta2
  x2 <- (1/n) * as.numeric(crossprod(X[,2], Y - X[,1]*beta[1]))
  beta[2] <- softthresh(x2, lambda)
  
  # Calculate objective
  fvec[i + 1] <- lassoobj(X, Y, beta, lambda)
}

beta1 <- beta

# [ToDo] Feel free to see how the results change with different set of parameters below
###################################################################################

# Second set of parameters
lambda = 0.5 # tuning parameter 1
beta_start = rep(2, 2) # starting point 2

# Create storage vectors:
fvec2 <- vector(length = niter + 1)
fvec2[1] <- lassoobj(X, Y, beta_start, lambda)

beta <- beta_start
# For loop
for(i in 1:niter){
  # update beta1
  x1 <- (1/n) * as.numeric(crossprod(X[,1], Y - X[,2]*beta[2]))
  beta[1] <- softthresh(x1, lambda)
  # update beta2
  x2 <- (1/n) * as.numeric(crossprod(X[,2], Y - X[,1]*beta[1]))
  beta[2] <- softthresh(x2, lambda)
  
  # Calculate objective
  fvec2[i + 1] <- lassoobj(X, Y, beta, lambda)
}

beta2 <- beta

# Third set of parameters
lambda = 1 # tuning parameter 2
beta_start = rep(0, 2) # starting point 1

# Create storage vectors:
fvec3 <- vector(length = niter + 1)
fvec3[1] <- lassoobj(X, Y, beta_start, lambda)

beta <- beta_start
# For loop
for(i in 1:niter){
  # update beta1
  x1 <- (1/n) * as.numeric(crossprod(X[,1], Y - X[,2]*beta[2]))
  beta[1] <- softthresh(x1, lambda)
  # update beta2
  x2 <- (1/n) * as.numeric(crossprod(X[,2], Y - X[,1]*beta[1]))
  beta[2] <- softthresh(x2, lambda)
  
  # Calculate objective
  fvec3[i + 1] <- lassoobj(X, Y, beta, lambda)
}

beta3 <- beta

# [ToDo] Feel free to see how the results change with different set of parameters below
###################################################################################

# Second set of parameters
lambda = 0.5 # tuning parameter 
beta_start = rep(2, 2) # starting point 

# Apply coordinate descent
out = coordinateLasso(X, Y, beta_start, lambda, niter = niter)

plot(0:niter, out$fobj_vec)

# Third set of parameters
lambda = 1 # tuning parameter 
beta_start = rep(0, 2) # starting point 1

# Apply coordinate descent
out = coordinateLasso(X, Y, beta_start, lambda, niter = niter)

plot(0:niter, out$fobj_vec)
