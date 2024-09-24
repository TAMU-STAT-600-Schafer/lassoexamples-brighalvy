#############################################################
# Generation of data 2 for Example 1 (beta0 = 0, p = 2)
#############################################################
n = 50 # sample size
beta = c(1, 0.5) # true coefficients vector
p = length(beta)
sigma = 0.4 # noise standard deviation

library(mnormt)
set.seed(234865) # set seed

# Generate matrix of covariates
Sigma = matrix(0.7, p, p) + diag(rep(1-0.7, p)) # covariance for design X
X = rmnorm(n, mean = rep(0, p), varcov = Sigma) # matrix of covariates

# Generate response Y
Y = X %*% beta + sigma * rnorm(n) 

# Verify has the right scaling
norm1 = as.numeric(crossprod(X[ , 1])/n) # should be 1
norm2 = as.numeric(crossprod(X[ , 2])/n) # should be 1

# Scaling of matrix of covariates X so that diag(X'X)/n has 1s
############################################################
# [ToDo] Create Xscaled from X that satisfies the above
Xscaled <- X
Xscaled[,1] <- X[,1]/sqrt(norm1)
Xscaled[,2] <- X[,2]/sqrt(norm2)

apply(Xscaled, 2, \(col) crossprod(col)/n)


# Coordinate descent implementation for the case p = 2, same as before but with scaled X
#############################################################
# Helper functions source
source("LassoFunctions.R")

niter = 50 # fixed number of iterations (for simplicity)
lambda = 0.2 # tuning parameter 
beta_start = rep(0, 2) # starting point 

# Apply coordinate descent
out = coordinateLasso(Xscaled, Y, beta_start, lambda, niter = niter)

plot(0:niter, out$fobj_vec)

# [ToDo] Transform out$beta so it can be applied to original X rather than Xscaled
beta_unscaled <- out$beta*c(sqrt(norm1), sqrt(norm2))

