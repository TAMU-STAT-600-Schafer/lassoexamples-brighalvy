# Generation of data for Example 3 (beta0 = 2, p = 2)
#############################################################
n = 50 # sample size
beta = c(1, 0.5) # true coefficients vector
beta0 = 2 # intercept, nonzero
p = length(beta)
sigma = 0.4 # noise standard deviation

library(mnormt)
set.seed(983645) # set seed

# Generate matrix of covariates
Sigma = matrix(0.7, p, p) + diag(rep(1-0.7, p)) # covariance for design X
X = rmnorm(n, mean = rep(0, p), varcov = Sigma) # matrix of covariates

# Generate response Y
Y = beta0 + X %*% beta + sigma * rnorm(n) 

# Centering and scaling
#############################################################
# [ToDo] Center both X and Y 
Xcenter <- t(t(X) - colMeans(X))
Ycenter <- Y - mean(Y)

# [ToDo] Scale centered X (so that n^{-1}X^tX has 1s on the diagonal)
normX <- diag(crossprod(Xcenter))/n
Xscaled <- Xcenter %*% diag(sqrt(1/normX))


# Coordinate descent implementation for the case p = 2
#############################################################
# Helper functions source
source("LassoFunctions.R")

niter = 50 # fixed number of iterations (for simplicity)
lambda = 0.5 # tuning parameter 
beta_start = rep(1, 2) # starting point 


# Apply coordinate descent on centered and scaled Xscaled, and centered Ycentered
out = coordinateLasso(Xscaled, Ycenter, beta_start, lambda, niter = niter)

plot(0:niter, out$fobj_vec)
# out2 = coordinateLasso2(Xscaled, Ycenter, beta_start, lambda, niter = niter)


# [ToDo] Perform back-centering and scaling to get intercept beta0 and vector of coefficients beta on the original scale
beta <- out$beta * sqrt(normX)

# Intercept according to the formula
beta0 <- mean(Y) - sum(colMeans(X)*beta)


