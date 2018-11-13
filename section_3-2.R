#-----------------------------------------------------------------------
# Code for the simulation study in Section 3.2 of the paper
# "Variance prior forms for high-dimensional Bayesian variable selection"
# https://arxiv.org/abs/1801.03019
#-----------------------------------------------------------------------

library(MCMCpack)
library(ggplot2)
library(reshape2)

#-----------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------

# Gibbs sampler to compute posterior means for the independent prior 
gibbs_rr <- function(iter = 1000, burn = 500, Y, X, tXY, tXX, tau) {
  tot <- iter + burn
  p <- ncol(X)
  n <- nrow(X)
  b <- matrix(0, nrow = p, ncol = iter)
  s <- rep(0, iter)
  beta <- rnorm(0, p)
  sigma2 <- 1
  for(i in 1:tot) {
    V <- chol2inv(chol(tXX / sigma2 + diag(p) / tau^2))
    beta_mean <- 1/sigma2 * crossprod(V, tXY)
    beta <- mvrnorm(1, mu = beta_mean, Sigma = V)
    rss <- sum((Y - X%*%beta)^2)
    sigma2 <- rinvgamma(1, n/2, rss/2)
    if (i > burn) {
      b[, i - burn] <- beta
      s[i - burn] <- sigma2
    }
  }
  
  beta <- apply(b, 1, mean)
  sigma2 <- mean(s)
  
  return(list(beta = beta, sigma2 = sigma2))
}

#-----------------------------------------------------------------------
# Parameters
#
# n: number of observations
# p: number of covariates
# Y: response vector (n x 1)
# X: matrix of predictors (n x p)
# beta: true coefficient vector (p x 1)
# q: number of non-zero beta
# sd.true: true error variance
# tau: hyper-parameter for prior on beta
#-----------------------------------------------------------------------

n <- 100
p <- 90
q <- 6
tau <- 10
sd.true <- sqrt(3)

beta_vals <- c(-2.5, -2, -1.5, 1.5, 2, 2.5)
index <- sample(1:p, q)
beta <- rep(0, p)
beta[index] <- beta_vals

# store values
nrep <- 100
sigmas <- matrix(NA, nrow = nrep, ncol = 4)
betas <- array(NA, dim = c(nrep, p, 4))

#-----------------------------------------------------------------------
# Simulation Study
#-----------------------------------------------------------------------

for(i in 1:nrep) {
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- X%*%beta + rnorm(n, sd = sd.true)
  Y <- as.vector(Y)
  tXX <- crossprod(X)
  tXY <- crossprod(X, Y)
  
  # LS
  ls <- summary(lm(Y ~ X - 1))
  rmse <- ls$sigma
  ls_b <- ls$coefficients[,1]
  
  # Conjugate
  V_conj <- chol2inv(chol(tXX + diag(p) / tau^2))
  beta_conj <- V_conj %*% tXY
  H_conj <-  X %*% V_conj %*% t(X)
  sigma2_conj <- (t(Y) %*% (diag(n) - H_conj) %*% Y)/(n - 2)
  
  # Conjugate with Zellner's prior
  V_zell <- tau^2/(1 + tau^2) * chol2inv(chol(tXX))
  beta_zell <- V_zell %*% tXY
  H_zell <- X %*% V_zell %*% t(X)
  sigma2_zell <- (t(Y) %*% (diag(n) - H_zell) %*% Y)/(n - 2)
  
  # Independent prior
  ind_out <- gibbs_rr(iter = 1000, burn = 500, Y, X, tXY, tXX, tau)
  beta_ind <- ind_out$beta
  sigma2_ind <- ind_out$sigma2
  
  # Store values
  betas[i,, 1] <- ls_b
  betas[i,, 2] <- beta_conj
  betas[i,, 3] <- beta_zell
  betas[i,, 4] <- beta_ind
  sigmas[i,] <- c(rmse^2, sigma2_conj, sigma2_zell, sigma2_ind) 
  
  print(i) 
  
}  

stop()


# median of sigmas

apply(sigmas, 2, median)

# plot sigmas
colnames(sigmas) <- c("Least Squares", "Conjugate", "Zellner", "Independent")
sigmas_melt <- melt(sigmas)

g <- ggplot(sigmas_melt, aes(Var2, value)) 
g <- g + geom_boxplot(aes(fill = Var2))
g <- g + labs(x = " ", y = expression(hat(sigma)^2))
g <- g + scale_fill_discrete(labels = colnames(sigmas))
g <- g + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), 
               text = element_text(size = 18))
g <- g + theme(legend.title=element_blank()) + labs(x = NULL)
g <- g + geom_hline(yintercept = sd.true^2, color = "red")
g


# plot mean of betas

beta_means <- apply(betas, c(2, 3) , mean)
beta_means <- cbind(beta, beta_means)
colnames(beta_means) <- c("True", "Least Squares", "Conjugate", "Zellner", "Independent")

betas_melt <- melt(beta_means)
colnames(betas_melt) <- c("j", "Method", "beta_j")

g <- ggplot(betas_melt, aes(j, beta_j))
g <- g + labs(x = "j", y = expression(beta[j]))
g <- g + geom_jitter(width = 0.05, height = 0.05, mapping = aes(color = Method, shape = Method))
g

