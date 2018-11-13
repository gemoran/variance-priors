#-----------------------------------------------------------------------
# Code for the analysis in Section 7 of the paper
# "Variance prior forms for high-dimensional Bayesian variable selection"
# https://arxiv.org/abs/1801.03019
#-----------------------------------------------------------------------

library(BAS)
library(SSLASSO)
library(ncvreg)
library(parcor)
library(scalreg)
library(xtable)
library(ggplot2)
library(reshape2)

set.seed(12345)

# function to compute Gaussian loss
gLoss = function(y, X, beta) {
  return(sum((y - X %*% beta)^2))
}

# --------------- Methods compared --------------- #
# SSL 
# MCP
# SCAD
# LASSO
# Scaled Lasso
# MCP hard (gamma = 1.0001)

#---------------Parameters --------------- #
# p: number of predictors
# n: number of observations
# y: response, n x 1 vector
# X: predictors, n x p matrix

#------------------------------------------
# Get data
#------------------------------------------

# load protein dataset from the BAS package
data(protein)

# get design matrix with interactions
protein.lm <- lm(prot.act4 ~
                   buf + buf*pH + buf*NaCl + buf*con + buf*ra +
                   buf*det + buf*MgCl2 + buf*temp +
                   pH + pH*NaCl + pH*con + pH*ra + pH*det +
                   pH*MgCl2 +pH*temp +
                   NaCl + NaCl*con + NaCl*ra + NaCl*det +
                   NaCl*MgCl2 + NaCl*temp +
                   con + con*ra + con*det +
                   con*MgCl2 +con*temp +
                   ra + ra*det + ra*MgCl2 + ra*temp +
                   det + det*MgCl2 + det*temp +
                   MgCl2 + MgCl2*temp + I(NaCl^2) + I(pH^2) + I(con^2) + I(temp^2),
                 data=protein,x=T)

protein.designmat <- protein.lm$x
proteindata <-  data.frame(cbind(protein.designmat[,-1],protein$prot.act4))
names(proteindata)[89] <- "prot.act4"

y <- proteindata$prot.act4
X <- protein.designmat[,-1]
# save variable names as X_names
X_names <- colnames(protein.designmat)[-1]

colnames(X) <- paste("X",1:88,sep="")

p <- ncol(X)
n <- nrow(X)

#------------------------------------------
# BAS results (benchmark)
#------------------------------------------

# Run BAS
# NOTE: this takes a while
bas_fit <- bas.lm(y~ X,
                    prior="JZS",
                    method = "BAS")


# Look at Posterior Inclusion Probabilities (PIP) from BAS 
pip <- summary(bas_fit)[c(2:(p+1)),1]

# Median Probability Model: variables with PIP > 0.5
bas_index <- which(pip > 0.5)

# plot PIP
plot(1:p, pip, pch = 2, col = "blue")

#------------------------------------------
# Spike-and-Slab LASSO with unknown variance
#------------------------------------------

# run SSL (unknown var)
fit <- SSLASSO(X, y, variance = "unknown")

# how many variables found by BAS did SSL (unknown var) find?
sum(fit$model %in% bas_index)

# which variables did SSL (unknown var) find?
fit$model
X_names[fit$model]

SSL_not_BAS <- fit$model[!(fit$model %in% bas_index)]
BAS_not_SSL <- bas_index[!(bas_index %in% fit$model)]

# predictor correlation matrix
X_cor <- cor(X)
X_cor[SSL_not_BAS, BAS_not_SSL]

# predictors 54 and 11 highly correlated (rho = 0.988)

# sigma estimate from SSL (unknown var)
fit$sigmas[length(fit$sigmas)]


#------------------------------------------
# Spike-and-Slab LASSO with fixed variance
#------------------------------------------

# get initial value for sigma
df <- 3
sigquant <- 0.9
sigest <- sd(y)
qchi <- qchisq(1 - sigquant, df)
ncp <- sigest^2 * qchi / df
sigma_est <- sqrt(df * ncp / (df - 2))

# run SSL (fixed var)
fit <- SSLASSO(X, y, variance = "fixed", sigma = sigma_est)

# how many variables found by BAS did SSL (unknown var) find?
sum(fit$model %in% bas_index)

# which variables did SSL (unknown var) find?
fit$model
X_names[fit$model]

SSL_not_BAS <- fit$model[!(fit$model %in% bas_index)]
BAS_not_SSL <- bas_index[!(bas_index %in% fit$model)]

# predictor correlation matrix
X_cor[SSL_not_BAS, BAS_not_SSL, drop = F]

# predictors 68 and 10 are highly correlated (rho = 0.735)


#------------------------------------------
# 8-fold Cross validation
#------------------------------------------

nfolds <- 8
nrep <- 100 
methods <- c("varSSL",
             "SSL", 
             "MCP",
             "SCAD",
             "LASSO",
             "scaledLASSO",
             "hardMCP",
             "adaLASSO")

# store results
result <- array(0, dim = c(length(methods), nfolds, nrep))
rownames(result) <- methods


# simulation study
for(i in 1:nrep) {
  cv.ind <- ceiling(sample(1:n)/(n+sqrt(.Machine$double.eps))*nfolds)
  
  for(k in 1:nfolds) {
    
    y.train <- y[cv.ind != k]
    X.train <- X[cv.ind != k, ]
    y.test <- y[cv.ind == k]
    X.test <- X[cv.ind == k, ]
    
    # Spike-and-Slab LASSO with unknown variance
    
    time <- system.time(fit <- SSLASSO(X.train, y.train, variance = "unknown"))
    
    beta_est <- fit$beta[, ncol(fit$beta)]
    result["varSSL", k, i] <- gLoss(y.test, X.test, beta_est)
    
    # Spike-and-Slab LASSO with fixed variance
    # get initial value for sigma
    df <- 3
    sigquant <- 0.9
    sigest <- sd(y.train)
    qchi <- qchisq(1 - sigquant, df)
    ncp <- sigest^2 * qchi / df
    sigma_est <- sqrt(df * ncp / (df - 2))
    
    
    time <- system.time(fit <- SSLASSO(X.train, y.train, variance = "fixed", sigma = sigma_est))
    
    beta_est <- fit$beta[, ncol(fit$beta)]
    result["SSL", k, i] <- gLoss(y.test, X.test, beta_est)
    
    
    # MCP
    
    time <- system.time(cvfit <- cv.ncvreg(X.train, y.train, family = c("gaussian"), penalty = "MCP", warn = FALSE))
    
    fit <- cvfit$fit
    beta_est <- fit$beta[-1, cvfit$min]
    result["MCP", k, i] <- gLoss(y.test, X.test, beta_est)
    
    
    # SCAD
    
    time <- system.time(cvfit <- cv.ncvreg(X.train, y.train, family = c("gaussian"), penalty = "SCAD", warn = FALSE))
    
    fit <- cvfit$fit
    beta_est <- fit$beta[-1, cvfit$min]
    result["SCAD", k, i] <- gLoss(y.test, X.test, beta_est)
    
    # LASSO
    
    time <- system.time(cvfit <- cv.ncvreg(X.train, y.train, family=c("gaussian"),penalty="lasso",warn=FALSE))
    
    fit <- cvfit$fit
    beta_est <- fit$beta[-1, cvfit$min]
    result["LASSO", k, i] <- gLoss(y.test, X.test, beta_est)
    
    # Scaled  LASSO
    
    time <- system.time(fit <- scalreg(X.train, y.train))
    
    beta_est <- fit$coefficients
    result["scaledLASSO", k, i] <- gLoss(y.test, X.test, beta_est)
    
    
    # MCP  (gamma = 1.0001)
    
    time <- system.time(cvfit <- cv.ncvreg(X.train, y.train, family=c("gaussian"),penalty="MCP",warn=FALSE,gamma=1.0001))
    
    fit <- cvfit$fit
    
    beta_est <- fit$beta[-1, cvfit$min]
    result["hardMCP", k, i] <- gLoss(y.test, X.test, beta_est)
    
    
    # Adaptive lasso
    
    time <- system.time(fit <- adalasso(X.train, y.train, k = 10, intercept=F))
    
    beta_est <- fit$coefficients.adalasso
    result["adaLASSO", k, i] <- gLoss(y.test, X.test, beta_est)
    
  } 
  
  print(i)
}

# calculate 8-fold cross-validation error for every replication
result_cv <- apply(result, c(1,3), mean)

# outliers in scaled LASSO
sum(result_cv["scaledLASSO",] > 100)

# remove scaled LASSO because of many outliers 
result_cv <- result_cv[-which(methods == "scaledLASSO"), ]

# mean 8-fold cross validation over 100 replications
result_mean <- apply(result_cv, 1, mean)

# order methods based on result_mean
ord <- order(result_mean)
result_cv <- result_cv[ord, ]


# Plot CV error

labels <- c("var-SSL",
            "fixed-SSL",
            "MCP*",
            "SCAD",
            "LASSO",
            "MCP**",
            "ada-LASSO"
)

labels <- labels[ord]

dat <- melt(result_cv)

g <- ggplot(data = dat, aes(Var1, value)) 
g <- g + geom_boxplot(aes(fill = Var1))
g <- g + labs(x = " ", y = "CV error")
g <- g + scale_x_discrete(labels = labels)
g <- g + theme(axis.text.x=element_text(color = "black", size=15, angle=0)) 
g <- g  + theme(legend.position="none")
g <- g + theme(axis.ticks = element_blank())
g <- g + theme(text = element_text(size=20))
g








