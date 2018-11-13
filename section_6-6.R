#-----------------------------------------------------------------------
# Code for the simulation study in Section 6.6 of the paper
# "Variance prior forms for high-dimensional Bayesian variable selection"
# https://arxiv.org/abs/1801.03019
#-----------------------------------------------------------------------

library(SSLASSO)
library(scaledSSL)
library(mvtnorm)
library(ncvreg)
library(scalreg)
library(parcor)
library(tidyverse)
library(reshape2)
library(xtable)

set.seed(12345)

# --------------- Keep track of the following quantities --------------- #
# PE ...prediction error
# FN .... false negatives
# FP .... false positives
# TP .... true positives
# TN .... true negatives
# DIM.... dimensionality
# COR....  found the correct model
# TIME... execution time
# HAM ... Hamming distance


# --------------- Methods compared --------------- #
# SSL 
# scaledSSL
# MCP
# SCAD
# LASSO
# Scaled Lasso
# MCP hard (gamma = 1.0001)

#---------------Parameters --------------- #
# p: number of predictors
# n: number of observations
# Y: response, n x 1 vector
# X: predictors, n x p matrix
# beta: regression coefficients, p x 1 vector
# sigma_true: true error standard deviation
# q: number of non-zero beta
# Sigma: covariance of predictors, X
# rho: correlation in blocks of Sigma
# nrep: number of replications

p <- 1000
n <- 100
q <- 6  
rho <- 0.9 
block <- matrix(rho, 50, 50) 
Sigma <- diag(20)%x%block  
diag(Sigma) <- 1
sigma_true <- sqrt(3)

nrep <- 100 

lambda1 <- 1
lambda0 <- seq(1, 100, by = 1)

#  Function to calculate predictive error
pred_error = function(tXX, beta, beta_hat) {
  diff = beta - beta_hat
  err = t(diff) %*% (tXX %*% diff)
  err
}

# ------------ Set up variables to store results --------#

methods <- c("varSSL",
             "oracleSSL",
             "SSL", 
             "scaledSSL",
             "MCP",
             "SCAD",
             "LASSO",
             "scaledLASSO",
             "hardMCP",
             "adaLASSO")

metrics <- c("PE", "FN", "FP", "TP", "TN", "DIM", "COR", "TIME", "HAM")

n_methods <- length(methods)
n_metrics <- length(metrics)

result <- array(0, c(n_methods, n_metrics, nrep))
colnames(result) <- metrics
rownames(result) <-  methods

sigmas <- matrix(0, nrow = nrep, ncol = length(methods) + 1)
colnames(sigmas)  <- c("Oracle", methods)

sigmas_adj <- matrix(0, nrow = nrep, ncol = length(methods) + 1)
colnames(sigmas_adj) <- c("Oracle", methods)

#---------------- Simulation Study --------------------------#

for(i in 1:nrep) {
  
  # Generate data
  
  X <- rmvnorm(n, numeric(p), Sigma)
  index <- c(1, 51, 101, 151, 201, 251)
  neg_index <- (1:p)[-c(index)]
  beta <- numeric(p)
  beta[index] <- c(-2.5,-2,-1.5,1.5,2,2.5) 
  model <- numeric(p)
  model[index] <- 1
  y <- X %*% beta + rnorm(n, sd = sigma_true)
  tXX <- crossprod(X)
  
  sigmas[i, "Oracle"]  <- sum((y - X%*%beta)^2)/n
  
  # Methods
  
  # Spike-and-Slab LASSO with unknown variance
  
  time <- system.time(fit <- SSLASSO(X, y, variance = "unknown"))
  lambda0_start <- which(fit$iter < 100)[1]
  lambda0 <- seq(lambda0_start, 100, by = 1)
  time <- system.time(fit <- SSLASSO(X, y, variance = "unknown", lambda1 = 1, lambda0 = lambda0))
  
  index_est <- fit$model
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  beta_est <- fit$beta[, ncol(fit$beta)]
  model_est <- numeric(p)
  model_est[index_est] <- 1
  sigma_est  <- (fit$sigmas[length(fit$sigmas)])^2
  
  sigmas[i ,"varSSL"] <- sigma_est 
  sigmas_adj[i, "varSSL"] <- sigma_est * (n + 2) / (n - length(index_est))
  
  result["varSSL", "PE", i] <- pred_error(tXX, beta, beta_est) 
  result["varSSL", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["varSSL", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["varSSL", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["varSSL", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["varSSL", "DIM", i] <- length(index_est)              
  result["varSSL", "COR", i] <- as.numeric(all(model == model_est))    
  result["varSSL", "TIME", i] <- as.numeric(time[3])          
  result["varSSL", "HAM", i] <- sum(abs(model_est - model)) 
  
  # Spike-and-Slab LASSO with fixed variance = true variance
  
  time <- system.time(fit <- SSLASSO(X, y, variance = "fixed", sigma = sigma_true, lambda1 = 1, lambda0 = lambda0))
  
  index_est <- fit$model
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  beta_est <- fit$beta[, ncol(fit$beta)]
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"oracleSSL"] <- sigma_est 
  sigmas_adj[i, "oracleSSL"] <- sigma_est * n / (n - length(index_est))
  
  result["oracleSSL", "PE", i] <- pred_error(tXX, beta, beta_est) 
  result["oracleSSL", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["oracleSSL", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["oracleSSL", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["oracleSSL", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["oracleSSL", "DIM", i] <- length(index_est)              
  result["oracleSSL", "COR", i] <- as.numeric(all(model == model_est))    
  result["oracleSSL", "TIME", i] <- as.numeric(time[3])          
  result["oracleSSL", "HAM", i] <- sum(abs(model_est - model)) 
  
  # Spike-and-slab lasso with fixed variance = 1
  
  time <- system.time(fit <- SSLASSO(X, y, variance = "fixed", sigma = 1, lambda1 = 1, lambda0 = lambda0))
  
  index_est <- fit$model
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  beta_est <- fit$beta[, ncol(fit$beta)]
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"SSL"] <- sigma_est 
  sigmas_adj[i, "SSL"] <- sigma_est * n / (n - length(index_est))
  
  result["SSL", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["SSL", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["SSL", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["SSL", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["SSL", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["SSL", "DIM", i] <- length(index_est)              
  result["SSL", "COR", i] <- as.numeric(all(model == model_est))    
  result["SSL", "TIME", i] <- as.numeric(time[3])          
  result["SSL", "HAM", i] <- sum(abs(model_est - model)) 
  
  # Scaled SSL
  
  time <- system.time(fit <- scaledSSL(X, y, variance = "unknown", lambda1 = 1, lambda0 = lambda0))
  
  index_est <- fit$model
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  beta_est <- fit$beta[, ncol(fit$beta)]
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"scaledSSL"] <- sigma_est 
  sigmas_adj[i, "scaledSSL"] <- sigma_est * n / (n - length(index_est))
  
  result["scaledSSL", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["scaledSSL", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["scaledSSL", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["scaledSSL", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["scaledSSL", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["scaledSSL", "DIM", i] <- length(index_est)              
  result["scaledSSL", "COR", i] <- as.numeric(all(model == model_est))    
  result["scaledSSL", "TIME", i] <- as.numeric(time[3])          
  result["scaledSSL", "HAM", i] <- sum(abs(model_est - model))
  
  # MCP
  
  time <- system.time(cvfit <- cv.ncvreg(X, y, family = c("gaussian"), penalty = "MCP", warn = FALSE))
  
  fit <- cvfit$fit
  beta_est <- fit$beta[-1, cvfit$min]
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"MCP"] <- sigma_est 
  sigmas_adj[i, "MCP"] <- sigma_est * n / (n - length(index_est))
  
  result["MCP", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["MCP", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["MCP", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["MCP", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["MCP", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["MCP", "DIM", i] <- length(index_est)              
  result["MCP", "COR", i] <- as.numeric(all(model == model_est))    
  result["MCP", "TIME", i] <- as.numeric(time[3])          
  result["MCP", "HAM", i] <- sum(abs(model_est - model))
  
  # SCAD
  
  time <- system.time(cvfit <- cv.ncvreg(X, y, family = c("gaussian"), penalty = "SCAD", warn = FALSE))
  
  fit <- cvfit$fit
  beta_est <- fit$beta[-1, cvfit$min]
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"SCAD"] <- sigma_est 
  sigmas_adj[i, "SCAD"] <- sigma_est * n / (n - length(index_est))
  
  result["SCAD", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["SCAD", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["SCAD", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["SCAD", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["SCAD", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["SCAD", "DIM", i] <- length(index_est)              
  result["SCAD", "COR", i] <- as.numeric(all(model == model_est))    
  result["SCAD", "TIME", i] <- as.numeric(time[3])          
  result["SCAD", "HAM", i] <- sum(abs(model_est - model)) 
  
  # LASSO
  
  time <- system.time(cvfit <- cv.ncvreg(X, y, family=c("gaussian"),penalty="lasso",warn=FALSE))
  
  fit <-cvfit$fit
  beta_est <- fit$beta[-1, cvfit$min]
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"LASSO"] <- sigma_est 
  sigmas_adj[i, "LASSO"] <- sigma_est * n / (n - length(index_est))
  
  result["LASSO", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["LASSO", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["LASSO", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["LASSO", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["LASSO", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["LASSO", "DIM", i] <- length(index_est)              
  result["LASSO", "COR", i] <- as.numeric(all(model == model_est))    
  result["LASSO", "TIME", i] <- as.numeric(time[3])          
  result["LASSO", "HAM", i] <- sum(abs(model_est - model))  
  
  # Scaled  LASSO
  
  time <- system.time(fit <- scalreg(X, y))
  
  beta_est <- fit$coefficients
  
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"scaledLASSO"] <- sigma_est 
  sigmas_adj[i, "scaledLASSO"] <- sigma_est * n / (n - length(index_est))
  
  result["scaledLASSO", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["scaledLASSO", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["scaledLASSO", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["scaledLASSO", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["scaledLASSO", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["scaledLASSO", "DIM", i] <- length(index_est)              
  result["scaledLASSO", "COR", i] <- as.numeric(all(model == model_est))    
  result["scaledLASSO", "TIME", i] <- as.numeric(time[3])          
  result["scaledLASSO", "HAM", i] <- sum(abs(model_est - model))
  
  
  
  # MCP  (gamma = 1.0001)
  
  time <- system.time(cvfit <- cv.ncvreg(X, y, family=c("gaussian"),penalty="MCP",warn=FALSE,gamma=1.0001))
  
  fit <- cvfit$fit
  
  beta_est <- fit$beta[-1, cvfit$min]
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"hardMCP"] <- sigma_est 
  sigmas_adj[i, "hardMCP"] <- sigma_est * n / (n - length(index_est))
  
  result["hardMCP", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["hardMCP", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["hardMCP", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["hardMCP", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["hardMCP", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["hardMCP", "DIM", i] <- length(index_est)              
  result["hardMCP", "COR", i] <- as.numeric(all(model == model_est))    
  result["hardMCP", "TIME", i] <- as.numeric(time[3])          
  result["hardMCP", "HAM", i] <- sum(abs(model_est - model)) 
  
  # Adaptive lasso
  
  time <- system.time(fit <- adalasso(X, y, k = 10, intercept=F))
  
  beta_est <- fit$coefficients.adalasso
  index_est <- which(beta_est != 0)
  
  if(length(index_est) == 0) {
    neg_index_est <- 1:p
  } else {
    neg_index_est <- (1:p)[-c(index_est)]
  }
  
  model_est <- numeric(p)
  model_est[index_est] <- 1
  
  sigma_est <- mean((y - X %*% beta_est)^2)
  
  sigmas[i ,"adaLASSO"] <- sigma_est 
  sigmas_adj[i, "adaLASSO"] <- sigma_est * n / (n - length(index_est))
  
  result["adaLASSO", "PE", i] <- pred_error(tXX, beta, beta_est)  
  result["adaLASSO", "FN", i] <- sum(index %in% index_est == FALSE) 
  result["adaLASSO", "FP", i] <- sum(index_est %in% index == FALSE)  
  result["adaLASSO", "TP", i] <- sum(index_est %in% index == TRUE)   
  result["adaLASSO", "TN", i] <- sum(neg_index_est %in% neg_index == TRUE)
  result["adaLASSO", "DIM", i] <- length(index_est)              
  result["adaLASSO", "COR", i] <- as.numeric(all(model == model_est))    
  result["adaLASSO", "TIME", i] <- as.numeric(time[3])          
  result["adaLASSO", "HAM", i] <- sum(abs(model_est - model))
  
  
  print(i)
  
}

stop()


# Set up data.frame for ggplot
sig_adj_melt <- melt(sigmas_adj[, c("varSSL", "scaledSSL", "scaledLASSO")])
colnames(sig_adj_melt) <- c("nrep", "method", "value")

# Plot sigmas (adjusted)
g <- ggplot(sig_adj_melt, aes(method, value))
g <- g + geom_boxplot(aes(fill = method))
g <- g + labs(x = " ", y = expression(hat(sigma)^2))
g <- g + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
g <- g + scale_fill_discrete(labels = c(expression(paste(" SSL (unknown ", sigma^2, ")")),
                                        "Scaled SSL",
                                        "Scaled LASSO"))
g <- g + theme(legend.title=element_blank()) + labs(x = NULL)
g <- g + theme(text = element_text(size=20))
g <- g + geom_hline(yintercept = sigma_true^2, color = "red")
g


# Generate table

# function to calculate Matthews Correlation Coefficient
calculate_MCC <- function(TP, FP, TN, FN) {
  MCC <- (TP * TN - FP * FN)/(sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))
  MCC[is.nan(MCC)] <- 0
  MCC
}

# convert array into tibble
result_tibble <- melt(result)
colnames(result_tibble) <- c("Method", "Metric", "Replication", "Value")
result_tibble <- as.tibble(result_tibble)
result_tibble <- spread(result_tibble, Metric, Value)

result_tibble <- result_tibble %>% 
  mutate(MCC = calculate_MCC(TP, FP, TN, FN)) %>%
  gather(key = "Metric", value = "Value", -c(Method, Replication)) %>%
  group_by(Method, Metric)

result_mean <- result_tibble %>%
  summarize(Mean = mean(Value)) %>%
  spread(Metric, Mean)

result_summary <- result_tibble %>% 
  summarize(Mean = paste(format(round(mean(Value), 1), nsmall = 1), " {\\scriptsize (", format(round(sd(Value)/10, 1), nsmall = 1), ")}", sep = "")) %>%
  spread(Metric, Mean) %>%
  select(-MCC, -TN, -COR, -DIM, -TIME)

result_MCC <- result_tibble %>%
  filter(Metric == "MCC") %>%
  summarize(Mean = paste(format(round(mean(Value), 2), nsmall = 2), " {\\scriptsize (", format(round(sd(Value)/10, 2), nsmall =2), ")}", sep = "")) 

result_COR <- result_tibble %>% 
  filter(Metric == "COR") %>%
  summarize(Sum = sum(Value))

result_TIME <- result_tibble %>%
  filter(Metric == "TIME") %>%
  summarize(Mean = paste(format(round(mean(Value), 2), nsmall = 2), " {\\scriptsize (", format(round(sd(Value)/10, 2), nsmall = 2), ")}", sep = ""))

result_summary <- as.data.frame(result_summary)
result_summary$COR <- result_COR$Sum
result_summary$MCC <- result_MCC$Mean
result_summary$TIME <- result_TIME$Mean

ord <- order(result_mean$HAM)
result_summary <- result_summary[ord, ]
result_summary <- result_summary[, c("Method", "HAM", "PE", "MCC", "TP", "FP", "FN", "COR", "TIME")]

methods_latex <- c("SSL (unknown $\\sigma^2$)",
                   "SSL (fixed $\\sigma^2 = 3$)", 
                   "SSL (fixed $\\sigma^2 = 1$)",
                   "Scaled SSL",
                   "MCP*",
                   "SCAD",
                   "LASSO",
                   "Scaled LASSO",
                   "MCP**",
                   "Adaptive LASSO")

methods_latex <- methods_latex[ord]

result_summary$Method <- methods_latex

print(xtable(result_summary, digits = rep(0, ncol(result_summary)+1)), sanitize.text.function = function(x) {x}, 
      include.rownames = F)

