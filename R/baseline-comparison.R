###############################################################################
# Comparison of Shapley values computed against different reference ("baseline")
# distributions:
#    A: Global, independent features
#    B: Global, correlated features
#    C: Local, independent features (drawn in a small neighborhood around the target)
#    D: Single baseline
###############################################################################

# Set seed and load necessary libraries
set.seed(1234)
library(xtable)

###############################################################################
# Function definitions

#' Convert 1d p-vectors to 2d 1xp matrices and leave matrices in original shape
#'
#' @param x vector or matrix to be converted to consistent shape
#' @return x converted to 2 dimensions (unchanged if already 2d)
convert.2d <- function(x){
  if (is.null(dim(x))) X <- t(as.matrix(x))
  else X <- x
  return(X)
}

#' Linear, additive function in three features
#'
#' @param x input data (vector or matrix form)
#' @return Column vector of function evaluations on X
fx.linear.additive <- function(x) {
  X <- convert.2d(x)
  return(-2*X[,1] + 1.5*X[,2] + 0.5*X[,3])
}

#' Linear, interactive function in three features
#'
#' @param x input data (vector or matrix form)
#' @return Column vector of function evaluations on X
fx.linear.interaction <- function(x) {
  X <- convert.2d(x)
  return(-2*X[,1] + 1.5*X[,2] + 0.5*X[,3] - 2*X[,2]*X[,3])
}

#' Non-linear, additive function in three features
#'
#' @param x input data (vector or matrix form)
#' @return Column vector of function evaluations on X
fx.nonlinear.additive <- function(x) {
  X <- convert.2d(x)
  return(-2*sin(X[,1]) + 1.5*abs(X[,2]) + 0.125*(X[,3]^2))
}

#' Non-linear, interactive function in three features
#'
#' @param x input data (vector or matrix form)
#' @return Column vector of function evaluations on X
fx.nonlinear.interaction <- function(x) {
  X <- convert.2d(x)
  return(-2*sin(X[,1]) + 1.5*abs(X[,2]) + 0.125*(X[,3]^2) + cos(X[,2]*X[,3]))
}

#' Create synthetic samples for correlated multivariate normal distribution
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param mean_vector The vector of means of the full MVN distribution
#' @param sigma_matrix Variance-covariance matrix of the full MVN distribution
#' @param cond_inds indices of the variables which will be conditioned on (fixed at their target value)
#' @param n Number of samples to draw
#' @return n by p matrix of synthetic samples
synthetic.samples.correlated <- function(x_target, mean_vector, sigma_matrix, 
                                         cond_inds, n){
  p <- ncol(convert.2d(x_target))
  if ((length(cond_inds) > 0) & length(cond_inds) < p){
    uncond_inds <- (1:p)[!(1:p %in% cond_inds)]
    Sigma.yy <- sigma_matrix[uncond_inds, uncond_inds]
    Sigma.yx <- sigma_matrix[uncond_inds, cond_inds]
    Sigma.xy <- sigma_matrix[cond_inds, uncond_inds]
    Sigma.xx <- sigma_matrix[cond_inds, cond_inds]
    mu.y <- mean_vector[uncond_inds]
    mu.x <- mean_vector[cond_inds]
    x <- x_target[cond_inds]
    cond_mean <- as.numeric(mu.y + (Sigma.yx %*% solve(Sigma.xx)) %*% (x - mu.x))
    cond_var <- Sigma.yy - Sigma.yx %*% solve(Sigma.xx) %*% Sigma.xy
    X_uncond <- rmvnorm(n, mean = cond_mean, sigma = cond_var)
    X_sample <- matrix(NA, nrow = n, ncol = p)
    X_sample[, uncond_inds] <- X_uncond
    for (i in 1:length(cond_inds)){
      X_sample[, cond_inds[i]] <- x_target[i]
    }
  } else if (length(cond_inds) == 0){
    X_sample <- rmvnorm(n, mean = mean_vector, sigma = sigma_matrix)
  } else if (length(cond_inds) == p){
    X_sample <- matrix(rep(x_target, n), nrow = n, byrow = T)
  }
  return(X_sample)
}

#' Create synthetic samples for independent normal distributions
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param indep_mean Mean of each independent normal distribution
#' @param indep_sd Standard deviation of each independent normal distribution
#' @param cond_inds indices of the variables which will be conditioned on (fixed at their target value)
#' @param n Number of samples to draw
#' @return n by p matrix of synthetic samples
synthetic.samples.independent <- function(x_target, indep_mean, indep_sd, 
                                          cond_inds, n){
  p <- ncol(convert.2d(x_target))
  uncond_inds <- (1:p)[!(1:p %in% cond_inds)]
  p_uncond <- length(uncond_inds)
  X_uncond <- matrix(rnorm(n*p_uncond, mean = indep_mean, sd = indep_sd), 
                     nrow = n, byrow = T)
  X_sample <- matrix(NA, nrow = n, ncol = p)
  X_sample[, uncond_inds] <- X_uncond
  for (i in 1:length(cond_inds)){
    X_sample[, cond_inds[i]] <- x_target[i]
  }
  return(X_sample)
}

#' Create synthetic samples for local normal distributions
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param indep_bandwidth Standard deviation of each independent normal distribution drawn around x_target
#' @param cond_inds indices of the variables which will be conditioned on (fixed at their target value)
#' @param n Number of samples to draw
#' @return n by p matrix of synthetic samples
synthetic.samples.local <- function(x_target, indep_bandwidth, 
                                    cond_inds, n){
  p <- ncol(convert.2d(x_target))
  X_sample <- rmvnorm(n, mean = x_target, sigma = (indep_bandwidth^2)*diag(p))
  for (i in 1:length(cond_inds)){
    X_sample[, cond_inds[i]] <- x_target[i]
  }
  return(X_sample)
}

#' Compute the conditional expectation across multiple baseline values drawn from some distribution
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param density_type The baseline distribution, must be one of "global independent", "global correlated", or "local independent"
#' @param indep_mean Mean of each independent normal distribution (only used if density_type = "global independent")
#' @param indep_sd Standard deviation of each independent normal distribution (only used if density_type = "global independent")
#' @param mean_vector The vector of means of the full MVN distribution (only used if density_type = "global correlated")
#' @param sigma_matrix Variance-covariance matrix of the full MVN distribution (only used if density_type = "global correlated")
#' @param indep_bandwidth Standard deviation of each independent normal distribution drawn around x_target (only used if density_type = "local independent")
#' @param cond_inds Indices of features to be conditioned / frozen
#' @param n Number of samples to draw
#' @return Expected value of f(X) integrating over all variables not in cond_inds
cond.expectation <- function(x_target, func, density_type, indep_mean, indep_sd, 
                             mean_vector, sigma_matrix, indep_bandwidth, 
                             cond_inds, n){
  if (density_type == "global independent"){
    X_sample <- synthetic.samples.independent(x_target, indep_mean, indep_sd, 
                                              cond_inds, n)
  } else if (density_type == "global correlated"){
    X_sample <- synthetic.samples.correlated(x_target, mean_vector, sigma_matrix, 
                                             cond_inds, n)
  } else if (density_type == "local independent"){
    X_sample <- synthetic.samples.local(x_target, indep_bandwidth, cond_inds, n)
  }
  mean(func(X_sample))
}

#' Compute the conditional expectation against a single baseline
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param x_baseline The value of x at the single "baseline" point used as a reference
#' @param cond_inds Indices of features to be conditioned / frozen
#' @return Function evaluation for the synthetic blend of x_target and x_baseline
cond.expectation.single.baseline <- function(x_target, func, x_baseline, cond_inds){
  x_sample <- x_baseline
  x_sample[cond_inds] <- x_target[cond_inds]
  return(func(x_sample))
}

#' Functional ANOVA effect for a given reference distribution
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param density_type The baseline distribution, must be one of "global independent", "global correlated", or "local independent"
#' @param indep_mean Mean of each independent normal distribution (only used if density_type = "global independent")
#' @param indep_sd Standard deviation of each independent normal distribution (only used if density_type = "global independent")
#' @param mean_vector The vector of means of the full MVN distribution (only used if density_type = "global correlated")
#' @param sigma_matrix Variance-covariance matrix of the full MVN distribution (only used if density_type = "global correlated")
#' @param indep_bandwidth Standard deviation of each independent normal distribution drawn around x_target (only used if density_type = "local independent")
#' @param u Indices of features to be conditioned / frozen
#' @param n Number of samples to draw
#' @return Contrast of conditional expectations f_u
f_u <- function(x_target, func, density_type, indep_mean, indep_sd, 
                mean_vector, sigma_matrix, indep_bandwidth, u, n){
  out <- 0
  n_l <- length(u)
  if (n_l > 0){
    for (i in 1:n_l){
      if (n_l == 1){
        inds <- list(u)
      } else{
        inds <- combn(u, i, simplify = F)
      }
      eff.coef <- (-1)^(n_l - i)
      out <- out + eff.coef*sum(sapply(
        inds, function(z) cond.expectation(
          x_target, func, density_type, indep_mean, indep_sd, 
          mean_vector, sigma_matrix, indep_bandwidth, z, n
        )))
    }
  }
  # Add / subtract the overall expectation
  eff.coef <- (-1)^(n_l)
  out <- out + eff.coef*cond.expectation(
    x_target, func, density_type, indep_mean, indep_sd, 
    mean_vector, sigma_matrix, indep_bandwidth, numeric(0), n)
  return(out)
}

#' Functional ANOVA effect for a single baseline
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param x_baseline The value of x at the single "baseline" point used as a reference
#' @param u Indices of features to be conditioned / frozen
#' @return Contrast of conditional expectations f_u
f_u_single_baseline <- function(
    x_target, func, x_baseline, u){
  out <- 0
  n_l <- length(u)
  if (n_l > 0){
    for (i in 1:n_l){
      if (n_l == 1){
        inds <- list(u)
      } else{
        inds <- combn(u, i, simplify = F)
      }
      eff.coef <- (-1)^(n_l - i)
      out <- out + eff.coef*sum(sapply(
        inds, function(z) cond.expectation.single.baseline(
          x_target, func, x_baseline, z
        )))
    }
  }
  # Add / subtract the overall expectation
  eff.coef <- (-1)^(n_l)
  out <- out + eff.coef*cond.expectation.single.baseline(
    x_target, func, x_baseline, numeric(0))
  return(out)
}

#' Enumerate the power set of feature indices for the number of features in x_target
#'
#' @param x_target The value of x at the "target" point to be explained
#' @return List of subsets of {1, ..., p} where p = length(x_target)
create.powerset.inds <- function(x_target){
  p <- ncol(convert.2d(x_target))
  shap_inds <- list(numeric(0))
  if (p > 0){
    for (i in 1:p){
      if (p == 1){
        inds <- list(1)
      } else{
        inds <- combn(1:p, i, simplify = F)
      }
      shap_inds <- append(shap_inds, inds)
    }
  }
  return(shap_inds)
}

#' Convert a list of features included in a set to a length p binary vector with 
#' 1 if an index is in x and 0 otherwise
#'
#' @param x Vector containing the elements included in a feature subset
#' @param p Total number of features under consideration
#' @return Binary vector of length p
inds.to.binary.vector <- function(x, p){
  out <- rep(0, p)
  out[x] <- 1
  return(out)
}

#' Convert a powerset list (output of create.powerset.inds) into a binary 
#' matrix of with as many rows as elements in the list shap_inds and p columns
#' indicating whether or not a feature is "active" in a given subset 
#'
#' @param shap_inds List of subsets of {1, ..., p} where p = length(x_target)
#' @param p Total number of features under consideration
#' @return Binary matrix with dimension length(shap_inds) x p
create.index.matrix <- function(shap_inds, p){
  t(sapply(shap_inds, function(x) inds.to.binary.vector(x, p)))
}

#' Shapley values for a given reference distribution
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param density_type The baseline distribution, must be one of "global independent", "global correlated", or "local independent"
#' @param indep_mean Mean of each independent normal distribution (only used if density_type = "global independent")
#' @param indep_sd Standard deviation of each independent normal distribution (only used if density_type = "global independent")
#' @param mean_vector The vector of means of the full MVN distribution (only used if density_type = "global correlated")
#' @param sigma_matrix Variance-covariance matrix of the full MVN distribution (only used if density_type = "global correlated")
#' @param indep_bandwidth Standard deviation of each independent normal distribution drawn around x_target (only used if density_type = "local independent")
#' @param n Number of samples to draw
#' @return Shapley values for all p features in x_target
shapley.values <- function(x_target, func, density_type, indep_mean, indep_sd, 
                           mean_vector, sigma_matrix, indep_bandwidth, n){
  p <- ncol(convert.2d(x_target))
  # Create a list of the powerset of Shapley subsets of 1:p
  shap_inds <- create.powerset.inds(x_target)
  # Create a matrix of binary indices of the Shapley vector
  shap_binary_matrix <- create.index.matrix(shap_inds, p)
  # Rescale to sum to 1
  shap_weight_matrix <- (shap_binary_matrix / 
                           ifelse(rowSums(shap_binary_matrix) > 0, rowSums(shap_binary_matrix), 1)
  )
  # Calculate the functional ANOVA effects
  fanova_effects <- sapply(shap_inds, function(x) f_u(
    x_target, func, density_type, indep_mean, indep_sd, 
    mean_vector, sigma_matrix, indep_bandwidth, x, n))
  # Calculate each of the Shapley values
  return(as.numeric(fanova_effects %*% shap_weight_matrix))
}

#' Shapley values for a single baseline
#'
#' @param x_target The value of x at the "target" point to be explained
#' @param func Function f(X) under evaluation
#' @param x_baseline The value of x at the single "baseline" point used as a reference
#' @return Shapley values for all p features in x_target
shapley.values.single.baseline <- function(
    x_target, func, x_baseline){
  p <- ncol(convert.2d(x_target))
  # Create a list of the powerset of Shapley subsets of 1:p
  shap_inds <- create.powerset.inds(x_target)
  # Create a matrix of binary indices of the Shapley vector
  shap_binary_matrix <- create.index.matrix(shap_inds, p)
  # Rescale to sum to 1
  shap_weight_matrix <- (
    shap_binary_matrix / ifelse(rowSums(shap_binary_matrix) > 0, rowSums(shap_binary_matrix), 1)
  )
  # Calculate the functional ANOVA effects
  fanova_effects <- sapply(shap_inds, function(x) f_u_single_baseline(
    x_target, func, x_baseline, x))
  # Calculate each of the Shapley values
  return(as.numeric(fanova_effects %*% shap_weight_matrix))
}

###############################################################################
# Example computations

n <- 1000000
x_target <- c(1,1,1)
x_baseline <- c(0,0,0)
p <- length(x_target)
sigma_mat <- matrix(c(1,0.9,0.5,0.9,1,0.75,0.5,0.75,1), byrow = T, nrow = p)

# 1. Linear additive shapley values 
global.indep.shap <- shapley.values(
  x_target, func = fx.linear.additive, density_type = "global independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
global.correl.shap <- shapley.values(
  x_target, func = fx.linear.additive, density_type = "global correlated", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
local.indep.shap <- shapley.values(
  x_target, func = fx.linear.additive, density_type = "local independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
single.baseline.shap <- shapley.values.single.baseline(
  x_target, fx.linear.additive, x_baseline)
xtable(rbind(global.indep.shap, global.correl.shap, local.indep.shap, single.baseline.shap))

# 2. Linear interaction shapley values
global.indep.shap <- shapley.values(
  x_target, func = fx.linear.interaction, density_type = "global independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
global.correl.shap <- shapley.values(
  x_target, func = fx.linear.interaction, density_type = "global correlated", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
local.indep.shap <- shapley.values(
  x_target, func = fx.linear.interaction, density_type = "local independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
single.baseline.shap <- shapley.values.single.baseline(
  x_target, fx.linear.interaction, x_baseline)
xtable(rbind(global.indep.shap, global.correl.shap, local.indep.shap, single.baseline.shap))

# 3. Nonlinear additive shapley values 
global.indep.shap <- shapley.values(
  x_target, func = fx.nonlinear.additive, density_type = "global independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
global.correl.shap <- shapley.values(
  x_target, func = fx.nonlinear.additive, density_type = "global correlated", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
local.indep.shap <- shapley.values(
  x_target, func = fx.nonlinear.additive, density_type = "local independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
single.baseline.shap <- shapley.values.single.baseline(
  x_target, fx.nonlinear.additive, x_baseline)
xtable(rbind(global.indep.shap, global.correl.shap, local.indep.shap, single.baseline.shap))

# 4. Nonlinear interaction shapley values
global.indep.shap <- shapley.values(
  x_target, func = fx.nonlinear.interaction, density_type = "global independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
global.correl.shap <- shapley.values(
  x_target, func = fx.nonlinear.interaction, density_type = "global correlated", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
local.indep.shap <- shapley.values(
  x_target, func = fx.nonlinear.interaction, density_type = "local independent", 
  indep_mean = 0, indep_sd = 1, mean_vector = rep(0, p), 
  sigma_matrix = sigma_mat, indep_bandwidth = 0.25, n = n)
single.baseline.shap <- shapley.values.single.baseline(
  x_target, fx.nonlinear.interaction, x_baseline)
xtable(rbind(global.indep.shap, global.correl.shap, local.indep.shap, single.baseline.shap))
