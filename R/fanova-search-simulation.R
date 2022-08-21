###############################################################################
# Implementation of the functional ANOVA interaction search algorithm of 
# Hooker (2004), "Discovering Additive Structure in Black Box Functions"
###############################################################################

# Load necessary libraries
pkg.list <- c("sensitivity", "lhs")
lapply(pkg.list, require, character.only = TRUE)

###############################################################################
# Function defintions

#' Convert 1d p-vectors to 2d 1xp matrices and leave matrices in original shape
#'
#' @param x Vector or matrix to be converted to consistent shape
#' @return x converted to 2 dimensions (unchanged if already 2d)
convert.2d <- function(x){
  if (is.null(dim(x))) X <- t(as.matrix(x))
  else X <- x
  return(X)
}

#' f(X) being studied and analyzed using the Hooker (2004) algorithm
#'
#' @param x Input covariates (vector or matrix form)
#' @return Column vector of F(X) values evaluated on `x`
y.func <- function(x) {
  X <- convert.2d(x)
  return(X[,1] + X[,2] + X[,3] + X[,2]*X[,3])
}

#' "Pick-freeze" approximation of the L2 Cost of Exclusion (L2COE) introduced in
#' Liu and Owen (2006) "Estimating Mean Dimensionality of Analysis of Variance
#' Decompositions"
#'
#' @param func Function for which sensitivities are to be assessed
#' @param u Subset of p features under evaluation
#' @param n_s Number of samples to be drawn at random
#' @param p Number of features in the dataset (and expected by func)
#' @return Estimate of the cost of exclusion of u features
L2COE.liu.owen <- function(func, u, n_s, p){
  # Draw sample data
  X <- matrix(runif(n_s*p), ncol = p)
  Z <- matrix(runif(n_s*p), ncol = p)
  
  # Construct each of the 2^(length(u)) expectations
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
      for (j in 1:length(inds)){
        # Construct the resampled dataset
        synth_data <- X
        synth_data[,inds[[j]]] <- Z[,inds[[j]]]
        out <- out + eff.coef*func(synth_data)
      }
    }
  }
  # Add / subtract the overall expectation
  eff.coef <- (-1)^(n_l)
  out <- out + eff.coef*func(X)
  
  # Sum and square for superset importance
  liu.owen.importance <- sum((out)^2)/((2^(n_l))*n_s)
  return(liu.owen.importance)
}

#' Estimate of Var(f_u) where f_u is the functional ANOVA for set u
#'
#' @param X_background Reference dataset against which to evaluate func
#' @param func Function for which sensitivities are to be assessed
#' @param u Size of of the subset of p features under evaluation
#' @return Estimate of the cost of exclusion of u features
sigma.f <- function(X_background, func, u){
  f_u_vals <- (
    apply(X_background, 1, function(z) f_u(
      z, X_background, func, u
    )))
  e_f_u <- mean(f_u_vals)
  return(mean((f_u_vals-e_f_u)^2))
}

#' Functional ANOVA effect for f(x_active) and set u
#'
#' @param x_active Specific covariate values that determine f(x)
#' @param X_background Reference dataset with which to compute expectations
#' @param func Function for which sensitivities are to be assessed
#' @param u Subset of p features under evaluation
#' @return Functional ANOVA term f_u
f_u <- function(x_active, X_background, func, u){
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
        inds, function(z) cond.expectation(x_active, X_background, func, z)
      ))
    }
  }
  # Add / subtract the overall expectation
  eff.coef <- (-1)^(n_l)
  out <- out + eff.coef*cond.expectation(x_active, X_background, func, numeric(0))
  return(out)
}

#' Compute the conditional expectation of a function with some variables in the conditioning set
#'
#' @param x_active The current value of x, at which conditioned features should be fixed
#' @param X_background A set of background samples from f(X) which are used to compute the expectation for unconditioned features
#' @param func Function for which sensitivities are to be assessed
#' @param cond_inds Indices of features to be fixed
#' @return Conditional expectation of f(X) | X[,cond_inds]
cond.expectation <- function(x_active, X_background, func, cond_inds){
  X_synth <- X_background
  X_synth[,cond_inds] <- x_active[cond_inds]
  mean(func(X_synth))
}

###############################################################################
##### Simulation study: recover feature interaction ranking

# Draw data 
# X1, X2, X3 ~ Unif(0,1)
# y = F(x) = X1 + X2 + X3 + X2 X3
n <- 1000
p <- 3
X <- matrix(runif(n*p), ncol = p)
y <- y.func(X)

# Estimate overall variance of F to scale the sigma.f terms in the algorithm below
V_f <- (sigma.f(X, y.func, numeric(0)) + sigma.f(X, y.func, c(1)) + 
          sigma.f(X, y.func, c(2)) + sigma.f(X, y.func, c(3)) + 
          sigma.f(X, y.func, c(1,2)) + sigma.f(X, y.func, c(1,3)) + 
          sigma.f(X, y.func, c(2,3)) + sigma.f(X, y.func, c(1,2,3))
)

# Run the breadth-first algorithm to select active elements of 
# the 2^p powerset of functional ANOVA terms
eps <- 0.01
n_sim <- 100
correct_ranking <- matrix(NA, nrow = n_sim)
for (sim in 1:n_sim){
  # Simulation parameters
  S <- 0
  U <- list()
  iter <- 1

  # Data 
  n <- 500
  n_s <- 500
  p <- 3
  X <- matrix(runif(n*p), ncol = p)
  y <- y.func(X)
  
  # Estimate overall variance of F to scale the sigma.f terms in the algorithm below
  V_f <- (sigma.f(X, y.func, numeric(0)) + sigma.f(X, y.func, c(1)) + 
            sigma.f(X, y.func, c(2)) + sigma.f(X, y.func, c(3)) + 
            sigma.f(X, y.func, c(1,2)) + sigma.f(X, y.func, c(1,3)) + 
            sigma.f(X, y.func, c(2,3)) + sigma.f(X, y.func, c(1,2,3))
  )
  # Initialize the set of candidates
  for (i in 1:p){
    K1 <- combn((1:p), i, simplify = F)
    coe.vec <- sapply(K1, function(term) L2COE.liu.owen(
      y.func, term, n_s, p
    ))
    sort.inds <- order(coe.vec, decreasing = T)
    for (j in sort.inds){
      U[[iter]] <- K1[[j]]
      S <- S + sigma.f(X, y.func, K1[[j]])
      iter <- iter + 1
      # Terminate loop if S crosses the variance threshold
      if (S > V_f*(1 - eps)){
        break
      }
    }
    # Terminate loop if S crosses the variance threshold
    if (S > V_f*(1 - eps)){
      break
    }
  }
  A <- as.list(sort(c(U[[1]], U[[2]])))
  sorted_ranking <- list(
    A[[1]], A[[2]], U[[3]], U[[4]]
  )
  correct_ranking[sim] <- (prod(c(
    (sorted_ranking[[1]] == 2)*1, 
    (sorted_ranking[[2]] == 3)*1, 
    (sorted_ranking[[3]] == 1)*1, 
    (mean(sorted_ranking[[4]] == c(2,3)))*1
  )) == 1)
}

# Observe how often the correct ranking is recovered across the simulations
mean(correct_ranking)
