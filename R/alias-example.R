###############################################################################
# This script implements the example given in Section 3.1.1, which analyzes
# the first-order interaction aliases of a 12-sample SHAP design
###############################################################################

# Load necessary packages
require(FrF2)

# Generate the full 2^p design matrix
p <- 6
A <- fac.design(2, p, randomize = F)
B <- ifelse(as.matrix(A) == "2", 1, 0)
E <- rowSums(B)
X.orig <- unname(B[order(E),])

# Remove the first and last entries, mirroring the constraints imposed by the SHAP library
X.red <- X.orig[-c(1, nrow(X.orig)),]

# Compute all higher-order interactions
X <- model.matrix(~0+(.)^5, data.frame(X.red))

# Compute SHAP regression weights
shap_weight <- function(X){
  p <- ncol(X)
  z <- apply(X, 1, sum)
  comb <- factorial(p)/(factorial(p - z)*factorial(z))
  return((p - 1)/(comb*z*(p - z)))
}
W <- diag(shap_weight(X.red))

# Remove the p-th column of the full interaction matrix and subtract its value
X.constrained <- X[,-p] - X[,p]
X.r.star <- X.constrained[,1:(p-1)]
X.ji.star <- X.constrained[,(p):ncol(X.constrained)]

# Aliasing with the full factorial
(alias.1 = solve(t(X.r.star) %*% W %*% X.r.star) %*% t(X.r.star) %*% W %*% X.ji.star)

# Aliasing with the partial factorial
# Paired sampling
paired_sample_inds <- ((rowSums(X.red) == (p-1)) | (rowSums(X.red) == (1)))
X.constrained.paired <- X.constrained[paired_sample_inds,]
W.paired <- W[paired_sample_inds,paired_sample_inds]
X.r.star <- X.constrained.paired[,1:(p-1)]
X.ji.star <- X.constrained.paired[,(p):ncol(X.constrained)]
(alias.2 = solve(t(X.r.star) %*% W.paired %*% X.r.star) %*% t(X.r.star) %*% W.paired %*% X.ji.star)

# Observe that for first and second order interactions,
# alias 1 and 2 are exactly the same
sum(abs(alias.1[,1:(2*p-3)] - alias.2[,1:(2*p-3)]))
