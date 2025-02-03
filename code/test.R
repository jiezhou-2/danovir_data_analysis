n<- 100; p <- 5; true_p <- 2
set.seed(950)
X <- matrix(rnorm(n * p), nrow = n)
beta <- matrix(c(rep(1, true_p), rep(0, p - true_p)), ncol = 1)
y <- X %*% beta + 3 * rnorm(n)


X_centered <- apply(X, 2, function(x) x - mean(x))
Xs <- apply(X_centered, 2, function(x) x / sqrt(sum(x^2) / n))

library(glmnet)
fit <- glmnet(Xs, y, standardize = TRUE)
