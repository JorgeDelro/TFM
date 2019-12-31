
library(rstan)
library(mvtnorm)
set.seed(123)

L <- 250 # length
M <- 2    # variables
p <- 2   # lags

alpha <- c(-0.1,0.2)
A.2 <- matrix(c(0.85,0.05,-0.1,0.35), nrow = M)
A.1 <- eigen(A.2)$vectors %*% diag(eigen(A.2)$values / -2) %*% solve(eigen(A.2)$vectors)
Sigma <- matrix(c(1, -0.2, -0.2, 0.25),2,2)

y <- array(NA, c(L, M))
y[1, ] <- rnorm(M, 0, 0.3)
y[2, ] <- rnorm(M, 0, 0.3)
y[3, ] <- rnorm(M, 0, 0.3)

for(i in 3:T) {
   y[i, ] <- alpha + A.2 %*% y[(i-1),]  + A.1 %*% y[(i-2),] + c(rmvnorm(1, mean = rep(0,M) , sigma = Sigma))
}
stan_data <- list(L = L,
                  M = M,
                  P = p,
                  Y = y)

fitted.model <- stan(file = 'VAR_p.stan')
model_w_data <- stan(file = 'VAR_p.stan', data=stan_data, iter=200, chains=4)

c(t(cbind(alpha,A.1)))
print(model_w_data)

library(vars)

summary(VAR(y,p = 2))
