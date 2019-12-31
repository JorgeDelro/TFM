#

library(tidyverse)          # data manipulation
library(vars)               # VAR function
library(rstan)
library(fpp2)               # for forecast and plot
library(BigVAR)             # for simulation - Simulate a VAR process
library(mcompanion)         # for simulation - Generate a Multi-companion matrix
library(clusterGeneration)  # for simulation - Generate a positive definite matrix/covariance matrix

# Establish random number generator seed
# for reproducible results
set.seed(123)

# Simulation design
# 5 - dimensional VAR(p) where p could be: 
# p = 5 i.e. a microcycle
# p = 20 i.e. a mesocycle
# p = 260 i.e. a season

# Number of series
k <- 5
# Lag order
p <- 5
# Number of time points
nT = 100
# Example eigenvalues
# eigenvalues <- c(0.8266574, 0.6765737, 0.2216213, 0.6765737, 0.2216213)


res <- genstaVAR(k = 5,
                 p = 5,
                 nT = 100,
                 Sigma = NULL,
                 eigenvalues = NULL)

# Check Eigenvalues
Mod(eigen(res$Phi)$values)
# Plot time-series
autoplot(ts(res$Z), facets = T)

# Data to fit the model
Z_train <- ts(res$Z[1:80,])
# Data to assess forecasting
Z_forecast <- ts(res$Z[81:100,], start = 81, end = 100)

## package vars
mod_VAR <- VAR(y = Z_train,
               p = p,
               type = "const")

summary(mod_VAR)
mod_VAR_forecast <- predict(mod_VAR, n.ahead = 20, ci = 0.95, dumvar = NULL)
# print(mod_VAR_forecast)
# plot(mod_VAR_forecast)
Serie1 <- mod_VAR_forecast$fcst$Series.1[,1]
#colnames(Serie1) <- c("Point Forecast", "Lo 95", "Hi 95")
#class(Serie1) <- "forecast"

accuracy(Serie1, Z_forecast[,1])



##
stan_data <- list(nT = dim(Y_train)[1], # Length
                  K = k, # Variables
                  P = p, # Lag
                  Z = Y_train,
                  nF = 20)


model_w_data <- stan(file = 'VAR_p.stan', data=stan_data, iter=200, chains=4)
summary(model_w_data)
# Names parameters
print(names(model_w_data))

ols_forecast <- forecast(mod_VAR, h = 20, level = 95)
ols_forecast$forecast$Series.1
accuracy(ols_forecast$forecast$Series.1, Z_forecast[,1])
Serie1 <- ols_forecast$forecast$Series.1
xSerie1 <- Z_forecast[,1]

accuracy(Serie1, xSerie1)

forecast(mod_VAR) %>%
  autoplot() + xlab("Year")





