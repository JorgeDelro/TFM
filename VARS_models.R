#

library(tidyverse)
library(mvtnorm)
library(MTS)

# simulate using loop
# https://faculty.washington.edu/ezivot/book/Ch11.VARexamples2ndEdition.ssc
nobs = 250
set.seed(301)
# e.var = rmvnorm(nobs,sd=c(1,1),rho=0.5)
pi1 = matrix(c(0.7,0.2,0.2,0.7),2,2)
mu.vec = c(1,5)
sigma <- matrix(c(1,0.5,0.5,1), ncol=2)
e.var <- rmvnorm(n=nobs, mean = rep(0, nrow(sigma)), sigma=sigma)
e.var = t(e.var)
c.vec = as.vector((diag(2)-pi1)%*%mu.vec)
y.var = matrix(0,2,nobs)
y.var[,1] = mu.vec
for (i in 2:nobs) {
  y.var[,i] = c.vec+pi1%*%y.var[,i-1]+e.var[,i]
}
y.var = t(y.var)
dimnames(y.var) = list(NULL,c("y1","y2"))
eigen(pi1,only.values=T)
c.vec
colMeans(y.var)
tsplot(y.var,lty=1:2,ylab="y1, y2")
legend(200,9,legend=colIds(y.var),lty=1:2)

# package MTS
library(MTS)
mod_MTS <- MTS::VAR(x = y.var.ts,
                    output = T,
                    include.mean = T,
                    fixed = NULL)

# package VARS
library(vars)
mod_vars <- vars::VAR(y = y.var.ts,
                      p = 1,
                      type = "none")

mod_vars_forecast <- vars::VARselect(y = y.var.ts,
                                     lag.max = 10,
                                     type = "const")

forecast(mod_vars) %>% autoplot() 

### Simulation VAR
install.packages("BigVAR")
library(BigVAR)
# k = Number of series
k=3
# p = Maximum lag order
p=6
# B = Time series
B=matrix(0,nrow=k,ncol=p*k)

A1<- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
B[,1:k]=A1
B[,(4*k+1):(5*k)]=A2
# VarptoVar1MC returns a kp \times kp coefficient matrix representing all coefficient matrices contained in Ai as a VAR(1).
A <- VarptoVar1MC(B,p,k)
# MultVarSim returns a T \times k of realizations from a VAR.
# k, number of series
# A, coefficient matrix
# p, maximum lag order
# Sigma, residual covariance matrix
# T, number of simulations
Y <-MultVarSim(k,A,p,.1*diag(k),100)
Y <- ts(Y)
autoplot(Y, facets = T)

########## Forecasting
# Yang 2018
# the evaluation metrics used in this paper are: sensitivity, specificity, PPV, NPV, F1 score, sMAPE, RMSE and MAD.

# Pendar 2017
# First sec- tion is 1352-1392 that is used for estimation models 
# and second section is 1393-1394 that is used for inflation forecast during these two years.
# To evaluate the forecasts in this article forecast error criteria RMSE and MAPE have been used.

# Adenomom 2015
# short term (T = 8, 16)
# medium term (T = 32, 64)
# long term (T = 128, 256).

# Nilsson 2016
# The VAR(1) model includes 50 simultaneous equations, one of which will be used as a reference
# The models will attempt to predict one, four, eight, twelve, twenty and thirty steps ahead in time
# The analysis will therefore be conducted in three different scenarios, 
# a scenario with a training set of 320 observations, 
# a scenario with 160 observations, 
# a scenario with 80 observations, 
# and an additional dataset of 30 observations for each scenario are set aside to be evaluated with RMSE. 
# The analysis will be repeated 5000 times and the measurements of error will be averaged to produce reliable results.

# Data to fir the model
Y_train <- Y[1:80,]
# Data to assess forecasting
Y_forecast <- Y[81:100,]

# Vector autoregression
library(vars)
#Estimates a VAR by OLS per equation. The model is of the following form:
#  \bold{y}_t = A_1 \bold{y}_{t-1} + … + A_p \bold{y}_{t-p} + CD_t + \bold{u}_t
# where \bold{y}_t is a K \times 1 vector of endogenous variables and u_t assigns a spherical disturbance term of the same dimension. 
# The coefficient matrices A_1, …, A_p are of dimension K \times K. In addition, either a constant and/or a trend can be included 
# as deterministic regressors as well as centered seasonal dummy variables and/or exogenous variables 
# (term CD_T, by setting the type argument to the corresponding value and/or setting season to the desired frequency (integer) 
# and/or providing a matrix object for exogen, respectively. The default for type is const and for season and exogen the default is set to NULL.
# If for lag.max an integer value is provided instead of NULL (the default), the lag length is determined by the selected information criteria in ic, the default is Akaike.
mod_VAR <- VAR(y = Y_train,
               p = 6,
               type = "const",
               season = NULL, 
               exogen = NULL, 
               lag.max = NULL)
# Diagnostic testing
# Heteroscedasticity: multivariate ARCH test
arch.test(mod_VAR, lags.single = 16, lags.multi = 5, multivariate.only = TRUE)
# The Jarque-Bera normality tests
normality.test(mod_VAR, multivariate.only = TRUE)
# Serial correlation: Portmanteau test and the Breusch-Godfrey LM test 
serial.test(mod_VAR, lags.pt = 16, lags.bg = 5,
            type = c("PT.asymptotic", "PT.adjusted", "BG", "ES"))
# Stability: Computes an empirical fluctuation process according to a specified method from the generalised fluctuation test framework
stability(mod_VAR, type = c("OLS-CUSUM", "Rec-CUSUM", "Rec-MOSUM",
                      "OLS-MOSUM", "RE", "ME", "Score-CUSUM", "Score-MOSUM", "fluctuation"),
          h = 0.15, dynamic = FALSE, rescale = TRUE)

## Forecasting with var package
# The n.ahead forecasts are computed recursively for the estimated VAR, beginning with h = 1, 2, …, n.ahead:
#   \bold{y}_{T+1 | T} = A_1 \bold{y}_T + … + A_p \bold{y}_{T+1-p} + C D_{T+1}
#The variance-covariance matrix of the forecast errors is a function of Σ_u and Φ_s.
mod_VAR_forecast <- predict(mod_VAR, n.ahead = 20, ci = 0.95, dumvar = NULL)
print(mod_VAR_forecast)
plot(mod_VAR_forecast)
mod_VAR_forecast_2 <- forecast(mod_VAR)





