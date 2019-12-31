# Helper functions

# https://faculty.washington.edu/ezivot/book/Ch11.VARexamples2ndEdition.ssc
# simulateVAR
sim_VAR <- function(nobs, ) {
  # number of observations
  nobs = 250
  # n x n coefficient matrix
  pi1 = matrix(c(0.7,0.2,0.2,0.7),2,2)
  # mean of time series
  mu.vec = c(1,5)
  # Covariance matrix
  sigma <- matrix(c(1,0.5,0.5,1), ncol=2)
  # unorservable zero mean white noise vector process
  e.var <- rmvnorm(n=nobs, mean = rep(0, nrow(sigma)), sigma=sigma)
  e.var = t(e.var)
  # Intercepts for time series
  c.vec = as.vector((diag(2)-pi1)%*%mu.vec)
  # Simulate VAR process
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
} 


genstaVAR <- function(k = k,             # Number time-series
                      p = p,             # Order of VAR process
                      nT = nT,           # number of time points
                      Sigma = NULL,      # Covariance matrix
                      eigenvalues = NULL   # Eigenvalues of multi-companion matrix
) {
  
  # Checks
  if(is.null(k) | is.null(p) | is.null(nT)){
    stop("k, p and nt must be specify")
  }
  
  # Positive-definite covariance matrix
  if(is.null(Sigma)) {
    Sigma <- genPositiveDefMat(dim = k, covMethod=c("eigen"))
  }
  # Random eigenvalues
  if(is.null(eigenvalues)) {
    eigenvalues <- runif(k, min = 0, max = 0.99)
  }
  # Multi-companion matrix with positive eigenvalues
  mc_matrix <- sim_mc(5,k, eigval = eigenvalues)
  
  # Coefficient matrix
  B <- matrix(0,nrow=k,ncol=p*k)
  B[,1:k] <- mc_matrix$mat
  Phi <- VarptoVar1MC(B,p,k)
  
  # Simulate VAR
  Z <-MultVarSim(k = k,
                 A1 = Phi,
                 p = p,
                 Sigma = Sigma$Sigma,
                 nT)
  
  return(list(Z = Z, 
              Sigma = Sigma$Sigma, 
              Phi = Phi))
}





