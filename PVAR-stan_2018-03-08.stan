data {
  int<lower=0> nT; // No. of time periods observed
  int<lower=1> nP; // Lag order
  int<lower=0> nY; // No. of variables in y_t
  int Y[nT, nY]; // Obs of state variables
  int nX; // No. of covariates
  vector[nX] X[nT]; // Covariates
}
transformed data {
  vector[nY] Y_real[nT]; // declare count as array with vectors (real) for lhs evaluation
  int Y_obs[nT - nP, nY]; // declare Y observations where lagged values are available
  for(t in 1:nT) {
    Y_real[t] = to_vector(Y[t]);
    }
  for(t in nP+1:nT) {
    Y_obs[t-nP] = Y[t];
    }
}
parameters {
  matrix[nY, nY] B[nP];
  cholesky_factor_corr[nY] L_corr_noise;
  vector<lower=0>[nY] sd_noise;
  vector[nY] A; // constant of each Y
  matrix[nY, nX] C; // intercepts for X
}
transformed parameters {
  matrix[nY,nY] L_sigma;
  L_sigma = diag_pre_multiply(sd_noise, L_corr_noise);
}
model {
//  vector[nY] Y_lat[nT-nP]; // declare latent Y
  vector[nY] lambda[nT-nP]; // declare latent Y
  real Y_lat[nT - nP, nY]; // declare latent Y
  vector[nY] mus[nT-nP]; // declare mus_t
  for (t in 1:(nT-nP)) { // for each t
    mus[t] = A + C * X[t+nP]; // mu = A + C X
    for (p in 1:nP) // for each lag
      mus[t] += B[p] * Y_real[t+nP-p]; // mu = mu + B_p Y_t-p
  }
  L_corr_noise ~ lkj_corr_cholesky(2.0); // prior for correlation noise
  sd_noise ~ normal(0, 1); // prior for standard deviation noise
  A ~ normal(0, 1); // prior vor intercepts

  for (p in 1:nP)
    to_vector(B[p]) ~ normal(0, 1); // priors for B
  lambda ~ multi_normal_cholesky(mus,L_sigma);

  for(t in 1:nT) {
    Y_obs[t] ~ poisson_log(Y_lat[t]); // Y_lat is 2 dim array of real (equivalent to int)
    lambda[t] = to_vector(Y_lat[t]); // transform to array
    }
}
generated quantities {
  matrix[nY,nY] Sigma;
  Sigma = L_sigma * L_sigma';
}
