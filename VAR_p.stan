
data {
	int<lower=0> nT; // Number of time periods observed
	int<lower=0> K; // Number of time series
	int<lower=1> P; // Lag order
	vector[K] Z[nT]; // k-dimensional multivariate time series
	int<lower=0> nF; // Number of time periods to forecast
}
transformed data {
	vector[K] Z_t[nT-P]; // 
	for (t in 1:(nT-P)) 
		Z_t[t] = Z[t + P];
	
}
parameters {
	matrix[K, K] Phi_mat[P]; // Estimated coefficient matrix
	cholesky_factor_corr[K] L_corr_noise; // 
	vector<lower=0>[K] sd_noise;
	vector[K] Phi_0;
}
transformed parameters {
	matrix[K,K] L_sigma;
	// Return the product of the diagonal matrix form the vector sd_noise
	// and the matrix L_corr_noise
  L_sigma = diag_pre_multiply(sd_noise, L_corr_noise);
}
model {
	vector[K] mus[nT-P]; 
	
	// Priors
	L_corr_noise ~ lkj_corr_cholesky(2.0);
	sd_noise ~ normal(0,1);
	Phi_0 ~ normal(0,1);
	for (p in 1:P)
		to_vector(Phi_mat[p]) ~ normal(0, 1);
		
	// k-dimensional mean vector
	for (t in 1:(nT-P)) {
		mus[t] = Phi_0;
		for (p in 1:P) 
			mus[t] = mus[t] + Phi_mat[p] * Z[t+p-1];
	}
	
	// Likelihood
	Z_t ~ multi_normal_cholesky(mus,L_sigma);
}
generated quantities {
	matrix[K,K] Sigma; // Covariance matrix
	//vector[K] Z_f[nF]; // Forecasting vector
	vector[K] Z_f_temp[P + nF];
	vector[K] mus_f[nF]; 
	int count1;
  int count2;
	
	// Covariance matrix
	Sigma = L_sigma * L_sigma'; 
	
	// 
  count1 = 0;
  for(t in (nT-P):(nT)){
    count1 = count1 + 1;
    Z_f_temp[count1] = Z[t];
  }
  
  count2 = P;
  for (f in 1:(nF)) {
    count2 = count2 + 1;
		mus_f[f] = Phi_0;
		for (p in 1:P) 
		  mus_f[f] = mus_f[f] + Phi_mat[p] * Z_f_temp[f+p-1];
	
	Z_f_temp[count2] = multi_normal_cholesky_rng(mus_f[f],L_sigma);
	}
    }
