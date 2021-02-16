functions{

  real omega_upper_tri(int a, int b, matrix omega){
     if(a<=b)
        return omega[a,b];
     else
        return omega[b,a];
  }

}

data {
  int<lower=1> n;
  int<lower=0> X[n, n];
  int<lower=1> k;
  int C[n,2];
  real<lower=0> rates_std_prior[2];
}

parameters {
 
  positive_ordered[2] mu;
  positive_ordered[2] phi;
  matrix<lower=0,upper=1>[k,k] omega; 
  
}

model {
  
  mu[1] ~ normal(1, rates_std_prior[1]);
  mu[2] ~ normal(1, rates_std_prior[2]);
  phi[1] ~ normal(1, rates_std_prior[1]);
  phi[2] ~ normal(1, rates_std_prior[1]); 

  for(i in 1:k){
     for(j in 1:k){
        
        omega[i,j] ~ beta(1,1);
     }
  }

   for (i in 1:n) {
    for (j in i + 1:n) {
      real log_mu_ij_0 = neg_binomial_2_lpmf(X[i, j] | mu[1], phi[1]);
      real log_mu_ij_1 = neg_binomial_2_lpmf(X[i, j] | mu[2], phi[2]);
      
     
     real     log_nu_ij_0 = bernoulli_lpmf(0 |  omega_upper_tri(C[i,2],C[j,2], omega));
     real     log_nu_ij_1 = bernoulli_lpmf(1 |  omega_upper_tri(C[i,2],C[j,2], omega));
     
      
      real z_ij_0 = log_mu_ij_0 + log_nu_ij_0;
      real z_ij_1 = log_mu_ij_1 + log_nu_ij_1;
      if (z_ij_0 > z_ij_1) {target += z_ij_0 + log1p_exp(z_ij_1 - z_ij_0);}
      else {target += z_ij_1 + log1p_exp(z_ij_0 - z_ij_1);}
    }
  }
}

generated quantities {

  real Q[n ,n];
  
  for (i in 1:n) {
    Q[i, i] = 0;
    for (j in i+1:n) {
      real log_mu_ij_0 = neg_binomial_2_lpmf(X[i, j] | mu[1], phi[1]);
      real log_mu_ij_1 = neg_binomial_2_lpmf(X[i, j] | mu[2], phi[2]);
     
      
       real   log_nu_ij_0 = bernoulli_lpmf(0 |  omega_upper_tri(C[i,2],C[j,2], omega));
       real   log_nu_ij_1 = bernoulli_lpmf(1 |  omega_upper_tri(C[i,2],C[j,2], omega));
      

      real z_ij_0 = log_mu_ij_0 + log_nu_ij_0;
      real z_ij_1 = log_mu_ij_1 + log_nu_ij_1;
      Q[i, j] = 1 / (1  + exp(z_ij_0 - z_ij_1));
      Q[j, i] = Q[i, j];
    }
  }
}


