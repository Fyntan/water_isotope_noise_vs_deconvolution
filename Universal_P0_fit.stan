//Universal P0 Fit - Power-law + white-noise

data {
  int<lower = 0> N;
  vector[N] freq;
  vector[N] spec;
  
  real alpha_mu;
  real<lower = 0> alpha_sigma;
  real beta_mu;
  real<lower = 0> beta_sigma;
  real epsilon_mu;
  real<lower = 0> epsilon_sigma;
  real phi_scale;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> beta;
  real<lower = 0> epsilon;
  real<lower = 0> phi;
}

transformed parameters {
  vector[N] y_hat;

  for(n in 1:N){
    y_hat[n] = alpha*freq[n]^-beta + epsilon;
  }
}

model {
  
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(beta_mu, beta_sigma);
  epsilon ~ normal(epsilon_mu, epsilon_sigma);
  
  phi ~ normal(0, phi_scale);
  
  spec ~ gamma(phi, phi ./ y_hat);
  
}
