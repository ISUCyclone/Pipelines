data {
  int<lower=0> N; // number of TMLs
  int<lower=0> K; // number of inspectations
  int<lower=0> n; // number of non-missing
  
  vector<lower=0>[n] y; // thickness of the pipes
  int<lower=0,upper=N> Rowy[n];
  int<lower=0,upper=K> Coly[n];
  
  vector<lower=0>[K] T;
  vector<lower=0,upper=1>[N] q1;
  vector<lower=0,upper=1>[N] q2;
  vector<lower=0,upper=1>[N] q3;
  vector<lower=0,upper=1>[N] q4;
}

parameters {
  real mubeta1;
  real mubeta2;
  real mubeta3;
  real mubeta4;
  
  real<lower=0,upper=1000> sigmabeta;
  real<lower=0,upper=1000> sigma;
  
  vector<lower=0>[N] beta;
}

transformed parameters  {
  vector[N] mu;
  for(i in 1:N) {
    mu[i] = q1[i] * mubeta1 +q2[i] * mubeta2 + q3[i] * mubeta3 + q4[i] * mubeta4;
  }
}

model {
  mubeta1 ~ normal(0,1000);
  mubeta2 ~ normal(0,1000);
  mubeta3 ~ normal(0,1000);
  mubeta4 ~ normal(0,1000);
  
  for(i in 1:N){
    beta[i] ~ lognormal(mu[i],sigmabeta);
  }
  
  for(i in 1:n){
    y[i] ~ normal(0.25-beta[Rowy[i]]*(T[Coly[i]] - 7347)/365, sigma);
  }
}
