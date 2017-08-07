data {
  int<lower=0> N;       // number of non-missing values
  int<lower=0> m; // number of TMLs
  int<lower=0> n; // number of inspections
  
  vector<lower=0>[N] Y; // The thickness of the pipelines
  int<lower=0,upper=m> jj[N];// Store the Row
  int<lower=0,upper=n> kk[N];// Store the Column
  
  vector<lower=0>[n] t; //Time of inspections
}

parameters  {
  real muBeta;
  real<lower=0,upper=5> sigmaBeta;
  
  real<lower=0,upper=5> sigmaEpsilon;
  vector<lower=0>[m] beta;
}

model {
  muBeta ~ normal(0, 1000);
  beta ~ lognormal(muBeta, sigmaBeta);
  
  for(nn in 1:N){
    Y[nn] ~ normal(0.25 - beta[jj[nn]]*(t[kk[nn]] - 7347)/365, sigmaEpsilon);
  }
}
