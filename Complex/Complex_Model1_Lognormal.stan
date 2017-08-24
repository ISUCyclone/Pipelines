data {
  int<lower=0> N;       // number of the TMLs
  int<lower=0> K;       // number of inspections
  int<lower=0> N_elbow;  // number of elbows
  int<lower=0> N_pipe;   // number of pipes
  matrix[N,K] y;
  vector<lower=0>[K] t;
}

parameters  {
  real<lower=0.4,upper=0.47> muElbow; //By default, Uniform
  real<lower=0.4,upper=0.47> muPipe;
  real<lower=0.5,upper=0.62> muTee;
  real<lower=0.00001,upper=5> sigmaY0;
  real<lower=0.000001,upper=0.022> medianBeta;
  real<lower=0.00001,upper=5> sigmaBeta;
  real<lower=9.31,upper=12> muT;
  real<lower=0.00001,upper=10> sigmaT;
  real<lower=0.00001,upper=0.25> sigmaEpsilon;
  
  vector<lower=0>[N] y_init;
  vector<lower=0>[N] beta;
  vector<lower=11048>[N] t_init;
}

transformed parameters {
  real mubeta;
  mubeta = log(medianBeta);
}

model {
//***************************************
// Initial Thickness
  for (n in 1:N_elbow) {
    y_init[n] ~ normal(muElbow,sigmaY0);
  }
  for (n in (N_elbow + 1):(N_elbow+N_pipe)) {
    y_init[n] ~ normal(muPipe,sigmaY0);
  }
  for (n in (N_elbow+N_pipe+1):N) {
    y_init[n] ~ normal(muTee,sigmaY0);
  }

//***************************************
// Beta(Corrosion Rate)
  beta ~lognormal(mubeta,sigmaBeta);

//***************************************  
// Ti(Initial Time)
  t_init ~ lognormal(muT, sigmaT);

//***************************************
// Model
  for (n in 1:N) {
    for (k in 1:K){
      y[n,k] ~ normal(y_init[n]-
      step(t[k] - t_init[n]) * beta[n]*(t[k]-t_init[n])/365,
      sigmaEpsilon);
    }
  }
}
