
data {
  int<lower=0> t;
  int y[t]; // observed cases
  real <lower = 0> inc_loc;
  real <lower = 0> inc_scale;
  real <lower = 0> delay_alpha;
  real <lower = 0> delay_beta;
  int tau; // window size
}

transformed data {
  vector[(2 * t) - 1] delayvec; // vector of discretised delay distribution
  vector[(2 * t) - 1] incvec; // vector of discretised incubation period
  matrix[(2 * t) - 1, (2 * t) - 1] incmat; // incubation period matrix (different convolutions)
  vector[(2 * t) - 1] inc_conv_delay; // convolution of 1 incubation period + 1 delay
  int indw; // index variables

  for(i in 1:((2 * t) - 1)){

    // Discretised Incubation Period distribution
    incvec[i] = (i < t) ? 0 : lognormal_cdf(i - t + 1, inc_loc, inc_scale) -
    lognormal_cdf(i - t, inc_loc, inc_scale);

    // Discretised Notification Delay distribution
    delayvec[i] = (i >= t) ? gamma_cdf(i - t + 1, delay_alpha, delay_beta) -
    gamma_cdf(i - t , delay_alpha, delay_beta) : 0;
  }

  for(i in 1:((2 * t) - 1)) {
    for(j in 1:((2 * t) - 1)) {
      indw = t + i - j;
      // Matrix for Incubation Period
      incmat[i, j] = (indw <= 0) ? 0: ((indw >= ((2 * t) -1)) ? 0 : incvec[indw]);
    }
  }

  // Incubation Period + Notification Delay convolution
  inc_conv_delay = incmat * delayvec;
}

parameters {
  real <lower = 0> lambda[t];
  real <lower = 0> epsilon;
}

transformed parameters {
  real mu[t];

  for(i in 1:t){
    mu[i] = 0;
    for(j in 1:i) {
      mu[i] += lambda[j] * inc_conv_delay[i - j + t];
    }
  }
}

model {
  epsilon ~ normal(0, 1) T[0, ];

  // for(i in 1:t) {
  //   lambda[i] ~ normal(0, 1 / sqrt(epsilon)) T[0, ];
  // }

  // for(s in (tau + 1):t) {
  //   for (i in (s-tau + 1):s) {
  //     target += poisson_lpmf(y[i] | mu[s]);
  //   }
  // }
  lambda[1] ~ normal(0, 10);
  for(k in 2:t) {
    lambda[k] ~ normal(lambda[k - 1], 1 / sqrt(epsilon)) T[0, ];
  }
  for(i in 1:t){
    y[i] ~ poisson(mu[i]);
  }
}

generated quantities {
  vector[(2 * t) - 1] icd;
  real muout[t];
  int infections[t];
  icd = inc_conv_delay;
  muout = mu;
  infections = poisson_rng(lambda);
}

