#include /include/trans_probs.stan

data {
  int smooth_cf;
  int smooth_inc;
  int remission; 
  int trend;
  int<lower=0> nage;
  int<lower=0> eqage;
  int<lower=0> mort_num[nage];
  int<lower=0> mort_denom[nage];
  int<lower=0> prev_num[nage];
  int<lower=0> prev_denom[nage];
  int<lower=0> inc_num[nage];
  int<lower=0> inc_denom[nage];
  int<lower=0> rem_num[nage];
  int<lower=0> rem_denom[nage];
  int<lower=0> nyr;
  
  // only in smoothed model 
  int<lower=0> K; // number of spline basis variables including the intercept
  matrix[nage,K] X;
  real<lower=0> sprior; 

  // alternative models
  int increasing_cf; // requires smooth_cf 
  int const_cf; // special case of increasing_cf 

  // Multiplier for incidence in each age and calendar year
  // first index is age, second is calendar year
  matrix<lower=0>[nage,nyr] inc_trend;
  matrix<lower=0>[nage,nyr] cf_trend;
}

parameters {
  vector<lower=0>[nage*(1-smooth_inc)] inc_par;
  vector<lower=0>[nage*(1-smooth_cf)] cf_par;
  vector[nage*remission] rem_par; // trick to exclude parameter if no remission
  vector[K*smooth_cf*(1-const_cf)] beta;
  vector<lower=0>[smooth_cf] lambda;
  vector[K*smooth_inc] beta_inc;
  vector<lower=0>[smooth_inc] lambda_inc;
  real<lower=0,upper=1> prevzero;

  vector[1*increasing_cf] lcfbase;

}

transformed parameters {
  //// range constraints on these cause problems due to floating point fuzz
  vector<lower=0>[nage] cf;
  vector<lower=0>[nage*increasing_cf] dcf;  // only in increasing model
  vector<lower=0>[nage] inc;
  vector<lower=0>[nage] inc_prob;
  vector[nage] rem;
  vector<lower=0>[nage] rem_prob;
  
  matrix[nage+1,3] state_probs; 
  row_vector[3] tmp;
  matrix[3,3] P;
  real prev[nage];
  real mort[nage];

  matrix<lower=0>[nage*trend,nyr*trend] cf_yr;
  matrix<lower=0>[nage*trend,nyr*trend] inc_yr;
  row_vector[3] state_probs_yr[(nage+1)*trend,nyr*trend];   
  
  /// Case fatality as smooth spline function of age
  /// Spline basis X passed from R
  if (smooth_inc) inc = exp(X*beta_inc); else inc = inc_par;
  if (remission) rem = rem_par; else rem = rep_vector(0, nage);
  rem_prob = 1 - exp(-rem);
  
  // Infer age zero prevalence from data if there are any data at age zero
  if (prev_denom[1]==0 || prev_num[1]==0)
    prev[1] = 0;
  else prev[1] = prevzero;

  if (increasing_cf) {
      // Baseline for eqage (e.g. age 50) is a random effect
      for (a in 1:(eqage-1)){
        cf[a] = exp(lcfbase[1]);
      }
      if (!const_cf){
	dcf = exp(X*beta);
      } else dcf = rep_vector(0, nage);
      for (a in eqage:nage){	
	cf[a] = cf[a-1] + dcf[a];
      }
  } else {
	if (smooth_cf) cf = exp(X*beta); else cf = cf_par;
  }

  if (trend) {
  // Define year-specific cf, inc, rem as function of year-indep versions 
    cf = exp(X*beta);
    for (b in 1:nyr){
      inc_yr[1:nage, b] = inc .* inc_trend[,b];  // note dot star .* for elementwise
      cf_yr[1:nage, b] = cf .* cf_trend[,b];
    }
    // state occupancy at age 0 (a=1)
    for (b in 1:nage) { 
      state_probs_yr[1,b,1] = 1;
      state_probs_yr[1,b,2] = 0;
      state_probs_yr[1,b,3] = 0;
    }
  } else { 
    state_probs[1,1] = 1;
    state_probs[1,2] = 0;
    state_probs[1,3] = 0;
  }

  for (a in 1:nage){
    if (trend) {
      if (a > 1) { 
	int y;
	for (b in 2:a){
	  y = nyr - a + b;  // y = nage-a+1 is birth.  y = nyr = nage is current year 
	  P = trans_probs(inc_yr[b-1, y-1], cf_yr[b-1, y-1], rem[b-1]);
	  tmp = state_probs_yr[b-1, y-1, 1:3] * P;
	  state_probs_yr[b, y, 1:3] = tmp;
	}
	if (a+1 <= nyr) { 
	  for (b in (a+1):nyr) {
	    state_probs_yr[b, y, 1:3] = rep_row_vector(0, 3);
	  }
	}
      }
      // data are the outcomes at the end of the current year
      P = trans_probs(inc_yr[a,nyr], cf_yr[a,nyr], rem[a]);
      inc_prob[a] = P[1,2] + P[1,3];
      prev[a] = state_probs_yr[a,nyr,2] /
	(state_probs_yr[a,nyr,1] + state_probs_yr[a,nyr,2]);
    } else { 
      P = trans_probs(inc[a], cf[a], rem[a]);
      inc_prob[a] = P[1,2] + P[1,3];
      if (a > 1)
	prev[a] = state_probs[a,2] / (state_probs[a,1] + state_probs[a,2]);
      tmp = state_probs[a,1:3] * P;  // temp variable to avoid warning
      state_probs[a+1,1:3] = tmp;
    }
    mort[a] = P[1,3]*(1 - prev[a]) + P[2,3]*prev[a];
    //// work around floating point fuzz
    if (mort[a] < 0) mort[a] = 0;
    if (mort[a] > 1) mort[a] = 1;
  }
}

model {
  mort_num ~ binomial(mort_denom, mort);
  inc_num ~ binomial(inc_denom, inc_prob);
  prev_num ~ binomial(prev_denom, prev);
  if (remission) rem_num ~ binomial(rem_denom, rem_prob);

  if (smooth_cf)  {
    if (!const_cf) {
      for (i in 1:(K-2)) {
	beta[i] ~ normal(0, lambda[1]);
      }
      for (i in (K-1):K){
	beta[i] ~ normal(0, 100);
      }
    }
    lambda[1] ~ exponential(sprior);
  }
  else {
    for (a in 1:nage){
      cf_par[a] ~ exponential(1);
    }
  }
  if (smooth_inc)  { 
    beta_inc[1] ~ normal(0, 100); 
    for (i in 1:(K-2)) {
      beta_inc[i] ~ normal(0, lambda_inc[1]);
    }
    for (i in (K-1):K){
      beta_inc[i] ~ normal(0, 100);
    }
    lambda_inc[1] ~ exponential(sprior);
  }
  else {
    for (a in 1:nage){
      inc_par[a] ~ exponential(1);
    }
  }
  
}

generated quantities {
  vector[nage] ll_mort;
  vector[nage] ll_inc;
  vector[nage] ll_prev;
  vector[nage*remission] ll_rem;
  for (a in 1:nage) {
      ll_mort[a] = binomial_lpmf(mort_num[a] | mort_denom[a], mort[a]);
      ll_inc[a] = binomial_lpmf(inc_num[a] | inc_denom[a], inc_prob[a]);
      ll_prev[a] = binomial_lpmf(prev_num[a] | prev_denom[a], prev[a]);
      if (remission) 
	  ll_rem[a] = binomial_lpmf(rem_num[a] | rem_denom[a], rem_prob[a]);
  }
}
