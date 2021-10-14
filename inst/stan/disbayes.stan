#include /include/trans_probs.stan

data {
  int smooth_cf;
  int smooth_inc;
  int remission; 
  int trend;
  int prev_zero;
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
  real<lower=0> sprior[2]; 

  // alternative models
  int increasing_cf; // requires smooth_cf 
  int const_cf; // special case of increasing_cf 
  int const_rem;

  int nbias; // 1 if no incidence bias modelled
  // which of the two alternative incidence values generates the incidence data 
  int incdata_ind;
  int prevdata_ind;  // .. and the prevalence data 
  
  // Multiplier for incidence in each age and calendar year
  // first index is age, second is calendar year
  matrix<lower=0>[nage,nyr] inc_trend;
  matrix<lower=0>[nage,nyr] cf_trend;

  real<lower=0> inc_prior[2]; 
  real<lower=0> cf_prior[2]; 
  real<lower=0> rem_prior[2]; 

  // Empirical Bayes method where smoothing parameters are fixed
  int scf_isfixed;
  int sinc_isfixed;
  real<lower=0> lambda_cf_fixed;
  real<lower=0> lambda_inc_fixed;
}

parameters {
  vector<lower=0>[nage*(1-smooth_inc)] inc_par;
  vector<lower=0>[nage*(1-smooth_cf)] cf_par;
  vector<lower=0>[remission*(nage*(1 - const_rem) + 1*const_rem)] rem_par;
  vector[K*smooth_cf*(1-const_cf)] beta;
  vector<lower=0>[smooth_cf*(1-scf_isfixed)] lambda_cf;
  vector<lower=0>[smooth_inc*(1-sinc_isfixed)] lambda_inc;
  vector[K*smooth_inc] beta_inc;
  vector<lower=0,upper=1>[prev_zero] prevzero;

  // Log HR between alternative incidence values in bias model, assumed common between ages
  vector[nbias==2] bias_loghr;

  vector<lower=0>[1*increasing_cf] cfbase;
}

transformed parameters {
  //// range constraints on these cause problems due to floating point fuzz
  vector<lower=0>[nage] cf;
  vector<lower=0>[nage*increasing_cf] dcf;  // only in increasing model
  matrix<lower=0>[nage,nbias] inc;
  matrix<lower=0,upper=1>[nage,nbias] inc_prob;
  vector<lower=0>[nage] rem;
  vector<lower=0,upper=1>[nage] rem_prob;
  
  row_vector[3] state_probs[(nage+1)*(1-trend),nbias]; 
  row_vector[3] tmp;
  matrix<lower=0,upper=1>[3,3] P;
  matrix<lower=0>[nage,nbias] prev;
  real<lower=0> mort[nage];

  matrix<lower=0>[nage*trend,nyr*trend] cf_yr;
  real<lower=0> inc_yr[nage*trend,nyr*trend,nbias];
  row_vector[3] state_probs_yr[(nage+1)*trend,nyr*trend,nbias];   

  real<lower=0> lambda_cf_use;
  real<lower=0> lambda_inc_use;
  if (scf_isfixed || !smooth_cf) lambda_cf_use = lambda_cf_fixed; else lambda_cf_use = lambda_cf[1];
  if (sinc_isfixed || !smooth_inc) lambda_inc_use = lambda_inc_fixed; else lambda_inc_use = lambda_inc[1];
  
  /// Case fatality as smooth spline function of age
  /// Spline basis X passed from R
  if (smooth_inc) inc[1:nage,1] = exp(X*beta_inc); else inc[1:nage,1] = inc_par;
  if (nbias > 1) { 
    inc[1:nage,2] = exp(log(inc[1:nage,1]) + bias_loghr[1]);
  }
  if (remission) {
    if (const_rem) {
      for (a in 1:nage)
	rem[a] = rem_par[1];
    }
    else rem = rem_par;
  } else rem = rep_vector(0, nage);
  
  // Infer age zero prevalence from data if there are any data at age zero, or if we asked it to
  for (k in 1:nbias) { 
    if (prev_denom[1] > 0 && (prev_num[1] > 0 || prev_zero))
      prev[1,k] = prevzero[1];
    else prev[1,k] = 0; 
  }
  
  if (increasing_cf) {
      // Baseline for eqage (e.g. age 50) is a random effect
      for (a in 1:(eqage-1)){
        cf[a] = cfbase[1];
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
      cf_yr[1:nage, b] = cf .* cf_trend[,b];
    }
    for (k in 1:nbias){ 
      for (b in 1:nyr){
	for (a in 1:nage){
	  inc_yr[a, b, k] = inc[a,k] * inc_trend[a,b];
	}
      }
      // state occupancy at age 0 (a=1)
      for (b in 1:nyr) { 
	state_probs_yr[1,b,k,1] = 1;
	state_probs_yr[1,b,k,2] = 0;
	state_probs_yr[1,b,k,3] = 0;
	// initialise state occupancy at other ages to keep Stan happy that all array elements are initialised
	// only the upper diagonal of this (year >= age) is needed
	for (a in 2:(nage+1)){
	  state_probs_yr[a,b,k,1:3] = rep_row_vector(0, 3);
	}
      }
    }
  } else {
    for (k in 1:nbias){ 
      state_probs[1,k,1] = 1;
      state_probs[1,k,2] = 0;
      state_probs[1,k,3] = 0;
    }
  }

  for (a in 1:nage){
    if (trend) {
      for (k in 1:nbias){
	if (a > 1) { 
	  int y;
	  for (b in 2:a){
	    y = nyr - a + b;  // y = nage-a+1 is birth.  y = nyr = nage is current year 
	    P = trans_probs(inc_yr[b-1, y-1, k], cf_yr[b-1, y-1], rem[b-1]);
	    tmp = state_probs_yr[b-1, y-1, k, 1:3] * P;
	    state_probs_yr[b, y, k, 1:3] = tmp;
	  }
	}
	// data are the outcomes at the end of the current year
	P = trans_probs(inc_yr[a,nyr,k], cf_yr[a,nyr], rem[a]);
	inc_prob[a,k] = P[1,2] + P[1,3];
	prev[a,k] = state_probs_yr[a,nyr,k,2] /
	  (state_probs_yr[a,nyr,k,1] + state_probs_yr[a,nyr,k,2]);
	if (k==1) mort[a] = P[1,3]*(1 - prev[a,1]) + P[2,3]*prev[a,1];
      }
    } else { 

      for (k in 1:nbias){ 
	P = trans_probs(inc[a,k], cf[a], rem[a]);
	inc_prob[a,k] = P[1,2] + P[1,3];
	if (a > 1)
	  prev[a,k] = state_probs[a,k,2] / (state_probs[a,k,1] + state_probs[a,k,2]);
	tmp = state_probs[a,k,1:3] * P;  // temp variable to avoid warning
	state_probs[a+1,k,1:3] = tmp;
	if (k==1) {
	  mort[a] = P[1,3]*(1 - prev[a,1]) + P[2,3]*prev[a,1];
	  rem_prob[a] = P[2,1];
	}
      }      

    }
    //// work around floating point fuzz
    if (mort[a] < 0) mort[a] = 0;
    if (mort[a] > 1) mort[a] = 1;
  }
}

model {
  mort_num ~ binomial(mort_denom, mort);
  inc_num ~ binomial(inc_denom, inc_prob[1:nage,incdata_ind]);
  prev_num ~ binomial(prev_denom, prev[1:nage,prevdata_ind]);
  if (remission) {
    rem_num ~ binomial(rem_denom, rem_prob);
    if (const_rem) rem_par[1] ~ gamma(rem_prior[1], rem_prior[2]);
    else { 
      for (a in 1:nage){
	rem_par[a] ~ gamma(rem_prior[1], rem_prior[2]);
      }
    }
  }
  
  if (smooth_cf)  {
    if (!const_cf) {
      for (i in 1:(K-2)) {
	beta[i] ~ normal(0, lambda_cf_use);
      }
      for (i in (K-1):K){
	beta[i] ~ normal(0, 100);
      }
    }
    if (!scf_isfixed) 
      lambda_cf[1] ~ gamma(2, sprior[2]);
  }
  else {
    for (a in 1:nage){
      cf_par[a] ~ gamma(cf_prior[1], cf_prior[2]); // boundary-avoiding with mode 2
    }
  }
  if (increasing_cf) {
    cfbase[1] ~ gamma(cf_prior[1], cf_prior[2]);
  }
  if (smooth_inc)  { 
    for (i in 1:(K-2)) {
      beta_inc[i] ~ normal(0, lambda_inc_use);
    }
    for (i in (K-1):K){
      beta_inc[i] ~ normal(0, 100);
    }
    if (!sinc_isfixed) 
      lambda_inc[1] ~ gamma(2, sprior[1]); 
  }
  else {
    for (a in 1:nage){
      inc_par[a] ~ gamma(inc_prior[1], inc_prior[2]);
    }
  }
  
  if (nbias==2){
    bias_loghr ~ normal(0,1);
  }

  if (prev_zero){
    prevzero[1] ~ beta(2,2); // boundary-avoiding
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
