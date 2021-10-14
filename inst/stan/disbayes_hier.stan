#include /include/trans_probs.stan

// "intercept" model: random intercepts, but identical slopes and deviations from linearity

data {
  int<lower=0> nage;
  int<lower=0> narea; 
  int<lower=0> ng;
  int<lower=0> eqage;
  int remission;
  int prev_zero;
  int<lower=0> mort_num[nage,narea,ng];
  int<lower=0> mort_denom[nage,narea,ng];
  int<lower=0> prev_num[nage,narea,ng];
  int<lower=0> prev_denom[nage,narea,ng];
  int<lower=0> inc_num[nage,narea,ng];
  int<lower=0> inc_denom[nage,narea,ng];
  int<lower=0> rem_num[nage,narea,ng];
  int<lower=0> rem_denom[nage,narea,ng];
  
  int<lower=0> K; // number of spline basis variables including the intercept
  matrix[nage,K] X;
  real<lower=0> sprior[2]; 
  real mipm;
  real<lower=0> mips;
  real mism;
  real<lower=0> miss;
  real<lower=0> gpint_a; 
  real<lower=0> gpint_b; 
  real<lower=0> gpslope_a; 
  real<lower=0> gpslope_b;
  real<lower=0> gender_int_priorsd;
  real<lower=0> gender_slope_priorsd;
  
  // alternative models
  int interceptonly;
  int increasing;
  int common;
  int const_cf; 
  int const_rem;
  int smooth_inc;

  // Empirical Bayes method where random effects hyperparameters are fixed
  int sd_int_isfixed;
  int sd_slope_isfixed;
  real<lower=0> sd_int_fixed;
  real<lower=0> sd_slope_fixed;
  real<lower=0> inc_prior[2]; 
  real<lower=0> rem_prior[2]; 

  // Empirical Bayes method where smoothing parameters are fixed
  int scf_isfixed;
  int scfmale_isfixed;
  int sinc_isfixed;
  real<lower=0> lambda_cf_fixed;
  real<lower=0> lambda_cf_male_fixed;
  real<lower=0> lambda_inc_fixed;
}

parameters {
  real<lower=0> inc_par[nage*(1 - smooth_inc),narea,ng];
  real<lower=0> rem_par[remission*(nage*(1-const_rem) + 1*const_rem),ng];

   // standard normal terms contributing to area-specific coefficients in non-centered parameterisation.
  matrix[(K-2)*(1-const_cf), narea*(1 - common) + 1*(common)] barea; 
  matrix[(1-interceptonly)*(1 - const_cf), narea*(1 - common) + 1*common] barea_slope; 
  matrix[1, narea*(1 - common)] barea_inter; 

  vector[K*(ng-1)] bmale;        // male effect on beta 
  
  // for model with increasing slopes
  vector[(narea*(1-common) + 1*common)*increasing] lcfbase;

  real beta_inc[K*smooth_inc,narea,ng];
  
  real mean_inter; // random effect mean intercept
  vector<lower=0>[1-sd_int_isfixed] sd_inter; 
  vector[1-const_cf] mean_slope;  
  vector<lower=0>[(1 - const_cf)*(1 - interceptonly)*(1 - increasing)*(1 - sd_slope_isfixed)] sd_slope;  // excluded if intercept-only, or increasing, or empirical Bayes
  vector<lower=0>[(1-const_cf)*(1-scf_isfixed)] lambda_cf;
  vector<lower=0>[(ng-1)*(1-scfmale_isfixed)] lambda_cf_male;
  vector<lower=0>[smooth_inc*(1-sinc_isfixed)] lambda_inc;
  real<lower=0,upper=1> prevzero[narea*prev_zero,ng];
}

transformed parameters {
  //// range constraints on these cause problems due to floating point fuzz
  real<lower=0> inc[nage,narea,ng];   // independent incidence for each area
  real<lower=0> cf[nage,narea,ng];
  real<lower=0> dcf[nage*increasing,narea,ng];  // only in increasing model
  real<lower=0> inc_prob[nage,narea,ng];
  real<lower=0> prev[nage,narea,ng];
  real<lower=0> mort[nage,narea,ng];
  real<lower=0> rem[nage,narea,ng];
  real<lower=0> rem_prob[nage,narea,ng];
  matrix[nage+1,3] state_probs; 
  row_vector[3] tmp;
  matrix[3,3] P;
  real<lower=0> sdint_use;
  real<lower=0> sdslope_use;
  
  matrix[K,narea] bareat;  // area-specific coefficients.
  real beta[K,narea,ng];
  real<lower=0> lambda_cf_use;
  real<lower=0> lambda_cf_male_use;
  real<lower=0> lambda_inc_use;
  vector[narea*increasing] lcfbase_use;

  if (sd_int_isfixed) sdint_use = sd_int_fixed; else sdint_use = sd_inter[1];
  if (sd_slope_isfixed || const_cf || interceptonly || increasing) sdslope_use = sd_slope_fixed; else sdslope_use = sd_slope[1];
  if (scf_isfixed || const_cf) lambda_cf_use = lambda_cf_fixed; else lambda_cf_use = lambda_cf[1];
  if (scfmale_isfixed || (ng==1)) lambda_cf_male_use = lambda_cf_male_fixed; else lambda_cf_male_use = lambda_cf_male[1];
  if (sinc_isfixed || !smooth_inc) lambda_inc_use = lambda_inc_fixed; else lambda_inc_use = lambda_inc[1];

  for (j in 1:narea){
    if (common) {  // no difference between areas. implemented to allow statistical model comparison
      if (increasing) { lcfbase_use[j] = lcfbase[1]; }
      if (const_cf){
	for (i in 1:(K-1)) {
	  bareat[i,j] = 0; // constant 
	}
      } else { 
	for (i in 1:(K-2)){
	  bareat[i,j] = barea[i,1] * lambda_cf_use;
	}
	if (increasing) { // increasing and smooth
	  bareat[K-1,j] = barea_slope[1,1] * lambda_cf_use; // slope for increments. shrunk. 
	} else { // unconstrained and smooth
	  bareat[K-1,j] = mean_slope[1];
	}
      }
      if (increasing) {
	bareat[K,j] = mean_slope[1]; // common slope (ie intercept for increments)
      } else { 
	bareat[K,j] = mean_inter; // area-level random intercept
      }
    }

    else {  // area-specific terms 
      if (increasing) { lcfbase_use[j] = lcfbase[j]; }
      if (const_cf){
	for (i in 1:(K-1)) {
	  bareat[i,j] = 0;
	}
      } else { 
	for (i in 1:(K-2)) {
	    bareat[i,j] = barea[i,j] * lambda_cf_use; // smoothing terms shared between areas
	}
	if (interceptonly) {
	    bareat[K-1,j] = mean_slope[1];   // common slope between areas 
	} else if (increasing) {
	  bareat[K-1,j] = barea_slope[1,j] * lambda_cf_use; // slope for increments. shrunk.
	} else { // default
	  bareat[K-1,j] = mean_slope[1] + barea_slope[1,j] * sdslope_use; // area-level random slope 
	}
      }
      if (increasing) {
	bareat[K,j] = mean_slope[1]; // common slope (ie intercept for increments)
      } else{ 
	bareat[K,j] = mean_inter + barea_inter[1, j] * sdint_use; // area-level random intercept
      }
    }

  }
    
  for (g in 1:ng)  {
    for (j in 1:narea) { 
	for (a in 1:nage){ 
	  if (smooth_inc) {
	    inc[a,j,g] = exp(X[a,]*to_vector(beta_inc[,j,g]));
	  } else { 
	    inc[a,j,g] = inc_par[a,j,g];
	  }
	}
	for (k in 1:K) {
	  if (ng > 1) {
	    beta[k,j,g] = bareat[k,j] + bmale[k]*(g-1); // additive gender and area effects
	  } else {
	    beta[k,j,g] = bareat[k,j];
	  }	    
      }

      // Infer age zero prevalence from data if there are any data at age zero, or if we asked it to
      if (prev_denom[1,j,g] > 0 && (prev_num[1,j,g] > 0 || prev_zero))
	prev[1,j,g] = prevzero[j,g];
      else prev[1,j,g] = 0; 
      state_probs[1,1] = 1;
      state_probs[1,2] = 0;
      state_probs[1,3] = 0;
    
      if (increasing){
	/// Annual increments in case fatality as smooth spline function of age
	for (a in 1:nage) { 
	  dcf[a,j,g] = exp(X[a,]*to_vector(beta[,j,g]));
	}
	// Baseline for eqage (e.g. age 50) is a random effect
	for (a in 1:(eqage-1)){
	  cf[a,j,g] = exp(lcfbase_use[j]);
	}
	for (a in eqage:nage){
	  cf[a,j,g] = cf[a-1,j,g] + dcf[a,j,g];
	}      
      } 
      else {
	/// Case fatality as smooth spline function of age
	/// Spline basis X is 
	/// nage x K matrix  *  K colvector =  nage  vector
	for (a in 1:nage){ 
	  cf[a,j,g] = exp(X[a,]*to_vector(beta[,j,g]));
	}
      }
	
      if (remission) {
	if (const_rem) {
	  for (a in 1:nage)
	    rem[a,j,g] = rem_par[1,g];
	} else 
	  rem[,j,g] = rem_par[,g];
      } else { 	for (a in 1:nage) { rem[a,j,g] = 0; } }
      for (a in 1:nage){
	P = trans_probs(inc[a,j,g], cf[a,j,g], rem[a,j,g]);
	inc_prob[a,j,g] = P[1,2] + P[1,3];
	rem_prob[a,j,g] = P[2,1];
	if (a > 1)
	  prev[a,j,g] = state_probs[a,2] / (state_probs[a,1] + state_probs[a,2]);
	tmp = state_probs[a,1:3] * P;  // temp variable to avoid warning
	state_probs[a+1,1:3] = tmp;
	mort[a,j,g] = P[1,3]*(1 - prev[a,j,g]) + P[2,3]*prev[a,j,g];
	//// work around floating point fuzz
	if (mort[a,j,g] < 0) mort[a,j,g] = 0;
	if (mort[a,j,g] > 1) mort[a,j,g] = 1;
      }
    }
  }
}

model {    
  mean_inter ~ normal(mipm, mips);

  // These all get transformed in different ways according to the model 
  if (!const_cf) { 
    mean_slope ~ normal(mism, miss); 
    if (common) {
      for (i in 1:(K-2)){ 
	barea[i,1] ~ normal(0, 1);
      }     
      if (!interceptonly && !const_cf){
	barea_slope[1,1] ~ normal(0, 1);
      }
    } else {
      for (j in 1:narea){
	for (i in 1:(K-2)){ 
	  barea[i,j] ~ normal(0, 1);
	}     
	barea_inter[1,j] ~ normal(0, 1);
	if (!interceptonly && !const_cf){
	  barea_slope[1,j] ~ normal(0, 1);
	}
      }
    }
  }

  if (smooth_inc) {
    // Random effects model or additive gender effects not used for incidence.
    // Just have independent smooth curves for each area 
    // Common smoothness variance 
    for (j in 1:narea){
      for (g in 1:ng){
	for (i in 1:(K-2)) {
	  beta_inc[i,j,g] ~ normal(0, lambda_inc_use);
	}
	for (i in (K-1):K){
	  beta_inc[i,j,g] ~ normal(0, 100);
	}
      }
    }
  }

  // Model for the data
  if (increasing && common) {
    lcfbase[1] ~ normal(0, 100);
  }
  for (j in 1:narea){
    for (g in 1:ng) {
      mort_num[,j,g] ~ binomial(mort_denom[,j,g], mort[,j,g]);
      inc_num[,j,g] ~ binomial(inc_denom[,j,g], inc_prob[,j,g]);
      prev_num[,j,g] ~ binomial(prev_denom[,j,g], prev[,j,g]);
      if (remission) {
	rem_num[,j,g] ~ binomial(rem_denom[,j,g], rem_prob[,j,g]);
      }
    }
    if (increasing && !common){
      lcfbase[j] ~ normal(mean_inter, sd_inter);
    }
  }
  if (remission){
    if (const_rem) {
      for (g in 1:ng) rem_par[1,g] ~ gamma(rem_prior[1], rem_prior[2]);
    } else {
      for (g in 1:ng) {
	for (a in 1:nage) {
	  rem_par[a,g] ~ gamma(rem_prior[1], rem_prior[2]);
	}
      }
    }
  }
  
    // Degree of smoothness
  if (!const_cf && !scf_isfixed){
    lambda_cf ~ gamma(2, sprior[2]);
  }
  if (smooth_inc && !sinc_isfixed){
    lambda_inc ~ gamma(2, sprior[1]);
  }
  if ((!interceptonly) && (!increasing) && (!const_cf)){
	// Variation in slopes 
      sd_slope ~ gamma(gpslope_a, gpslope_b);
    }
    // Variation in intercepts
  sd_inter ~ gamma(gpint_a, gpint_b);
  
  if (ng > 1) {
    if (!const_cf && !scfmale_isfixed){
      lambda_cf_male ~ gamma(2, sprior[2]);
    }
    // effect of being male, assumed common between areas
    //... on spline coefficients governing deviation from linearity
    for (i in 1:(K-2)){
      bmale[i] ~ normal(0, lambda_cf_male_use); // smoothness variance par shared 
    }
    // ...on slopes 
    bmale[K-1] ~ normal(0, gender_slope_priorsd);
    // ...on intercepts
    bmale[K] ~ normal(0, gender_int_priorsd); 
    // five fold difference between men and women is log(5)=1.6 on log scale 
    // want sd such that 95% prob that bmale is between -log(5) and log(5)
    // sd = log(5) / qnorm(0.975) = 0.82 
  }

  if (!smooth_inc){
    for (a in 1:nage){
      for (j in 1:narea){
	for (g in 1:ng){
	  inc_par[a,j,g] ~ gamma(inc_prior[1], inc_prior[2]);
	}
      }
    }
  }

  if (prev_zero){
    for (j in 1:narea){
      for (g in 1:ng){
	prevzero[j,g] ~ beta(2,2); // boundary-avoiding
      }
    }
  }
}

generated quantities {
  vector[nage*narea*ng] ll_mort;
  vector[nage*narea*ng] ll_inc;
  vector[nage*narea*ng] ll_prev;
  vector[nage*narea*ng*remission] ll_rem;
  int i = 1;
  for (a in 1:nage) {
    for (j in 1:narea) {
      for (g in 1:ng) {
	ll_mort[i] = binomial_lpmf(mort_num[a,j,g] | mort_denom[a,j,g], mort[a,j,g]);
	ll_inc[i] = binomial_lpmf(inc_num[a,j,g] | inc_denom[a,j,g], inc_prob[a,j,g]);
	ll_prev[i] = binomial_lpmf(prev_num[a,j,g] | prev_denom[a,j,g], prev[a,j,g]);
	if (remission) 
	  ll_rem[i] = binomial_lpmf(rem_num[a,j,g] | rem_denom[a,j,g], rem_prob[a,j,g]);
	i = i + 1;
      }
    }
  }
}
