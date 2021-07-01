#include /include/trans_probs.stan

// "intercept" model: random intercepts, but identical slopes and deviations from linearity

data {
  int<lower=0> nage;
  int<lower=0> narea; 
  int<lower=0> ng;
  int<lower=0> eqage;
  int remission;
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
  real<lower=0> sprior; 
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
}

// CHANGES IN GENDER ADAPTATION
// extend dim of data arrays.  add gender constant
// rename beta to barea 
// extend dim of bmale to matrix to allow g>2 
// seven transformed pars changed form matrix to real array 
// tp, loop over areas: rename betat to bareat 
// tp : add a g in loop
// define beta=barea, then if statement to only add on bmale   
// then rest of transformed pars is the same 

parameters {
  real<lower=0> inc[nage,narea,ng];   // independent incidence for each area
  real rem_par[nage*remission,narea,ng]; // trick to exclude parameter if no remission
  
  matrix[K,narea] barea;  // standard normal terms contributing to area-specific coefficients in non-centered parameterisation.

  vector[K*(ng-1)] bmale;        // male effect on beta 
  
  // for model with increasing slopes
  vector[narea*increasing] lcfbase;
  
  real mean_inter; // random effect mean intercept
  real sd_inter; 
  real mean_slope;  
  vector[(1 - interceptonly)*(1 - increasing)] sd_slope;  // excluded if intercept-only or increasing
  real lambda_smooth;      // don't bound these to avoid unnecessary validation. priors already defined to be positive. 
  vector[ng-1] lambda_smooth_male;
  real<lower=0,upper=1> prevzero[narea,ng];
}

transformed parameters {
  //// range constraints on these cause problems due to floating point fuzz
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

  matrix[K,narea] bareat;  // area-specific coefficients.
  real beta[K,narea,ng];
    
    for (j in 1:narea){
	for (i in 1:(K-2)) {
	    bareat[i,j] = barea[i,j] * lambda_smooth; // smoothing terms shared between areas
	}

	if (interceptonly) {
	    bareat[K-1,j] = mean_slope;   // common slope between areas 
	} else if (increasing) {
	    bareat[K-1,j] = barea[K-1,j] * lambda_smooth;
	} else { // default
	    bareat[K-1,j] = mean_slope + barea[K-1, j] * sd_slope[1]; // area-level random slope 
	}

	if (increasing) {
	    bareat[K,j] = mean_slope; // common slope (ie intercept for increments)
	} else{ 
	    bareat[K,j] = mean_inter + barea[K, j] * sd_inter; // area-level random intercept
	}
  } 

  for (g in 1:ng)  {
    for (j in 1:narea) { 
      for (k in 1:K) {
	  if (ng > 1) {
	    beta[k,j,g] = bareat[k,j] + bmale[k]*(g-1); // additive gender and area effects
	  } else {
	    beta[k,j,g] = bareat[k,j];
	  }	    
      }

      // Infer age zero prevalence from data if there are any data at age zero
      if (prev_denom[1,j,g]==0 || prev_num[1,j,g]==0)
	prev[1,j,g] = 0;
      else prev[1,j,g] = prevzero[j,g];
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
	  cf[a,j,g] = exp(lcfbase[j]);
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
	
      if (remission) rem[,j,g] = rem_par[,j,g]; else { 	for (a in 1:nage) { rem[a,j,g] = 0; } }
      for (a in 1:nage){
	rem_prob[a,j,g] = 1 - exp(-rem[a,j,g]);
	P = trans_probs(inc[a,j,g], cf[a,j,g], rem[a,j,g]);
	inc_prob[a,j,g] = P[1,2] + P[1,3];
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
  mean_slope ~ normal(mism, miss); 

    // These all get transformed in different ways according to the model 
    for (j in 1:narea){
	for (i in 1:K){ 
	    barea[i,j] ~ normal(0, 1);
	}     
    }
  
  // Model for the data
  for (j in 1:narea){
    for (g in 1:ng) {
      mort_num[,j,g] ~ binomial(mort_denom[,j,g], mort[,j,g]);
      inc_num[,j,g] ~ binomial(inc_denom[,j,g], inc_prob[,j,g]);
      prev_num[,j,g] ~ binomial(prev_denom[,j,g], prev[,j,g]);
      if (remission) {
	rem_num[,j,g] ~ binomial(rem_denom[,j,g], rem_prob[,j,g]);
      }
    }
    if (increasing){
      lcfbase[j] ~ normal(mean_inter, sd_inter);  
    }
  }
  
    // Degree of smoothness
    lambda_smooth ~ exponential(sprior);
    if ((!interceptonly) && (!increasing)){
	// Variation in slopes 
	sd_slope ~ gamma(gpslope_a, gpslope_b);
    }
    // Variation in intercepts
    sd_inter ~ gamma(gpint_a, gpint_b);

  if (ng > 1) {
    lambda_smooth_male ~ exponential(sprior);
    // effect of being male, assumed common between areas
    //... on spline coefficients governing deviation from linearity
    for (i in 1:(K-2)){
      bmale[i] ~ normal(0, lambda_smooth_male); // smoothness variance par shared 
    }
    // ...on slopes 
    bmale[K-1] ~ normal(0, gender_slope_priorsd);
    // ...on intercepts
    bmale[K] ~ normal(0, gender_int_priorsd); 
    // five fold difference between men and women is log(5)=1.6 on log scale 
    // want sd such that 95% prob that bmale is between -log(5) and log(5)
    // sd = log(5) / qnorm(0.975) = 0.82 
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
