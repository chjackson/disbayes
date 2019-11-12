// TODO separate files for smooth/unsmooth for now. 
// Until get smooth model stable 

#include /include/trans_probs.stan

data {
    int smooth;
    int<lower=0> nage;
    int<lower=0> mort_num[nage];
    int<lower=0> mort_denom[nage];
    int<lower=0> prev_num[nage];
    int<lower=0> prev_denom[nage];
    int<lower=0> inc_num[nage];
    int<lower=0> inc_denom[nage];
    vector<lower=0>[nage] rem;
}

parameters {
    vector<lower=0>[nage] inc;
    vector<lower=0>[nage] cf;
}

transformed parameters {
    //// range constraints on these cause problems due to floating point fuzz
    vector<lower=0>[nage] inc_prob;
    matrix[nage+1,3] state_probs; 
    row_vector[3] tmp;
    matrix[3,3] P;
    real prev[nage];
    real mort[nage];
    
    state_probs[1,1] = 1;
    state_probs[1,2] = 0;
    state_probs[1,3] = 0;
    
    for (a in 1:nage){
	P = trans_probs(inc[a], cf[a], rem[a]);
	inc_prob[a] = P[1,2];
	prev[a] = state_probs[a,2] / (state_probs[a,1] + state_probs[a,2]);
	tmp = state_probs[a,1:3] * P;  // temp variable to avoid warning
	state_probs[a+1,1:3] = tmp;
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
    for (a in 1:nage){
	cf[a] ~ exponential(1);
    }
}
