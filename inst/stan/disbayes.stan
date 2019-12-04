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

    // only in smoothed model 
    int<lower=0> K; // number of spline basis variables including the intercept
//    vector[K] beta_fix; 
    matrix[nage,K] X;
    /* matrix[K-1,2*(K-1)] S1; */
    real<lower=0> sprior; 

}

parameters {
    vector<lower=0>[nage] inc;
    vector[K] beta;
    vector<lower=0>[2] lambda;
}

transformed parameters {
    //// range constraints on these cause problems due to floating point fuzz
    vector<lower=0>[nage] cf;
    vector<lower=0>[nage] inc_prob;
//    vector[K] beta;
    matrix[nage+1,3] state_probs; 
    row_vector[3] tmp;
    matrix[3,3] P;
    real prev[nage];
    real mort[nage];
    /* matrix[K-1,K-1] K1; */
    /* vector[K] zero; */
    /* for (i in 1:K) { */
    /* 	zero[i] = 0; */
    /* } */
    /* K1 = S1[1:(K-1),1:(K-1)] * lambda[1]  + S1[1:(K-1),K:(2*(K-1))] * lambda[2]; */
    
    /// Case fatality as smooth spline function of age
    /// Spline basis X passed from R
    cf = exp(X*beta);
    
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

    // use jagam defaults for now 
    beta[1] ~ normal(0, 100); 
//	beta[2:10] ~ multi_normal(zero[2:K],K1);

    // diagonalised version
    /* beta ~ normal(beta_fix, 100); // lambda[1]); */
    for (i in 2:(K-1)) {
	beta[i] ~ normal(0, lambda[1]);
    }
    for (i in K:K) {
        beta[i] ~ normal(0, lambda[2]);
    }
    for (i in 1:2) {
        lambda[i] ~ exponential(sprior);
    }
////
}
