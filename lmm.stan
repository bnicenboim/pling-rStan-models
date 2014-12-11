## See Sorensen & Vasishth (under review) A tutorial on fitting Bayesian linear mixed models using Stan
# The code relies heavily on the code of Tanner Sorensen
# see: http://www.ling.uni-potsdam.de/~vasishth/statistics/BayesLMMs.html

#by saving in x, x_u, and x_w the output of model.matrix(formula), the model can be used for different lmms

data {
      int<lower=0> N_obs; 
      int<lower=0> N_coef;  //predictors +intercept
      int<lower=0> N_coef_u;  //predictors +intercept
      int<lower=0> N_coef_w;  //predictors +intercept

      int<lower=1> N_subj;                 //number of subjects
      int<lower=1> N_item;                 //number of items
      int<lower=1> subj[N_obs];    //subject id
      int<lower=1> item[N_obs];    //item id

      matrix[N_obs,N_coef] x;
      matrix[N_obs,N_coef_u] x_u;
      matrix[N_obs,N_coef_w] x_w;
      vector[N_obs] rt;

}

parameters {
    	vector[N_coef] beta;
    	real<lower=0,upper=10000> sigma;

      //subj
    	vector<lower=0> [N_coef_u]  sigma_u;     // subj sd
      cholesky_factor_corr[N_coef_u] L_u;      // correlation matrix for random intercepts and slopes subj
      matrix[N_coef_u,N_subj] z_u;

      //items
      vector<lower=0> [N_coef_w]  sigma_w;     // subj sd
      cholesky_factor_corr[N_coef_w] L_w;      // correlation matrix for random intercepts and slopes item
      matrix[N_coef_w,N_item] z_w;
}


transformed parameters {
     
      vector[N_obs] fixeff;
      vector[N_obs] raneff;
      vector[N_obs] mu;
      matrix[N_coef_u,N_subj]  u;         // random intercept and slopes subj
      matrix[N_coef_w,N_item]  w;         // random intercept and slopes subj

      fixeff <- x * beta;   
      u <- (diag_pre_multiply(sigma_u,L_u) * z_u); // subj random effects
      w <- (diag_pre_multiply(sigma_w,L_w) * z_w); // item random effects
    
        //  for (coef in 1:N_coef) 
        //     for (i in 1:N_obs) 
        //        raneff[i] <-  x_u[i,coef] * u[coef,subj[i]] +  x_w[i,coef] * w[coef,item[i]];     
      
      for (i in 1:N_obs)
        raneff[i] <-   row(x_u,i) * col(u,subj[i])+ row(x_w,i) * col(w,item[i]);

      mu <- fixeff+raneff;
} 

model {

      # priors:
      # beta ~ normal(0,10);  //the priors should be change for the data, or ommited
      # sigma ~ normal(0,10);
      # sigma_u ~ normal(0,10);
      # sigma_w ~ normal(0,10);
      L_u ~ lkj_corr_cholesky(4.0);
      L_w ~ lkj_corr_cholesky(4.0);
      # to_vector(z_u) ~ normal(0,1);
      # to_vector(z_w) ~ normal(0,1);
      rt ~ normal(mu, sigma);  // in the glmm with binomial link normal has to be replaced by binomial, and sigma has to be removed


}

generated quantities {
      matrix[N_coef_u,N_coef_u] Cor_u;
      matrix[N_coef_w,N_coef_w] Cor_w;

      Cor_u <- L_u * L_u';  //Correlations between random effects by subj
      Cor_w <- L_w * L_w';  //Correlations between random effects by item
}



