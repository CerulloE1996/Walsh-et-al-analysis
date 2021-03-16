// Direct comparisons only - code modified from Nyaga et al  https://journals.sagepub.com/doi/10.1177/0962280216669182

data {
    int N; // number of obs. / rows 
    int total_tests; // total number of tests 
    int NS; // total number of studies 
    int TP[N];
    int pos[N];
    int TN[N];
    int neg[N];
    int Study[N];
    int Test[N];
    int<lower=0,upper=1> holdout[N]; //index whether the observation should be used (for K-fold CV)
}
 
parameters {
    vector<lower=0>[2] tau;   
    vector[2] mu[total_tests];  
    cholesky_factor_corr[2] L_Omega; 
    vector<lower=0>[2] sigma; 
    matrix[NS, total_tests] z2[2];
    vector[2] z[NS];
}

transformed parameters {
    matrix[NS, 2] logit_pi[total_tests];
    matrix[NS, 2] log_lik[total_tests];
    matrix[NS, 2] log_lik_study;

    for (n in 1:N) {
        logit_pi[Test[n],Study[n],] = to_row_vector( mu[Test[n],] + 
                                       diag_pre_multiply(sigma, L_Omega) * to_vector(z2[,Study[n],Test[n]]) ); // for non-NMA

  // Pointwise (i.e. observation level) Likelihood 
    log_lik[Test[n], Study[n], 1] =  binomial_logit_lpmf(TP[n]  | pos[n], logit_pi[Test[n], Study[n], 1]);
    log_lik[Test[n], Study[n], 2] =  binomial_logit_lpmf(TN[n]  | neg[n], logit_pi[Test[n], Study[n], 2]);
    }
 for (s in 1:NS) {
    for (j in 1:2)
      log_lik_study[s,j] = sum(log_lik[ , s, j]); 
   }
}

model {
    //Prior Model
  for (j in 1:2) {
       mu[,j] ~ normal(0, 2);
       z[,j]  ~ std_normal();
    for (s in 1:NS)
       z2[j,s,] ~ std_normal();
      }

       tau ~ std_normal();
       sigma ~ std_normal();
       L_Omega ~ lkj_corr_cholesky(2);

  // Likelihood Model
    for (n in 1:N) {
      if(holdout[n] == 0) {
        target += log_lik[Test[n], Study[n], 1];
        target += log_lik[Test[n], Study[n], 2];
       }
    }
}

generated quantities { 
     vector[total_tests] Se; 
     vector[total_tests] Sp;  
     vector[total_tests] lSe; 
     vector[total_tests] lSp; 
     vector[choose(total_tests,2)] diff_Se;
     vector[choose(total_tests,2)] diff_Sp;
     vector[choose(total_tests,2)] ratio_Se;
     vector[choose(total_tests,2)] ratio_Sp;
     corr_matrix[2] Omega; 
     matrix[2,2] Sigma; 
     vector[total_tests] S; // Superiority index (Deutsch et al.) 
     matrix[total_tests, total_tests] A;
     matrix[total_tests, total_tests] B;
     matrix[total_tests, total_tests] C;
     vector<lower=0>[2] tausq;
     vector<lower=0>[2] sigmabsq;
     matrix[total_tests, total_tests] sigmasq[2];
     matrix[total_tests, total_tests] rho[2];
     real rho12;
     matrix[2,2] Sigma_bs; 
     vector[2] pred[total_tests]; 
     vector[total_tests] Se_pred; 
     vector[total_tests] Sp_pred;  
     vector[total_tests] lSe_pred; 
     vector[total_tests] lSp_pred; 

  Omega = multiply_lower_tri_self_transpose(L_Omega); 
  Sigma = quad_form_diag(Omega, sigma);

    sigmabsq[1] = Sigma[1,1];
    sigmabsq[2] = Sigma[2,2];

    for (j in 1:2) {
        tausq[j] = tau[j] * tau[j];
        for (k in 1:total_tests){
            for (l in 1:total_tests){
                sigmasq[j,k,l] = (sigmabsq[j] + tausq[j])*((sigmabsq[j] + tausq[j]));
                rho[j,k,l] = sigmabsq[j]/sqrt(sigmasq[j,k,l]);
            }
        }
    }
   
    for (j in 1:2) {
    //  for (l in 1:total_tests){
     //   for(m in 1:total_tests){

 rho12 =      Omega[1,1]*sqrt(Sigma[1,1])*sqrt(Sigma[2,2]) /
                   sqrt( (Sigma[1,1] + tausq[j]) * (Sigma[2,2] + tausq[j]) ); 
     }

    for (l in 1:total_tests){
        for(m in 1:total_tests){
            A[l, m] = if_else((mu[l,1] >  mu[m,2]) && ( mu[l,1] >  mu[m,2]), 1, 0);
            B[l, m] = if_else((mu[l,1] <  mu[m,2]) && ( mu[l,1] <  mu[m,2]), 1, 0);
            C[l, m] = if_else((mu[l,1] == mu[m,2]) && ( mu[l,1] == mu[m,2]), 1, 0);
        }

        S[l] = (2*sum(row(A, l)) + sum(row(C, l)))/(2*sum(row(B, l)) + sum(row(C, l)));
    }
   

  // construct the (homog) between-study var-cov. matrix, needed to generate prediction intervals 
   Sigma_bs[1,1] =  sigma[1]^2 + tau[1]^2; // diseased
   Sigma_bs[2,2] =  sigma[2]^2 + tau[2]^2; // non-diseased
   Sigma_bs[1,2] = rho12*sqrt(Sigma_bs[1,1])*sqrt(Sigma_bs[2,2]);
   Sigma_bs[2,1] = rho12*sqrt(Sigma_bs[1,1])*sqrt(Sigma_bs[2,2]);

 for (t in 1:total_tests) { 
  Se[t] = inv_logit(mu[t,1]); 
  Sp[t] = inv_logit(mu[t,2]); 
  lSe[t] = mu[t,1];
  lSp[t] = mu[t,2];
  pred[t, ] = multi_normal_rng( mu[t, ] , Sigma_bs );
  lSe_pred[t] = pred[t, 1];
  lSp_pred[t] = pred[t, 2];
  Se_pred[t] = inv_logit(pred[t, 1]);
  Sp_pred[t] = inv_logit(pred[t, 2]);
  }

  // calculate pairwise accuracy differences (5C2 = 10 for Se and 10 for Sp)
for (i in 1:2) {
 diff_Se[i] = Se[1] - Se[1 + i];
 diff_Sp[i] = Sp[1] - Sp[1 + i]; 
 ratio_Se[i] = Se[1] / Se[1 + i];
 ratio_Sp[i] = Sp[1] / Sp[1 + i]; 
}

 diff_Se[3] = Se[2] - Se[3];
 diff_Sp[3] = Sp[2] - Sp[3];
 ratio_Se[3] = Se[2] / Se[3];
 ratio_Sp[3] = Sp[2] / Sp[3];

}

