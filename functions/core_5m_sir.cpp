#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;

static IntegerVector sort(IntegerVector v);
static double mean_double_cpp(NumericVector x);
static double mean_cpp(IntegerVector x);
static int max_cpp(int a, int b);
static int min_cpp(int a, int b);

// [[Rcpp::export]]
Rcpp::List sir_c(IntegerVector C, int N, double gamma0, IntegerVector z, int T_fin, bool store) {
  // Read data information
  int T = C.length();
  int K = max(z) + 1;
  bool gamma_unknown = (gamma0 == -1.0);
  IntegerVector S(T), DeltaC(T, 0), flag(T, 1);
  int t;
  DeltaC(0) = C(0);
  for (t = 0; t < T; t++)
  {
    S(t) = N - C(t);
    if (t > 0)
    {
      DeltaC(t) = C(t) - C(t - 1);
      if (DeltaC(t) <= 0)
      {
        flag(t) = 0;
      }
    }
  }
  flag(0) = 0;
  flag(1) = 0;
  
  // Set algorithm settings
  int iter = 100000;
  int burn = iter*0.5; 
  int error = 1; // 1 for NB variance and 0 for Poisson variance
  double phi_start = 10, gamma_start, R0_start = 1;
  if (error == 0)
  {
    phi_start = -1.0;
  }
  if (gamma_unknown)
  {
    gamma_start = 0.1;
  }
  else
  {
    gamma_start = gamma0; 
  }
  double tau_logphi_ds = 1.0, tau_logphi_dr = 1.0, tau_logbeta = 0.1, tau_loggamma = 0.1, tau_logR0 = 0.1;
  
  // Set hyperparameters
  double a_phi = 0.001, b_phi = 0.001; 
  double a_R0 = 1, b_R0 = 1;
  double a_gamma = 0.01, b_gamma = 0.2;
  
  // Set temporary variables
  int it, k, count = 0, count_2 = 0;
  double phi_temp, beta_temp, beta_temp_2, gamma_temp, gamma_temp_2, R0_temp;
  double hastings;
  double accept_phi_ds = 0, accept_phi_dr = 0, accept_R0 = 0, accept_gamma = 0;
  NumericMatrix phi_ds_store(iter - burn, K);
  NumericMatrix phi_dr_store(iter - burn, K);
  NumericMatrix gamma_store(iter - burn, K);
  NumericMatrix R0_store(iter - burn, K);
  IntegerVector I(T, 0);
  IntegerVector R(T, 0);
  IntegerMatrix I_store(iter - burn, T);
  IntegerMatrix R_store(iter - burn, T);
  IntegerMatrix S_predict(200, T_fin);
  IntegerMatrix I_predict(200, T_fin);
  IntegerMatrix R_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  IntegerMatrix C_predict(200, T_fin);
  NumericMatrix S_fit(200, T + T_fin);
  NumericMatrix I_fit(200, T + T_fin);
  NumericMatrix R_fit(200, T + T_fin);
  NumericMatrix N_fit(200, T + T_fin);
  NumericMatrix C_fit(200, T + T_fin);
  
  // Initialization
  NumericVector phi_ds(K, phi_start);
  NumericVector phi_dr(K, phi_start);
  NumericVector gamma(K, gamma_start);
  NumericVector R0(K, R0_start);
  IntegerVector DeltaS(T - 1);
  IntegerVector DeltaR(T - 1);
  // Estimate I and R
  I(1) = C(1);
  for (t = 2; t < T; t++)
  {
    if (error == 1)
    {
      R(t) = min_cpp(R(t - 1) + rnbinom_mu(1, phi_dr(z(t)), I(t - 1)*gamma(z(t)))(0), C(t));
    }
    else
    {
      R(t) = min_cpp(R(t - 1) + rpois(1, I(t - 1)*gamma(z(t)))(0), C(t));
    }
    I(t) = I(t - 1) + DeltaC(t) - (R(t) - R(t - 1));
  }
  for (t = 0; t < T - 1; t++)
  {
    DeltaS(t) = S(t + 1) - S(t);
    DeltaR(t) = R(t + 1) - R(t);
  }
  
  // MCMC
  for (it = 0; it < iter; it++)
  {
    // Update phi_ds
    if (error == 1)
    {
      for (k = 0; k < K; k++)
      {
        phi_temp = exp(r_truncnorm(log(phi_ds(k)), tau_logphi_ds, log(1), log(100)));
        // beta_temp = beta(k);
        beta_temp = R0(k)*gamma(k);
        hastings = 0;
        for(t = 0; t < T - 1; t++)
        {
          if (flag(t) == 1 && z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
          {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp - DeltaS(t)) - (phi_temp - DeltaS(t))*log(phi_temp + beta_temp*S(t)*I(t)/N);
            hastings = hastings - (phi_ds(k)*log(phi_ds(k)) - lgamma(phi_ds(k)) + lgamma(phi_ds(k) - DeltaS(t)) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi_ds(k)) - b_phi*phi_ds(k));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi_ds(k) = phi_temp;
          if (it > burn) {
            accept_phi_ds++;
          }
        }
      }
    }
    
    // Update phi_dr
    if (error == 1)
    {
      for (k = 0; k < K; k++)
      {
        phi_temp = exp(r_truncnorm(log(phi_dr(k)), tau_logphi_dr, log(1), log(100)));
        gamma_temp = gamma(k);
        hastings = 0;
        for(t = 0; t < T - 1; t++)
        {
          if (flag(t) == 1 && z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
          {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + DeltaR(t)) - (phi_temp + DeltaR(t))*log(phi_temp + gamma_temp*I(t));
            hastings = hastings - (phi_dr(k)*log(phi_dr(k)) - lgamma(phi_dr(k)) + lgamma(phi_dr(k) + DeltaR(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp*I(t)));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi_dr(k)) - b_phi*phi_dr(k));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi_dr(k) = phi_temp;
          if (it > burn) {
            accept_phi_dr++;
          }
        }
      }
    }
    
    // Update R0
    for (k = 0; k < K; k++)
    {
      R0_temp = exp(r_truncnorm(log(R0(k)), tau_logR0, log(0.1), log(20)));
      beta_temp = R0(k)*gamma(k);
      beta_temp_2 = R0_temp*gamma(k);
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (flag(t) == 1 && z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
        {
          if (error == 0)
          {
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - beta_temp_2*S(t)*I(t)/N);
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - beta_temp*S(t)*I(t)/N);
          }
          else
          {
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp_2*S(t)*I(t)/N));
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
          }
        }
      }
      hastings = hastings + (a_R0 - 1)*log(R0_temp) - b_R0*R0_temp;
      hastings = hastings - ((a_R0 - 1)*log(R0(k)) - b_R0*R0(k));
      if(hastings >= log(double(rand()%10001)/10000))
      {
        R0(k) = R0_temp;
        if (it > burn) {
          accept_R0++;
        }
      }
    }
    
    // Update gamma
    if(gamma_unknown)
    {
      for (k = 0; k < K; k++)
      {
        gamma_temp = gamma(k);
        gamma_temp_2 = exp(r_truncnorm(log(gamma(k)), tau_loggamma, log(a_gamma), log(b_gamma)));
        beta_temp = R0(k)*gamma(k);
        beta_temp_2 = R0(k)*gamma_temp_2;
        hastings = 0;
        for(t = 0; t < T - 1; t++)
        {
          if (flag(t) == 1 && z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
          {
            if (error == 0)
            {
              hastings = hastings + (DeltaR(t)*log(gamma_temp_2*I(t)) - gamma_temp_2*I(t));
              hastings = hastings - (DeltaR(t)*log(gamma_temp*I(t)) - gamma_temp*I(t));
              hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - beta_temp_2*S(t)*I(t)/N);
              hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - beta_temp*S(t)*I(t)/N);
            }
            else
            {
              hastings = hastings + (DeltaR(t)*log(gamma_temp_2*I(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp_2*I(t)));
              hastings = hastings - (DeltaR(t)*log(gamma_temp*I(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp*I(t)));
              hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp_2*S(t)*I(t)/N));
              hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
            }
          }
        }
        // hastings = hastings + (a_gamma - 1)*log(gamma_temp_2) - (b_gamma - 1)*log(1 - gamma_temp_2);
        // hastings = hastings - ((a_gamma - 1)*log(gamma_temp) - (b_gamma - 1)*log(1 - gamma_temp));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          gamma(k) = gamma_temp_2;
          if (it > burn) {
            accept_gamma++;
          }
        }
      }
    }
    
    // Update R and I
    I(1) = C(1);
    for (t = 2; t < T; t++)
    {
      if (error == 1)
      {
        R(t) = min_cpp(R(t - 1) + rnbinom_mu(1, phi_dr(z(t)), I(t - 1)*gamma(z(t)))(0), C(t));
      }
      else
      {
        R(t) = min_cpp(R(t - 1) + rpois(1, I(t - 1)*gamma(z(t)))(0), C(t));
      }
      I(t) = I(t - 1) + DeltaC(t) - (R(t) - R(t - 1));
    }
    for (t = 0; t < T - 1; t++)
    {
      DeltaR(t) = R(t + 1) - R(t);
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if (it >= burn)
    {
      R0_store(it - burn, _) = R0;
      phi_ds_store(it - burn, _) = phi_ds;
      phi_dr_store(it - burn, _) = phi_dr;
      gamma_store(it - burn, _) = gamma;
      R_store(it - burn, _) = R;
      I_store(it - burn, _) = I;
      
      if (it % ((iter - burn)/200) == 0)
      {
        // Fitting
        for (t = 0; t < T + T_fin; t++)
        {
          if (t < T)
          {
            k = z(t);
          }
          else
          {
            k = z(T - 1);
          }
          if (t == 0)
          {
            S_fit(count_2, t) = S(t);
            I_fit(count_2, t) = I(t);
            R_fit(count_2, t) = R(t);
            C_fit(count_2, t) = C(t);
            N_fit(count_2, t) = 0;
          }
          else if (t <= T)
          {
            S_fit(count_2, t) = S(t - 1) - R0(k)*gamma(k)*S(t - 1)*I(t - 1)/N;
            R_fit(count_2, t) = R(t - 1) + gamma(k)*I(t - 1);
            I_fit(count_2, t) = N - S_fit(count_2, t) - R_fit(count_2, t);
            C_fit(count_2, t) = I_fit(count_2, t) + R_fit(count_2, t);
            N_fit(count_2, t) = C_fit(count_2, t) - C_fit(count_2, t - 1);
          }
          else
          {
            S_fit(count_2, t) = S_fit(count_2, t - 1) - R0(k)*gamma(k)*S_fit(count_2, t - 1)*I_fit(count_2, t - 1)/N;
            R_fit(count_2, t) = R_fit(count_2, t - 1) + gamma(k)*I_fit(count_2, t - 1);
            I_fit(count_2, t) = N - S_fit(count_2, t) - R_fit(count_2, t);
            C_fit(count_2, t) = I_fit(count_2, t) + R_fit(count_2, t);
            N_fit(count_2, t) = C_fit(count_2, t) - C_fit(count_2, t - 1);
          }
        }
        // Calibration
        for (t = 1; t < T + T_fin; t++)
        {
          if (C_fit(count_2, t) - C_fit(count_2, t - 1) <= 0)
          {
            N_fit(count_2, t) = N_fit(count_2, t - 1);
          }
        }
        
        // Prediction
        beta_temp = R0(K - 1)*gamma(K - 1);
        gamma_temp = gamma(K - 1);
        for (t = 0; t < T_fin; t++) {
          if (t == 0)
          {
            if (error == 1)
            {
              S_predict(count_2, t) = S(T - 1) - rnbinom_mu(1, phi_ds(K - 1), beta_temp*S(T - 1)*I(T - 1)/N)(0);
              R_predict(count_2, t) = R(T - 1) + rnbinom_mu(1, phi_dr(K - 1), gamma_temp*I(T - 1))(0);
            }
            else
            {
              S_predict(count_2, t) = S(T - 1) - rpois(1, beta_temp*S(T - 1)*I(T - 1)/N)(0);
              R_predict(count_2, t) = R(T - 1) + rpois(1, gamma_temp*I(T - 1))(0);
            }
          }
          else
          {
            if (error == 1)
            {
              S_predict(count_2, t) = S_predict(count_2, t - 1) - rnbinom_mu(1, phi_ds(K - 1), beta_temp*S_predict(count_2, t - 1)*I_predict(count_2, t - 1)/N)(0);
              R_predict(count_2, t) = R_predict(count_2, t - 1) + rnbinom_mu(1, phi_dr(K - 1), gamma_temp*I_predict(count_2, t - 1))(0);
            }
            else
            {
              S_predict(count_2, t) = S_predict(count_2, t - 1) - rpois(1, beta_temp*S_predict(count_2, t - 1)*I_predict(count_2, t - 1)/N)(0);
              R_predict(count_2, t) = R_predict(count_2, t - 1) + rpois(1, gamma_temp*I_predict(count_2, t - 1))(0);
            }
          }
          I_predict(count_2, t) = N - S_predict(count_2, t) - R_predict(count_2, t);
          if (I_predict(count_2, t) < 0)
          {
            I_predict(count_2, t) = 0;
          }
          C_predict(count_2, t) = I_predict(count_2, t) + R_predict(count_2, t);
          if (t == 0)
          {
            N_predict(count_2, t) = C_predict(count_2, t) - (I(T - 1) + R(T - 1));
          }
          else
          {
            N_predict(count_2, t) = C_predict(count_2, t) - C_predict(count_2, t - 1);
          }
        }
        count_2 = count_2 + 1;
      }
    }
  }
  
  // Prediction
  NumericVector S_fit_mean(T + T_fin);
  NumericVector I_fit_mean(T + T_fin);
  NumericVector R_fit_mean(T + T_fin);
  NumericVector C_fit_mean(T + T_fin);
  NumericVector N_fit_mean(T + T_fin);
  NumericVector S_mean(T_fin);
  IntegerVector S_upp(T_fin);
  IntegerVector S_lwr(T_fin);
  NumericVector I_mean(T_fin);
  IntegerVector I_upp(T_fin);
  IntegerVector I_lwr(T_fin);
  NumericVector R_mean(T_fin);
  IntegerVector R_upp(T_fin);
  IntegerVector R_lwr(T_fin);
  NumericVector C_mean(T_fin);
  IntegerVector C_upp(T_fin);
  IntegerVector C_lwr(T_fin);
  NumericVector N_mean(T_fin);
  IntegerVector N_upp(T_fin);
  IntegerVector N_lwr(T_fin);
  IntegerVector vtemp;
  for (t = 0; t < T + T_fin; t++) 
  {
    S_fit_mean(t) = mean_double_cpp(S_fit.column(t));
    I_fit_mean(t) = mean_double_cpp(I_fit.column(t));
    R_fit_mean(t) = mean_double_cpp(R_fit.column(t));
    C_fit_mean(t) = mean_double_cpp(C_fit.column(t));
    N_fit_mean(t) = mean_double_cpp(N_fit.column(t));
  }
  for (t = 0; t < T_fin; t++) 
  {
    vtemp = sort(S_predict.column(t));
    S_lwr(t) = vtemp(4);
    S_upp(t) = vtemp(194);
    S_mean(t) = mean_cpp(vtemp);
    vtemp = sort(I_predict.column(t));
    I_lwr(t) = vtemp(4);
    I_upp(t) = vtemp(194);
    I_mean(t) = mean_cpp(vtemp);
    vtemp = sort(R_predict.column(t));
    R_lwr(t) = vtemp(4);
    R_upp(t) = vtemp(194);
    R_mean(t) = mean_cpp(vtemp);
    vtemp = sort(C_predict.column(t));
    C_lwr(t) = vtemp(4);
    C_upp(t) = vtemp(194);
    C_mean(t) = mean_cpp(vtemp);
    vtemp = sort(N_predict.column(t));
    N_lwr(t) = vtemp(4);
    N_upp(t) = vtemp(194);
    N_mean(t) = mean_cpp(vtemp);
  }
  
  // Result wrap-up
  NumericVector accept(4);
  accept(0) = accept_R0/(iter - burn)/K;
  accept(1) = accept_gamma/(iter - burn)/K;
  accept(2) = accept_phi_ds/(iter - burn)/K;
  accept(3) = accept_phi_dr/(iter - burn)/K;
  
  if (store)
  {
    return Rcpp::List::create(Rcpp::Named("S_pred") = S_predict, 
                              Rcpp::Named("I_pred") = I_predict, 
                              Rcpp::Named("R_pred") = R_predict,
                              Rcpp::Named("C_pred") = C_predict, 
                              Rcpp::Named("N_pred") = N_predict,
                              Rcpp::Named("S_fit") = S_fit, 
                              Rcpp::Named("I_fit") = I_fit, 
                              Rcpp::Named("R_fit") = R_fit,
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit,
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("accept") = accept,
                              Rcpp::Named("R0_store") = R0_store,
                              Rcpp::Named("phi_ds_store") = phi_ds_store, 
                              Rcpp::Named("phi_dr_store") = phi_dr_store, 
                              Rcpp::Named("gamma_store") = gamma_store,
                              Rcpp::Named("R_store") = R_store, 
                              Rcpp::Named("I_store") = I_store);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("C_pred_mean") = C_mean, 
                              Rcpp::Named("C_pred_upp") = C_upp, 
                              Rcpp::Named("C_pred_lwr") = C_lwr, 
                              Rcpp::Named("N_pred_mean") = N_mean, 
                              Rcpp::Named("N_pred_upp") = N_upp, 
                              Rcpp::Named("N_pred_lwr") = N_lwr,
                              Rcpp::Named("S_fit") = S_fit_mean, 
                              Rcpp::Named("I_fit") = I_fit_mean, 
                              Rcpp::Named("R_fit") = R_fit_mean,
                              Rcpp::Named("C_fit") = C_fit_mean, 
                              Rcpp::Named("N_fit") = N_fit_mean, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("R0_store") = R0_store,
                              Rcpp::Named("gamma_store") = gamma_store,
                              Rcpp::Named("R_store") = R_store, 
                              Rcpp::Named("I_store") = I_store);
  }
}

// [[Rcpp::export]]
Rcpp::List sir(IntegerVector I, IntegerVector R, int N, IntegerVector z, int T_fin, bool store) {
  // Read data information
  int T = I.length();
  IntegerVector S(T);
  int t;
  for (t = 0; t < T; t++)
  {
    S(t) = N - I(t) - R(t);
  }
  int K = max(z) + 1;
  
  // Set algorithm settings
  int iter = 100000;
  int burn = iter*0.5; 
  int error = 1; // 1 for NB variance and 0 for Poisson variance
  double phi_start = 10, gamma_start = 0.01, R0_start = 1;
  double tau_logphi_ds = 1.0, tau_logphi_dr = 1.0, tau_logbeta = 0.1, tau_loggamma = 0.1, tau_logR0 = 0.1;
  
  // Set hyperparameters
  double a_phi = 0.001, b_phi = 0.001; 
  double a_R0 = 1, b_R0 = 1;
  double a_gamma = 0.2, b_gamma = 1.8;
  
  // Set temporary variables
  int it, k, count = 0, count_2 = 0;;
  double phi_temp, beta_temp, beta_temp_2, gamma_temp, gamma_temp_2, R0_temp;
  double hastings, logposterior = 0, logposterior_map;
  double accept_phi_ds = 0, accept_phi_dr = 0, accept_R0 = 0, accept_gamma = 0;
  NumericVector phi_ds_map(K, phi_start);
  NumericVector phi_dr_map(K, phi_start);
  NumericVector gamma_map(K, gamma_start);
  NumericVector R0_map(K, R0_start);
  NumericVector logposterior_store(iter - burn);
  NumericMatrix phi_ds_store(iter - burn, K);
  NumericMatrix phi_dr_store(iter - burn, K);
  NumericMatrix gamma_store(iter - burn, K);
  NumericMatrix R0_store(iter - burn, K);
  IntegerMatrix S_predict(200, T_fin);
  IntegerMatrix I_predict(200, T_fin);
  IntegerMatrix R_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  IntegerMatrix C_predict(200, T_fin);
  
  // Initialization
  NumericVector phi_ds(K, phi_start);
  NumericVector phi_dr(K, phi_start);
  NumericVector gamma(K, gamma_start);
  NumericVector R0(K, R0_start);
  IntegerVector DeltaS(T - 1);
  IntegerVector DeltaR(T - 1);
  for (t = 0; t < T - 1; t++)
  {
    DeltaS(t) = S(t + 1) - S(t);
    DeltaR(t) = R(t + 1) - R(t);
    if (I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
    {
      k = z(t);
      beta_temp = gamma(k)*R0(k);
      gamma_temp = gamma(k);
      if (error == 0)
      {
        logposterior = logposterior + (-DeltaS(t))*log(beta_temp*S(t)*I(t)/N) - beta_temp*S(t)*I(t)/N;
        logposterior = logposterior + DeltaR(t)*log(gamma_temp*I(t)) - gamma_temp*I(t);
      }
      else
      {
        logposterior = logposterior + lgamma((-DeltaS(t)) + phi_ds(k)) - lgamma(phi_ds(k)) + phi_ds(k)*(log(phi_ds(k)) - log(beta_temp*S(t)*I(t)/N + phi_ds(k))) + (-DeltaS(t))*(log(beta_temp*S(t)*I(t)/N) - log(beta_temp*S(t)*I(t)/N + phi_ds(k)));
        logposterior = logposterior + lgamma(DeltaR(t) + phi_dr(k)) - lgamma(phi_dr(k)) + phi_dr(k)*(log(phi_dr(k)) - log(gamma_temp*I(t) + phi_dr(k))) + DeltaR(t)*(log(gamma_temp*I(t)) - log(gamma_temp*I(t) + phi_dr(k)));
      }
    }
  }
  for (k = 0; k < K; k++)
  {
    if (error == 1) {
      logposterior = logposterior + (a_phi - 1)*log(phi_ds(k)) - b_phi*phi_ds(k);
      logposterior = logposterior + (a_phi - 1)*log(phi_dr(k)) - b_phi*phi_dr(k);
    }
    else
    {
      if (k == 0)
      {
        phi_ds = -1;
        phi_dr = -1;
      }
    }
    logposterior = logposterior + (a_gamma - 1)*log(gamma(k)) - (b_gamma - 1)*log(1 - gamma(k));
    logposterior = logposterior + (a_R0 - 1)*log(R0(k)) - b_R0*R0(k);
  }
  
  // MCMC
  for (it = 0; it < iter; it++)
  {
    // Update phi_ds
    if (error == 1)
    {
      for (k = 0; k < K; k++)
      {
        phi_temp = exp(r_truncnorm(log(phi_ds(k)), tau_logphi_ds, log(1), log(100)));
        // beta_temp = beta(k);
        beta_temp = R0(k)*gamma(k);
        hastings = 0;
        for(t = 0; t < T - 1; t++)
        {
          if (z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
          {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp - DeltaS(t)) - (phi_temp - DeltaS(t))*log(phi_temp + beta_temp*S(t)*I(t)/N);
            hastings = hastings - (phi_ds(k)*log(phi_ds(k)) - lgamma(phi_ds(k)) + lgamma(phi_ds(k) - DeltaS(t)) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi_ds(k)) - b_phi*phi_ds(k));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi_ds(k) = phi_temp;
          logposterior = logposterior + hastings;
          if (it > burn) {
            accept_phi_ds++;
          }
        }
      }
    }
    
    // Update phi_dr
    if (error == 1)
    {
      for (k = 0; k < K; k++)
      {
        phi_temp = exp(r_truncnorm(log(phi_dr(k)), tau_logphi_dr, log(1), log(100)));
        gamma_temp = gamma(k);
        hastings = 0;
        for(t = 0; t < T - 1; t++)
        {
          if (z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
          {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + DeltaR(t)) - (phi_temp + DeltaR(t))*log(phi_temp + gamma_temp*I(t));
            hastings = hastings - (phi_dr(k)*log(phi_dr(k)) - lgamma(phi_dr(k)) + lgamma(phi_dr(k) + DeltaR(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp*I(t)));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi_dr(k)) - b_phi*phi_dr(k));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi_dr(k) = phi_temp;
          logposterior = logposterior + hastings;
          if (it > burn) {
            accept_phi_dr++;
          }
        }
      }
    }
    
    // Update R0
    for (k = 0; k < K; k++)
    {
      R0_temp = exp(r_truncnorm(log(R0(k)), tau_logR0, log(0.1), log(20)));
      beta_temp = R0(k)*gamma(k);
      beta_temp_2 = R0_temp*gamma(k);
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
        {
          if (error == 0)
          {
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - beta_temp_2*S(t)*I(t)/N);
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - beta_temp*S(t)*I(t)/N);
          }
          else
          {
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp_2*S(t)*I(t)/N));
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
          }
        }
      }
      hastings = hastings + (a_R0 - 1)*log(R0_temp) - b_R0*R0_temp;
      hastings = hastings - ((a_R0 - 1)*log(R0(k)) - b_R0*R0(k));
      if(hastings >= log(double(rand()%10001)/10000))
      {
        R0(k) = R0_temp;
        //beta(k) = beta_temp_2;
        logposterior = logposterior + hastings;
        if (it > burn) {
          accept_R0++;
        }
      }
    }
    
    // Update gamma
    for (k = 0; k < K; k++)
    {
      gamma_temp = gamma(k);
      gamma_temp_2 = exp(r_truncnorm(log(gamma(k)), tau_loggamma, log(10e-9), log(1 - 1e-9)));
      beta_temp = R0(k)*gamma(k);
      beta_temp_2 = R0(k)*gamma_temp_2;
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (z(t) == k && I(t) != 0 && DeltaR(t) > 0 && DeltaS(t) < 0)
        {
          if (error == 0)
          {
            hastings = hastings + (DeltaR(t)*log(gamma_temp_2*I(t)) - gamma_temp_2*I(t));
            hastings = hastings - (DeltaR(t)*log(gamma_temp*I(t)) - gamma_temp*I(t));
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - beta_temp_2*S(t)*I(t)/N);
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - beta_temp*S(t)*I(t)/N);
          }
          else
          {
            hastings = hastings + (DeltaR(t)*log(gamma_temp_2*I(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp_2*I(t)));
            hastings = hastings - (DeltaR(t)*log(gamma_temp*I(t)) - (phi_dr(k) + DeltaR(t))*log(phi_dr(k) + gamma_temp*I(t)));
            hastings = hastings + (-DeltaS(t)*log(beta_temp_2*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp_2*S(t)*I(t)/N));
            hastings = hastings - (-DeltaS(t)*log(beta_temp*S(t)*I(t)/N) - (phi_ds(k) - DeltaS(t))*log(phi_ds(k) + beta_temp*S(t)*I(t)/N));
          }
        }
      }
      hastings = hastings + (a_gamma - 1)*log(gamma_temp_2) - (b_gamma - 1)*log(1 - gamma_temp_2);
      hastings = hastings - ((a_gamma - 1)*log(gamma_temp) - (b_gamma - 1)*log(1 - gamma_temp));
      if(hastings >= log(double(rand()%10001)/10000))
      {
        gamma(k) = gamma_temp_2;
        logposterior = logposterior + hastings;
        //Rcout<<k<<": "<<gamma_temp<<"->"<<gamma_temp_2<<": "<<hastings<<" logp:"<<"\n";
        if (it > burn) {
          accept_gamma++;
        }
      }
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if (it >= burn)
    {
      R0_store(it - burn, _) = R0;
      gamma_store(it - burn, _) = gamma;
      if (it == burn)
      {
        logposterior_map = logposterior;
        phi_ds_map = phi_ds;
        phi_dr_map = phi_dr;
        gamma_map = gamma;
        R0_map = R0;
      }
      else
      {
        if(logposterior > logposterior_map)
        {
          logposterior_map = logposterior;
          phi_ds_map = phi_ds;
          phi_dr_map = phi_dr;
          gamma_map = gamma;
          R0_map = R0;
        }
      }
      logposterior_store(it - burn) = logposterior;
      if (store)
      {
        phi_ds_store(it - burn, _) = phi_ds;
        phi_dr_store(it - burn, _) = phi_dr;
      }
      //Prediction
      beta_temp = R0(K - 1)*gamma(K - 1);
      gamma_temp = gamma(K - 1);
      if (it % ((iter - burn)/200) == 0)
      {
        for (t = 0; t < T_fin; t++) {
          if (t == 0)
          {
            if (error == 1)
            {
              S_predict(count_2, t) = S(T - 1) - rnbinom_mu(1, phi_ds(K - 1), beta_temp*S(T - 1)*I(T - 1)/N)(0);
              R_predict(count_2, t) = R(T - 1) + rnbinom_mu(1, phi_dr(K - 1), gamma_temp*I(T - 1))(0);
            }
            else
            {
              S_predict(count_2, t) = S(T - 1) - rpois(1, beta_temp*S(T - 1)*I(T - 1)/N)(0);
              R_predict(count_2, t) = R(T - 1) + rpois(1, gamma_temp*I(T - 1))(0);
            }
          }
          else
          {
            if (error == 1)
            {
              S_predict(count_2, t) = S_predict(count_2, t - 1) - rnbinom_mu(1, phi_ds(K - 1), beta_temp*S_predict(count_2, t - 1)*I_predict(count_2, t - 1)/N)(0);
              R_predict(count_2, t) = R_predict(count_2, t - 1) + rnbinom_mu(1, phi_dr(K - 1), gamma_temp*I_predict(count_2, t - 1))(0);
            }
            else
            {
              S_predict(count_2, t) = S_predict(count_2, t - 1) - rpois(1, beta_temp*S_predict(count_2, t - 1)*I_predict(count_2, t - 1)/N)(0);
              R_predict(count_2, t) = R_predict(count_2, t - 1) + rpois(1, gamma_temp*I_predict(count_2, t - 1))(0);
            }
          }
          I_predict(count_2, t) = N - S_predict(count_2, t) - R_predict(count_2, t);
          if (I_predict(count_2, t) < 0)
          {
            I_predict(count_2, t) = 0;
          }
          C_predict(count_2, t) = I_predict(count_2, t) + R_predict(count_2, t);
          if (t == 0)
          {
            N_predict(count_2, t) = C_predict(count_2, t) - (I(T - 1) + R(T - 1));
          }
          else
          {
            N_predict(count_2, t) = C_predict(count_2, t) - C_predict(count_2, t - 1);
          }
        }
        count_2 = count_2 + 1;
      }
    }
  }
  
  // Fit the MAP model
  NumericVector S_fit(T + T_fin, 0.0);
  NumericVector I_fit(T + T_fin, 0.0);
  NumericVector R_fit(T + T_fin, 0.0);
  NumericVector C_fit(T + T_fin, 0.0);
  NumericVector N_fit(T + T_fin, 0.0);
  for (t = 0; t < T + T_fin; t++)
  {
    if (t < T)
    {
      k = z(t);
    }
    else
    {
      k = z(T - 1);
    }
    if (t == 0)
    {
      S_fit(t) = S(t);
      I_fit(t) = I(t);
      R_fit(t) = R(t);
      C_fit(t) = I(t) + R(t);
      N_fit(t) = C_fit(t);
    }
    else if (t <= T)
    {
      S_fit(t) = S(t - 1) - R0_map(k)*gamma_map(k)*S(t - 1)*I(t - 1)/N;
      R_fit(t) = R(t - 1) + gamma_map(k)*I(t - 1);
      I_fit(t) = N - S_fit(t) - R_fit(t);
      C_fit(t) = I_fit(t) + R_fit(t);
      N_fit(t) = C_fit(t) - C_fit(t - 1);
    }
    else
    {
      S_fit(t) = S_fit(t - 1) - R0_map(k)*gamma_map(k)*S_fit(t - 1)*I_fit(t - 1)/N;
      R_fit(t) = R_fit(t - 1) + gamma_map(k)*I_fit(t - 1);
      I_fit(t) = N - S_fit(t) - R_fit(t);
      C_fit(t) = I_fit(t) + R_fit(t);
      N_fit(t) = C_fit(t) - C_fit(t - 1);
    }
  }
  // Calibration
  for (t = 1; t < T + T_fin; t++)
  {
    if (C_fit(t) - C_fit(t - 1) <= 0)
    {
      // C_fit(t) = C_fit(t - 1);
      N_fit(t) = N_fit(t - 1);
    }
  }

  // Prediction
  NumericVector S_mean(T_fin);
  IntegerVector S_upp(T_fin);
  IntegerVector S_lwr(T_fin);
  NumericVector I_mean(T_fin);
  IntegerVector I_upp(T_fin);
  IntegerVector I_lwr(T_fin);
  NumericVector R_mean(T_fin);
  IntegerVector R_upp(T_fin);
  IntegerVector R_lwr(T_fin);
  NumericVector C_mean(T_fin);
  IntegerVector C_upp(T_fin);
  IntegerVector C_lwr(T_fin);
  NumericVector N_mean(T_fin);
  IntegerVector N_upp(T_fin);
  IntegerVector N_lwr(T_fin);
  IntegerVector vtemp;
  for (t = 0; t < T_fin; t++) 
  {
    vtemp = sort(S_predict.column(t));
    S_lwr(t) = vtemp(4);
    S_upp(t) = vtemp(194);
    S_mean(t) = mean_cpp(vtemp);
    vtemp = sort(I_predict.column(t));
    I_lwr(t) = vtemp(4);
    I_upp(t) = vtemp(194);
    I_mean(t) = mean_cpp(vtemp);
    vtemp = sort(R_predict.column(t));
    R_lwr(t) = vtemp(4);
    R_upp(t) = vtemp(194);
    R_mean(t) = mean_cpp(vtemp);
    vtemp = sort(C_predict.column(t));
    C_lwr(t) = vtemp(4);
    C_upp(t) = vtemp(194);
    C_mean(t) = mean_cpp(vtemp);
    vtemp = sort(N_predict.column(t));
    N_lwr(t) = vtemp(4);
    N_upp(t) = vtemp(194);
    N_mean(t) = mean_cpp(vtemp);
  }
  
  // Result wrap-up
  NumericVector accept(4);
  accept(0) = accept_R0/(iter - burn)/K;
  accept(1) = accept_gamma/(iter - burn)/K;
  accept(2) = accept_phi_ds/(iter - burn)/K;
  accept(3) = accept_phi_dr/(iter - burn)/K;

  if (store)
  {
    return Rcpp::List::create(Rcpp::Named("I_fit") = I_fit, 
                              Rcpp::Named("R_fit") = R_fit,
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("S_pred") = S_predict, 
                              Rcpp::Named("I_pred") = I_predict, 
                              Rcpp::Named("R_pred") = R_predict, 
                              Rcpp::Named("C_pred") = C_predict, 
                              Rcpp::Named("N_pred") = N_predict, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("accept") = accept,
                              Rcpp::Named("R0_map") = R0_map,
                              Rcpp::Named("gamma_map") = gamma_map, 
                              Rcpp::Named("logposterior_store") = logposterior_store,
                              Rcpp::Named("R0_store") = R0_store,
                              Rcpp::Named("phi_ds_store") = phi_ds_store, 
                              Rcpp::Named("phi_dr_store") = phi_dr_store, 
                              Rcpp::Named("gamma_store") = gamma_store);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("S_pred_mean") = S_mean, 
                              Rcpp::Named("I_pred_mean") = I_mean, 
                              Rcpp::Named("R_pred_mean") = R_mean, 
                              Rcpp::Named("C_pred_mean") = C_mean, 
                              Rcpp::Named("C_pred_upp") = C_upp, 
                              Rcpp::Named("C_pred_lwr") = C_lwr, 
                              Rcpp::Named("N_pred_mean") = N_mean, 
                              Rcpp::Named("N_pred_upp") = N_upp, 
                              Rcpp::Named("N_pred_lwr") = N_lwr,
                              Rcpp::Named("S_fit") = S_fit, 
                              Rcpp::Named("I_fit") = I_fit, 
                              Rcpp::Named("R_fit") = R_fit,
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("gamma_map") = gamma_map,
                              //Rcpp::Named("logposterior_store") = logposterior_store,
                              Rcpp::Named("R0_map") = R0_map,
                              Rcpp::Named("R0_store") = R0_store,
                              Rcpp::Named("gamma_store") = gamma_store);
  }
}

IntegerVector sort(IntegerVector v) {
  std::sort(v.begin(), v.end());
  return v;
}

double mean_double_cpp(NumericVector x) {
  int n = x.size(); // Size of vector
  double sum = 0; // Sum value
  // For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    // Shorthand for sum = sum + x[i]
    sum = sum + x(i);
  }
  return sum/n; // Obtain and return the Mean
}

double mean_cpp(IntegerVector x) {
  int n = x.size(); // Size of vector
  double sum = 0; // Sum value
  // For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    // Shorthand for sum = sum + x[i]
    sum = sum + x(i);
  }
  return sum/n; // Obtain and return the Mean
}

int max_cpp(int a, int b) {
  if (b > a) 
  {
    return (b);
  } 
  else
  {
    return (a);
  }
}

int min_cpp(int a, int b) {
  if (b <= a) 
  {
    return (b);
  } 
  else
  {
    return (a);
  }
}