#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc(double mu, double sigma, double lower, double upper);
double mn_approx(arma::rowvec theta, 
                 arma::colvec time, 
                 double sigma2, 
                 double h0, double h1, int c_seg);
  
// [[Rcpp::export]]
// main function
List Poi_fast(NumericMatrix count_mat, 
             NumericVector t_vec,
             NumericVector gamma_vec,
             NumericVector population,
             NumericVector h_prior,
             bool update_gamma, 
             double sig_error = 0.01,
             double omega = 0.01,
             int Niter = 20000){
  

  // disable printing of error messages
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  
  if(!update_gamma){
    Rcout<< "Run MCMC without updating change point location \n";
  }
  // model parameters //
  int p = count_mat.nrow();
  int V = count_mat.ncol();
  
  // other parameters //
  double theta_lim = 20;
  double theta_old;
  double theta_prop;
  int seg_size = 0; 
  double h0 = h_prior[0] / sig_error;
  double h1 = h_prior[1] / sig_error;
  
  int v_ind = 0;
  int count = 0;
  double log_mh_theta = 0.0;
  double ll_r = 0.0;
  double lpri_r = 0.0; 
  int acc_theta = 0;
  double val_temp = 0.0;
  double val_temp2 = 0.0;
  
  int idx_old = 0;
  int idx_new = 0;
  int v_cand = 1;

  int seg_size1 = 0;
  int seg_size2 = 0;
  arma::colvec c_temp(V, arma::fill::zeros);
  
  int idx_old_1 = 0;
  int idx_old_2 = 0;
  int idx_new_1 = 0;
  int idx_new_2 = 0;
  
  arma::rowvec theta_ori;
  arma::rowvec theta_ori2;
  arma::rowvec theta_star;
  arma::rowvec theta_star2;
  
  int K_star = 0; 
  double ll_r_temp;
  double log_mh_gamma = 0.0;
  int acc_gamma_add = 0;
  int acc_gamma_del = 0;
  int acc_gamma_swap = 0;
  
  int count_swap = 10; 
  int swap_idx = 0;
  int swap_star = 0;
  int swap_dir = 0; //swap_dir = 1 means swap with the left, and 2 means swap with the right
  int seg_size3 = 0; 

  
  // MAP //
  NumericVector posterior_nm(Niter / 2);
  arma::rowvec theta_post(V, arma::fill::zeros);
  
  int K_num = 0;
  double lg_theta = 0.0;
  arma::rowvec theta_seg;
  
  // PPM //
  arma::mat ppm(V, V, arma::fill::zeros);
  int ppm_cumu = 0;
  arma::colvec v_vec(Niter, arma::fill::zeros);
  
  // initialize the storage //
  arma::mat theta_store(Niter/2, V, arma::fill::zeros);
  arma::mat c_mat(V, Niter, arma::fill::zeros);
  arma::mat gamma_mat(Niter/2, V, arma::fill::zeros);
  arma::mat ll_store(Niter/2, V, arma::fill::zeros);
  
  
  // count matrix and theta matrix //
  arma::mat Y_mat(p, V, arma::fill::zeros);
  arma::mat theta_mat(p, V, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < V; j++){
      theta_mat(i, j) = log( (count_mat(i, j )+1) / population(i) );
      Y_mat(i, j) = count_mat(i, j );
    }
  }
  
  // c vector //
  int swap_flag = 0;
  arma::rowvec gamma_temp(V, arma::fill::zeros);
  arma::colvec c_vec(V, arma::fill::zeros);
  int acc = 0;
  for(int j = 0; j < V; j++){
    acc += gamma_vec[j];
    c_vec[j] = acc;
    gamma_temp[j] = gamma_vec[j];
  }
  
  // theta //
  arma::rowvec theta_vec_old;
  arma::rowvec theta_vec_prop;
  arma::rowvec theta_seg_old;
  arma::rowvec theta_seg_prop;
  arma::mat y_temp;
  arma::colvec t_temp;
  arma::colvec t_temp2;
  arma::colvec t_star;
  arma::colvec t_star2;
  arma::rowvec y_vec;
  
  // covariate information //
  arma::colvec time_vec; 
  time_vec = as<arma::colvec>(t_vec);

  for(int iter = 0; iter < Niter; iter ++){
    
    // monitor mcmc process //
    if(iter * 100/Niter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
   
    for(int i = 0; i < p; i ++ ){
      ////// 1.1 update alpha (exp (theta) ) in each segment //////
      // in the same state //
      y_vec = Y_mat.row(i);
      theta_vec_old = theta_mat.row(i);
      
      for(int j = 0; j < V; j++){
        // get the segment index of the current v //
        v_ind = c_vec[j];
        
        // extract the data in the same segment //
        y_temp = y_vec.cols(arma::find(c_vec == v_ind));
        
        // extract the covariate information for the segment //
        t_temp = time_vec.rows(arma::find(c_vec == v_ind));
        
        // update theta sequentially //
        theta_seg_old = theta_vec_old.cols(arma::find(c_vec == v_ind));
        seg_size = theta_seg_old.size();
        
        // propose a new theta using random walk //
        theta_old = theta_vec_old(j);
        theta_prop = rnorm_trunc(theta_old, 1, -theta_lim, theta_lim);
        
        theta_vec_prop = theta_vec_old;
        theta_vec_prop(j) = theta_prop;
        theta_seg_prop = theta_vec_prop.cols(arma::find(c_vec == v_ind));
        
        // likelihood ratio //
        val_temp = log(population(i)) + theta_prop;
        val_temp2 = log(population(i)) + theta_old;
        ll_r = y_vec(j) *(theta_prop - theta_old) - exp(val_temp) + exp(val_temp2); 
        
        // prior ratio //
        if(seg_size == 0){
          theta_vec_old(j) = theta_old;
          continue;
          
        }else{
          lpri_r = -mn_approx(theta_seg_old, t_temp, sig_error, h0, h1, seg_size) + 
            mn_approx(theta_seg_prop, t_temp, sig_error, h0, h1, seg_size);
         
          // MH for updating theta(j, v) //
          log_mh_theta = ll_r + lpri_r;
          
          if(log_mh_theta > log(unif_rand())){
            theta_vec_old(j) = theta_prop;
            acc_theta = acc_theta + 1;
          }else{
            theta_vec_old(j) = theta_old;
          }
        }
        
        // end for the current state (taxon) //
      }
      theta_mat.row(i) = theta_vec_old;
      
      if(update_gamma){
        ////// 2.update gamma vector //////
        v_cand = (int)(Rf_runif(1, V));
        while(gamma_temp[v_cand - 1] == 1 | gamma_temp[v_cand + 1] == 1 ) {
          v_cand = (int)(Rf_runif(1, V));
        }
        gamma_temp[v_cand] = 1 - gamma_temp[v_cand];
        v_vec[iter] = v_cand + 1;
        if(gamma_temp[v_cand] == 1) {
          // case 1: add //
          //seg_idx = c_vec[v_cand];
          c_temp = c_vec;
          acc = 0;
          for(int v = 0; v < V; v++){
            acc += gamma_temp[v];
            c_temp[v] = acc;
          }
          idx_old = c_vec[v_cand];
          idx_new = c_temp[v_cand];
          
          t_temp = time_vec.rows(arma::find(c_vec == idx_old));
          t_star = time_vec.rows(arma::find(c_temp == idx_new));
          t_star2 = time_vec.rows(arma::find(c_temp == (idx_new -1)));
          
          theta_ori = theta_mat.row(i);
          theta_ori = theta_ori.cols(arma::find(c_vec == idx_old));
          theta_star = theta_mat.row(i);
          theta_star = theta_star.cols(arma::find(c_temp == idx_new));
          theta_star2 = theta_mat.row(i);
          theta_star2 = theta_star2.cols(arma::find(c_temp == (idx_new -1)));
          
          seg_size = theta_ori.size();
          seg_size1 = theta_star.size();
          seg_size2 = theta_star2.size();
          if(seg_size != (seg_size1 + seg_size2 )){
            throw std::range_error("error length in add");
          }else{
            ll_r = mn_approx(theta_star, t_star, sig_error, h0, h1, seg_size1) + 
              mn_approx(theta_star2, t_star2, sig_error, h0, h1, seg_size2) -
              mn_approx(theta_ori, t_temp, sig_error, h0, h1, seg_size);
            
            lpri_r = log(omega) - log(1.0 - omega);
            log_mh_gamma = ll_r + lpri_r ;
            if(log_mh_gamma > log(unif_rand())){
              c_vec = c_temp;
              acc_gamma_add = acc_gamma_add + 1;
              gamma_temp(v_cand) = 1;
            }else{
              c_vec = c_vec;
              gamma_temp(v_cand) = 0;
            }
          }
          
          
        }else{
          // case 2: delete //
          //seg_idx = c_vec[v_cand];
          c_temp = c_vec;
          acc = 0;
          for(int vv = 0; vv < V; vv++){
            acc += gamma_temp[vv];
            c_temp[vv] = acc;
          }
          idx_old = c_vec[v_cand];
          idx_new = c_temp[v_cand];
          
          t_temp = time_vec.rows(arma::find(c_vec == idx_old));
          t_temp2 = time_vec.rows(arma::find(c_vec == (idx_old - 1)));
          t_star = time_vec.rows(arma::find(c_temp == idx_new));
          
          theta_ori = theta_mat.row(i);
          theta_ori = theta_ori.cols(arma::find(c_vec == idx_old));
          theta_ori2 = theta_mat.row(i);
          theta_ori2 = theta_ori2.cols(arma::find(c_vec == (idx_old - 1)));
          theta_star = theta_mat.row(i);
          theta_star = theta_star.cols(arma::find(c_temp == idx_new));
          
          seg_size = theta_star.size();
          seg_size1 = theta_ori.size();
          seg_size2 = theta_ori2.size();
          
          if(seg_size != (seg_size1 + seg_size2 )){
            throw std::range_error("error length in delete");
          }else{
            ll_r = mn_approx(theta_star, t_star, sig_error, h0, h1, seg_size) - 
              mn_approx(theta_ori, t_temp, sig_error, h0, h1, seg_size1) -
              mn_approx(theta_ori2, t_temp2, sig_error, h0, h1, seg_size2);
            lpri_r = log(1.0 - omega) - log(omega) ;
            log_mh_gamma = ll_r + lpri_r ;
            if(log_mh_gamma > log(unif_rand())){
              c_vec = c_temp;
              acc_gamma_del = acc_gamma_del + 1;
              gamma_temp(v_cand) = 0;
              
            }else{
              c_vec = c_vec;
              gamma_temp(v_cand) = 1;
            }
          }
          
          
        }
        // case 3: swap //
        if(iter % Niter == count_swap && iter != (Niter - 1 )){
          // check if we are able to perform swap //
          //Rcout << "swap" << std::endl;
          acc = 0;
          for(int vv = 1; vv < V; vv++){
            acc += gamma_temp[vv];
          }
          if(acc != 0){
            // perform swap for every several iterations //
            c_temp = c_vec;
            K_star = arma::max(c_temp);
            acc = 0;
            arma::rowvec cp_location(K_star - 1, arma::fill::zeros );
            for(int vv = 1; vv < V; vv ++){
              if(gamma_temp[vv] == 1){
                cp_location[acc] = vv;
                acc = acc + 1; 
              } 
            }
            
            swap_idx = (int)(Rf_runif(0, K_star - 1));
            swap_idx = cp_location[swap_idx];
            
            swap_dir = (int)(Rf_runif(0, 2)) + 1;
            while(gamma_temp[swap_idx - 1] == 1 && gamma_temp[swap_idx + 1] == 1){
              // proposed a new swap location //
              swap_idx = (int)(Rf_runif(0, K_star - 1));
              swap_idx = cp_location[swap_idx];
            }
            
            if(swap_idx == V - 1){
              swap_dir = 1;
            }else if(swap_idx == 1){
              swap_dir = 2;
            }else{
              if(swap_dir == 1 && gamma_temp[swap_idx - 1] == 1){
                swap_dir = 2;
              }else if(swap_dir == 2 && gamma_temp[swap_idx + 1] == 1){
                swap_dir = 1;
              }
            }
            
            
            if(swap_dir == 1){
              swap_star = swap_idx - 1;
            }else{
              swap_star = swap_idx + 1;
            }
            
            gamma_temp[swap_idx] = 1 - gamma_temp[swap_idx];
            gamma_temp[swap_star] = 1 - gamma_temp[swap_star];
            
            acc = 0;
            for(int vv = 0; vv < V; vv++){
              acc += gamma_temp[vv];
              c_temp[vv] = acc;
            }
            //seg_idx = c_vec[swap_idx];
            
            idx_old_1 = c_temp[swap_idx];
            idx_old_2 = c_vec[swap_star];
            idx_new_1 = c_temp[swap_idx];
            idx_new_2 = c_vec[swap_star];
            
            t_temp = time_vec.rows(arma::find(c_vec == idx_old_1));
            t_temp2 = time_vec.rows(arma::find(c_vec == idx_old_2));
            t_star = time_vec.rows(arma::find(c_temp == idx_new_1));
            t_star2 = time_vec.rows(arma::find(c_temp == idx_new_2));
            
            theta_ori = theta_mat.row(i);
            theta_ori = theta_ori.cols(arma::find(c_vec == idx_old_1));
            theta_ori2 = theta_mat.row(i);
            theta_ori2 = theta_ori2.cols(arma::find(c_vec == idx_old_2));
            theta_star = theta_mat.row(i);
            theta_star = theta_star.cols(arma::find(c_temp == idx_new_1));
            theta_star2 = theta_mat.row(i);
            theta_star2 = theta_star2.cols(arma::find(c_temp == idx_new_2));
            
            seg_size = theta_star.size();
            seg_size3 = theta_star2.size();
            seg_size1 = theta_ori.size();
            seg_size2 = theta_ori2.size();
            
            // check if we can perform swap //
            if((seg_size + seg_size3) != (seg_size1 + seg_size2 )){
              throw std::range_error("error length in swap");
            }else{
              ll_r = mn_approx(theta_star, t_star, sig_error, h0, h1, seg_size) + 
                mn_approx(theta_star2, t_star2, sig_error, h0, h1, seg_size3)-
                mn_approx(theta_ori, t_temp, sig_error, h0, h1, seg_size1)-
                mn_approx(theta_ori2, t_temp2, sig_error, h0, h1, seg_size2);
              
              log_mh_gamma = ll_r;
              
              // check group size after the swap 
              swap_flag = 0;
              for(int vv = 0; vv < V - 1; vv ++){
                if(gamma_temp[vv + 1] == gamma_temp[vv] && gamma_temp[vv + 1] == 1 ){
                  swap_flag = 1;
                }
              }
              
              if(log_mh_gamma > log(unif_rand()) && swap_flag != 1){
                c_vec = c_temp;
                acc_gamma_swap = acc_gamma_swap + 1;
              }else{
                c_vec = c_vec;
                gamma_temp[swap_idx] = 1 - gamma_temp[swap_idx];
                gamma_temp[swap_star] = 1 - gamma_temp[swap_star];
              }
              // monitor the process //
              count_swap = count_swap + 10;
              
              
            }
            
          }else{
            count_swap = count_swap + 10;
          }
        }
        
        
        
      }
      
    }
    
    if(iter >= Niter / 2){
      gamma_mat.row(iter - (Niter / 2)) = gamma_temp;
    }
    
    // summarize the result for the current iteration //
    c_mat.col(iter) = c_vec;
    if(iter >= Niter / 2){
      for(int i = 0; i < p; i++){
        theta_store.row(iter - (Niter / 2)) = theta_mat.row(i);
        
        // for calculate the importance rate 
        if(!update_gamma){
          theta_post = theta_mat.row(i);
          for(int v = 0; v < V; v ++){
            val_temp = log(population(i)) + theta_post(v);
            ll_store(iter - (Niter / 2), v) = y_vec(v) * val_temp - exp(val_temp) - lgamma(y_vec(v) + 1);
          }
          
        }
      }
      
      
      // caculate the unnormalized log-posterior //
      ll_r = 0.0;
      lg_theta = 0.0;
      for(int i = 0; i < p; i ++){
        y_vec = Y_mat.row(i);
        theta_post = theta_mat.row(i);
        for(int v = 0; v < V; v ++){
          val_temp = log(population(i)) + theta_post(v);
          // likelihood //
          ll_r = ll_r + y_vec(v) * val_temp - exp(val_temp);
          // prior for gamma // 
          ll_r = ll_r + gamma_temp(v) * log(omega) + (1 - gamma_temp(v)) * log(1 - omega);
        }
        // remove the prior for gamma[0]
        ll_r = ll_r - log(omega);
        K_num = arma::max(c_vec);
        // prior for theta (i.e. log (alpha)) // 
        for(int s = 1; s < K_num; s++){
          theta_seg = theta_post.cols(arma::find(c_vec == s));
          seg_size = theta_seg.size();
          t_temp = time_vec.rows(arma::find(c_vec == s));

          ll_r_temp = mn_approx(theta_seg, t_temp, sig_error, h0, h1, seg_size);
          lg_theta = lg_theta + ll_r_temp;
        }
      }
      
      posterior_nm(iter - (Niter / 2)) = lg_theta + ll_r;
    }
    
    // PPM //
    if(update_gamma){
      if(iter >= Niter / 2){
        for(int v = 0; v < V - 1; v ++){
          for(int vv = v + 1; vv < V; vv ++){
            if(c_vec(v) == c_vec(vv)){
              ppm(v, vv) = ppm(v, vv) + 1;
              ppm(vv, v) = ppm(vv, v) + 1;
            }
          }
        }
        ppm_cumu = ppm_cumu + 1;
      } 
    }
    // end mcmc //
  }
  
  if(update_gamma){
    // final PPM matrix //
    ppm = ppm / ppm_cumu;
    for(int v = 0; v < V; v ++)
    {
      ppm(v, v) = 1;
    }
  }
  
  Rcout<<count<< "MCMC done \n";
  
  // return //
  List result;
  if(update_gamma){
    result["gamma_mat"] = gamma_mat;
    result["z_mat"] = c_mat;
    result["ppm"] = ppm;
    result["posterior_trace"] = posterior_nm;
    result["theta_mcmc"] = theta_store;
    result["theta_accept_count"] = acc_theta ;
    result["gamma_add_count"] = acc_gamma_add ;
    result["gamma_delete_count"] = acc_gamma_del ;
    result["gamma_swap_count"] = acc_gamma_swap;
    result["update_idx"] = v_vec;
  }else{
    //result["z_mat"] = c_mat;
    result["ll_mcmc"] = ll_store;
    result["posterior_trace"] = posterior_nm;
    result["theta_mcmc"] = theta_store;
    result["theta_accept_count"] = acc_theta ;
    result["update_idx"] = v_vec;

  }
  return(result);
  
}


// utility functions //
// [[Rcpp::export]]
double mn_approx(arma::rowvec theta, 
                 arma::colvec time, 
                 double sigma2, 
                 double h0, double h1, int c_seg)
{
  
  double Stt = 0;
  double Stheta2 = 0;
  double Sthetat = 0;
  double theta_mean = mean(theta);
  double t_mean = mean(time);
  double p_dens = 0;
  
  
  for(int i = 0; i < c_seg; i ++){
    Stheta2 = Stheta2 + (theta[i] - theta_mean) * (theta[i] - theta_mean);
    Stt = Stt + (time[i] - t_mean) * (time[i] - t_mean);
    Sthetat = Sthetat + (theta[i] - theta_mean) * (time[i] - t_mean);
  }
  
  p_dens = -(c_seg / 2.0) * log(sigma2) - 0.5 * log(c_seg * h0 + 1.0) - 0.5 * log(h1 * Stt + 1.0) -
    (1.0 / (2.0 *sigma2) )* (Stheta2 -Sthetat*Sthetat / (Stt + 1.0 / h1) );
  
  return (p_dens);
}



double rnorm_trunc(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf)||(b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  
  // Second scenario
  else if((a*b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2)
    {
      z = unif_rs(a,b);
    }
    else if((lograt > logt1)&&(a < t3))
    {
      z = half_norm_rs(a,b);
    }
    else
    {
      z = exp_rs(a,b);
    }
    if(change)
    {
      z = -z;
    }
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

double exp_rs(double a, double b)
{
  double  z, u, rate;
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a))
  {
    z = R::rexp(rate);
  }
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a))
    {
      z = R::rexp(rate);
    }
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0)
  {
    xstar = 0.0;
  }
  else
  {
    xstar = a;
  }
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

double half_norm_rs(double a, double b)
{
  double x;
  x = fabs(norm_rand());
  while((x<a)||(x>b))
  {
    x = fabs(norm_rand());
  }
  return x;
}

double norm_rs(double a, double b)
{
  double x;
  x = Rf_rnorm(0.0, 1.0);
  while((x < a)||(x > b))
  {
    x = norm_rand();
  }
  return x;
}


