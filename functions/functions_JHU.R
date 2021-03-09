detect_changepoint = function(input.df, 
                              state_idx, 
                              h_vec,
                              NITER = 40000,
                              seed = seed
                              ){
  # extract state information #
  GRID_LEN = 10
  SIGMA2_GRID = exp(seq(log(0.0001), log(0.005), length.out = GRID_LEN))
 
  OMEGA = 0.05
  SELECT_NUM = 2
  
  all_states = unique(input.df$state)
  if(state_idx > 51 | state_idx < 0){
    stop("The state_idx should be in [1, 51]. \n ")
  }
  state_nam = all_states[state_idx]
  cat("Dectecting optimal change point for", state_nam, "\n")
  processed_res = input.df[input.df$state == state_nam, ]
  rownames(processed_res) = NULL
  population = unique(processed_res$population)
  
  # prepare model input #
  single_state = processed_res$I_obs
  
  Ymat = matrix(single_state, nrow = 1)
  t_vec = matrix(processed_res$time, ncol = 1)
  
  time_info = processed_res$date
  gamma_vector = rep(0, length(single_state))
  gamma_vector[1] = 1
  ini.cp = seq(1, length(single_state), by = 30)
  gamma_vector[ini.cp] = 1
  
  h_prior = h_vec
  cat("Starting",GRID_LEN, "round of model fitting... \n")
  psis_scores = rep(NA, GRID_LEN)
  ret_list = vector("list", GRID_LEN)
  
  for(rr in 1:GRID_LEN){
    #### fit step 1 ####
    cat("Step 1 in round:", rr ,"/", GRID_LEN, "\n")
    set.seed(seed)
    sigma2 = SIGMA2_GRID[rr]
    res_step1 = Poi_fast(Ymat, 
                         t_vec,
                         gamma_vec = gamma_vector,
                         population = population,
                         Niter = NITER,
                         omega = OMEGA,
                         h_prior = h_prior,
                         update_gamma = TRUE,
                         sig_error = sigma2)
    
    
    cat("Step 1 finished in round", rr ,"/", GRID_LEN, "\n")
    
    # select change points by PPM #
    ppm = res_step1$ppm
    hclust_ward = hclust(as.dist(1-ppm), method="ward.D")
    cls_ward = t(apply(matrix(1:20),1, function(k) cutree(hclust_ward,k=k)))
    
    cat("finding the optimal z vector \n")
    mbind_lg = minbinder(ppm,cls_ward, method="all",
                          include.lg=TRUE)
    
    z_best = mbind_lg$cl[1, ]
    gamma_ppm = rep(0, length(z_best))
    gamma_ppm[1] = 1
    for(ii in 2:length(z_best)){
      gamma_ppm[ii] = ifelse(z_best[ii] == z_best[ii - 1], 0, 1) 
    }
    gamma_ppm[length(gamma_ppm)] = 0
    for(ii in 2:length(gamma_ppm)){
      if(gamma_ppm[ii-1] == 1 & gamma_ppm[ii] == 1 ){
        gamma_ppm[ii] = 0
      }else{
        gamma_ppm[ii] = gamma_ppm[ii]
      }
    }
    
    if(sum(gamma_ppm) <= 1){
      res_list = list(psis_score = Inf,
                      gamma_ppm = gamma_ppm, 
                      population = population,
                      sigam2 = sigma2 
      )
      ret_list[[rr]] = res_list
      rm(res_step1)
      cat("Skip step 2 in round", rr ,"/", GRID_LEN, "\n")
      next
    }else{
      # find gamma with credible interval #
      gamma_mcmc = res_step1$gamma_mat
      gamma_ppi = apply(gamma_mcmc, 2, mean)
      cp_interval = get.gamma.ci(gamma.ppm = gamma_ppm, gamma.mcmc = gamma_mcmc)
      
      cp.window = lapply(cp_interval$details, detect.overlap )
      cp.pos = unlist(lapply(cp_interval$details, function(x){x$cp.idx}))
      names(cp.window) = cp.pos
      
      ii = 1
      cp.update = list()
      current.cps = cp.window[[ii]]
      while(ii <= (length(cp.pos) - 1)){
        #print(ii)
        next.cps = cp.window[[ii + 1]]
        check.overlap = unlist(lapply(current.cps, function(x){x %in% next.cps}))
        if(any(check.overlap)){
          current.cps = unique(c(current.cps, next.cps))
          cp.update[[ii]] = current.cps 
          ii = ii + 1
        }else{
          cp.update[[ii]] = current.cps 
          cp.update[[ii + 1]] = next.cps
          ii = ii + 1
          current.cps = cp.window[[ii]]
        }
      }
      
      dup.idx = duplicated(cp.update)
      if(any(dup.idx)){
        cp.update = cp.update[!dup.idx]
      }
      cp.update.list = list()
      for(ele in 1:length(cp.update)){
        current.ele = cp.update[[ele]]
        bool.idx = rep(NA, length(cp.window))
        for(ll in 1:length(cp.window) ){
          bool.idx[ll] = all(unlist(lapply(current.ele, function(x){x %in% cp.window[[ll]]})))
        }
        if(any(bool.idx)){
          idx = which(bool.idx)
          cp.update.list[[ele]]  = list(cp.idx = as.numeric(names(cp.window)[[idx]]),
                                        cp.ci = cp.window[[idx]])
        }else{
          cp.update.list[[ele]]  = list(cp.idx = median(cp.update[[ele]]),
                                        cp.ci = cp.update[[ele]])
        }
      }
      
      gamma.ppm.idx = unique(c(1, unlist(lapply(cp.update.list, function(x){x$cp.idx}))))
      gamma.ppm.update = rep(0, length(gamma_ppm))
      gamma.ppm.update[gamma.ppm.idx ] =1
      
      # gamma map #
      post_trace = res_step1$posterior_trace
      map_index = which.max(post_trace)
      gamma_map = gamma_mcmc[map_index, ]
      
      rm(res_step1)
      rm(gamma_mcmc)
      
      
      #### fit step 2 ####
      cat("Step 2 in round", rr ,"/", GRID_LEN, "\n")
      res_step2 = Poi_fast(Ymat, 
                           t_vec,
                           gamma_vec = gamma.ppm.update,
                           population = population,
                           Niter = NITER,
                           h_prior = h_prior,
                           update_gamma = FALSE,
                           sig_error = sigma2)
      
      cat("Step 2 finihsed in round", rr ,"/", GRID_LEN, "\n")
      
      # PSIS score #
      ll_mcmc = res_step2$ll_mcmc
      p_waic = sum(apply(ll_mcmc, 2, var))
      
      # fit the generalized Pareto distribution #
      psis_mcmc = apply(ll_mcmc, 2, get_refit_weight)
      ll_mcmc = apply(ll_mcmc, 2, function(x){x - max(x)})
      l_mat = exp(ll_mcmc)
      wt_mat = do.call(cbind, lapply(psis_mcmc, function(x){x$weights})) 
      idx_vec = seq(1,length(single_state))
      psis_single = lapply(idx_vec, cal_psis_loo, l_mat, wt_mat )
      psis_approx = -unlist(psis_single)
      psis_score = round(sum(psis_approx), 3)
      psis_scores[rr] = psis_score
      
      
      # estimate parameters in the regression model # 
      theta_mcmc = res_step2$theta_mcmc
      theta_res = t(apply(theta_mcmc, 2, 
                          function(x)c(mean(x), quantile(x, c(0.025, 0.975)))) )
      lga_df = as.data.frame(theta_res, stringsAsFactors = F)
      lga_df$DataQuality = processed_res$dataQualityGrade
      rownames(lga_df) = time_info
      colnames(lga_df)[1] = c("log_alpha")
      rm(res_step2)
      Xmat_input = cbind(1, t_vec)
      coeff_fit = lm.nb.new(gamma.v = gamma.ppm.update,
                            y = lga_df$log_alpha ,
                            h0 = h_prior[1],
                            h1 = h_prior[2],
                            sigma2 = sigma2,
                            X = Xmat_input)
      beta_est = lapply(coeff_fit$beta.ci,function(x){round(x, 3)})
      intercept_df = as.data.frame(do.call(rbind, lapply(beta_est, function(x){x[1, ]})))
      time_df = as.data.frame(do.call(rbind, lapply(beta_est, function(x){x[2, ]})))
      
      res_list = list(psis_score = psis_score,
                      gamma_ppm = gamma.ppm.update, 
                      gamma_map = gamma_map,
                      gamma_ppm_old = gamma_ppm,
                      gamma_ppi = gamma_ppi, 
                      change_points = gamma.ppm.idx,
                      change_point_details = cp.update.list,
                      loga_est = lga_df,
                      intercept = intercept_df,
                      time_effect = time_df,
                      population = population,
                      sigam2 = sigma2 
      )
      ret_list[[rr]] = res_list
    }
  } 
  names(ret_list) = names(psis_scores) = round(SIGMA2_GRID, 4)
  
  # select the top optimal settings #
  opt_score = base::sort(psis_scores)[1:SELECT_NUM]
  opt_settings = unlist(lapply(opt_score, function(x){which(psis_scores == x) } ))
  ret_opt = ret_list[opt_settings]
  rm(ret_list)
  
  # append the data quality and other information #
  ret_opt[[SELECT_NUM + 1]] = processed_res
  names(ret_opt)[SELECT_NUM + 1] = "data"
  
  # setting information #
  setting_list = list(sigma_grid = SIGMA2_GRID,
                      psis_socre = as.numeric(psis_scores),
                      beta_variance = h_prior, 
                      omega = OMEGA,
                      ret_number = SELECT_NUM)
  ret_opt[[SELECT_NUM + 2]] = setting_list
  names(ret_opt)[SELECT_NUM + 2] = "setting_info"
  
  return(ret_opt)
}

fit_SIR = function(ret_opt, T_pred = 3){
  dat = ret_opt$data
  recovery.rate = unique(dat$recovery_rate)
  region = as.character(unique(dat$state))
  
  pop = state.pop[state.pop$state == region,2 ]
  # prepare data for SIR model #
  ret_number = 1
  cp_list = ret_opt[1]
  sir_input_list = lapply(cp_list, get.sir.df, dat = dat)
  
  # fit SIR for each gamma vector #
  case_idx = 1
  sir_R0_segment = lapply(case_idx, 
                          fit.sir.segment, 
                          recovery.rate = recovery.rate,
                          sir_input_list = sir_input_list, 
                          pop = pop, 
                          T_pred = T_pred)
  
  names(sir_R0_segment) = names(ret_opt)[1:ret_number]
  
  return(sir_R0_segment )
  
}



#### other utility function ####
get_refit_weight = function(x){
  psis_obj = suppressWarnings(psis(x))
  ks = pareto_k_values(psis_obj)
  lw = weights(psis_obj, normalize = TRUE)[, 1]
  ret.list = list(k_val = ks, 
                  weights = exp(lw) )
  return(ret.list)
}

cal_psis_loo = function(idx, l_mat, wt_mat){
  el = l_mat[, idx]
  wt = wt_mat[, idx]
  log(sum(el * wt) / sum(wt))
}

get.gamma.ci = function(gamma.ppm, gamma.mcmc){
  gamma.idx = which(gamma.ppm == 1)[-1]
  wd.check = c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)
  gamma.ci = vector("list", length = length(gamma.idx))
  if(length(gamma.idx) != 0){
    for(ii in 1:length(gamma.idx) ){
      check.idx = gamma.idx[ii]
      gamma.check = gamma.mcmc[, check.idx]
      if(length(unique(gamma.check)) == 1){
        gamma.ci[[ii]] = list(idx = check.idx,
                              result = rep(99, 10))
        ii = ii + 1
      }else{
        res.t = rep(0, 10)
        count.idx = 1
        for(idx in wd.check){
          compare.idx = check.idx + idx
          #cat(idx, count.idx, "\n")
          if(compare.idx <= 1 | compare.idx > ncol(gamma.mcmc)){
            res.t[count.idx] = 99.0
            count.idx = count.idx + 1
          }else if(gamma.ppm[compare.idx] == 1){
            if(idx < 0){
              res.t[c(1:count.idx)] = 99.0
              count.idx = count.idx + 1
            }else{
              res.t[c(count.idx:10)] = 99.0
              break
            }
          }else{
            gamma.compare = gamma.mcmc[, compare.idx ]
            if(length(unique(gamma.compare)) == 1 ){
              res.t[count.idx] = 99.0
              count.idx = count.idx + 1
            }else{
              single.corr.test = cor.test(gamma.check, gamma.compare, "less")
              res.t[count.idx] = single.corr.test$p.value
              count.idx = count.idx + 1
            }
          }
        }
        gamma.ci[[ii]] = list(idx = check.idx, result = res.t)
      }
    }
    
  }else{
    # only select time = 1 as the change point #
    gamma.ci[[1]] = list(idx = 1, result = rep(99, length(wd.check)))
  }
  gamma.ci.list = lapply(gamma.ci,find.cp.window )
  
  cp.all = c()
  for(ii in 1:length(gamma.idx)){
    cp.idx = gamma.ci.list[[ii]][1]$cp.idx
    cp.window = gamma.ci.list[[ii]][2]$window.res
    wd.vec = wd.check[which(cp.window == 1)]
    wd.t.idx = cp.idx + wd.vec
    wd.t.idx = base::sort(c(wd.t.idx, cp.idx))
    cp.all = c(cp.all, wd.t.idx)
  }
  cp.all = c(1, unique(cp.all))
  cp.ppm = which(gamma.ppm == 1)
  cp.list = list(main = cp.ppm, 
                 invertal = cp.all, 
                 details = gamma.ci.list)
  return(cp.list)
}

find.cp.window = function(x, p.val = 0.05, adj = TRUE){
  wd.check = c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)
  cp.idx = x[[1]]
  full.window = x[[2]]
  if(min(full.window) == max(full.window)){
    window.res = rep(0, 10)
  }else{
    window.res = rep(0, 10)
    singular.idx =  which(full.window!=99)
    singular.idx2 = singular.idx+ 1
    rm.idx = unlist(sapply(singular.idx2, function(x){x %in% c(singular.idx, max(singular.idx) +1) }) )
    if(adj){
      keep.idx = which(full.window!=99)
      keep.idx2 = keep.idx[rm.idx]
      full.window2 = full.window[keep.idx]
      full.window2 = full.window2[rm.idx]
      p.adj = p.adjust( full.window2, method = "BH")
      sel.idx = keep.idx2[which(p.adj< p.val)]
      window.res[sel.idx] = 1
    }else{
      keep.idx = which(full.window!=99)
      keep.idx2 = keep.idx[rm.idx]
      full.window2 = full.window[keep.idx]
      full.window2 = full.window2[rm.idx]
      sel.idx = keep.idx[which(full.window2 < p.val)]
      window.res[sel.idx] = 1
    }
  }
  names(window.res) = wd.check
  return(list(cp.idx = cp.idx,
              window.res = window.res))
}

detect.overlap = function(x){
  cp.idx = x$cp.idx
  window.info = as.numeric(names(x$window.res))
  window.idx = window.info[which(x$window.res == 1)]
  if(length(window.idx) > 0){
    all.idx = c(cp.idx + window.idx, cp.idx)
  }else{
    all.idx = cp.idx
  }
  all.idx = base::sort(all.idx)
  return(all.idx)
  
}

lm.nb.new = function(gamma.v, y, h0,h1, sigma2, X){
  n.cov = ncol(X)
  seg.num = sum(gamma.v)
  seg.label = cumsum(gamma.v)
  coeff.sets = var.mat = intv.mat = vector("list",length = seg.num)
  for(i in 1:seg.num){
    intval.mat = matrix(0, nrow = n.cov, ncol = 2)
    y.seg = y[seg.label == i]
    y.seg = matrix(y.seg, ncol = 1)
    Xsub = X[seg.label == i, ]
    Xsub = matrix(Xsub, ncol = n.cov)
    beta_prior_var = diag(c(h0, h1))
    beta.bayes = solve(solve(beta_prior_var) + t(Xsub) %*% Xsub / sigma2,
                       t(Xsub) %*% y.seg / sigma2)
    var.mat[[i]] = solve(solve(beta_prior_var) + t(Xsub) %*% Xsub / sigma2)
    beta.bayes = as.numeric(beta.bayes)
    coeff.sets[[i]] = beta.bayes
    
    for(bb in 1:length(beta.bayes)){
      intval.mat[bb,] = c(qnorm(0.025, mean = beta.bayes[bb], sd = var.mat[[i]][bb, bb]),
                          qnorm(0.975, mean = beta.bayes[bb], sd = var.mat[[i]][bb, bb]) )
    }
    temp.mat = cbind(beta.bayes, intval.mat)
    colnames(temp.mat) = c("mean","2.5%", "97.5%")
    rownames(temp.mat) = c("intercept", "time")
    intv.mat[[i]] = temp.mat
  }
  ret.list = list(beta.est = coeff.sets,
                  beta.covar = var.mat,
                  beta.ci = intv.mat)
  return(ret.list)
}

get.sir.df = function(cp_res, dat){
  gamma_ppm = cp_res$gamma_ppm
  # get z vector #
  z_vec = cumsum(gamma_ppm)
  ret_df = data.frame(segment = z_vec, 
                      date = dat$date,
                      time = dat$time,
                      confirm = dat$C_obs,
                      removed = dat$R_obs,
                      active_infected = dat$I_obs,
                      stringsAsFactors = F )
  return(ret_df)
}

# function for run SIR model on each segement for a fixed gamma vector #
fit.sir.segment = function(rank_idx, 
                           sir_input_list = sir_input_list, 
                           pop = pop, 
                           recovery.rate = 0.1,
                           T_pred = 3){
  df = sir_input_list[[rank_idx]]
  #I_obs = df$active_infected
  #R_obs = df$removed
  z_vec = df$segment - 1
  
  cat("estimating R0 by segment \n")
  C_obs = df$confirm
  seg_sir = sir_c(C_obs, 
                  N = pop,
                  gamma0 = recovery.rate,
                  z = z_vec, 
                  T_fin = T_pred,
                  store = FALSE);
  
  # prediction using only the last segment #
  last_seg_idx = z_vec == max(z_vec)
  C_last = C_obs[last_seg_idx]
  z_last = rep(0, sum(last_seg_idx))
  
  last_sir = sir_c(C_last, 
                  N = pop,
                  gamma0 = recovery.rate,
                  z = z_last, 
                  T_fin = T_pred,
                  store = FALSE);
  
  pred_df_last = data.frame(N_pred_mean = last_sir[["N_pred_mean"]], 
                       N_pred_lower = last_sir[["N_pred_lwr"]], 
                       N_pred_upper = last_sir[["N_pred_upp"]])
  rm(last_sir)
  
  all_sir = sir_c(C_obs, 
                  N = pop,
                  gamma0 = recovery.rate,
                  z = rep(0, length(z_vec)), 
                  T_fin = T_pred,
                  store = FALSE);
  
  pred_df_full = data.frame(N_pred_mean = all_sir[["N_pred_mean"]], 
                            N_pred_lower = all_sir[["N_pred_lwr"]], 
                            N_pred_upper = all_sir[["N_pred_upp"]])
  rm(all_sir)
  
  # put together the SIR estimated R0 #
  R0_mcmc = seg_sir$R0_store
  R0_est = round(apply(R0_mcmc, 2, mean), 3)
  R0_CI = round(apply(R0_mcmc, 2, quantile, c(0.025, 0.975)), 3)
  
  # T-test for R0 between segments #
  z_idx = unique(z_vec) + 1
  compare_idx = cbind(z_idx[-length(z_idx)], z_idx[-1])
  pval_paird = apply(compare_idx, 1, 
                    function(x){wil.test = wilcox.test(x = R0_mcmc[, x[1]], y = R0_mcmc[, x[2]], paired = T);
                                pval = format(round(wil.test$p.value, 6), nsmall = 6); return(pval)})
  
  full_sir = data.frame(mean = R0_est, lower = R0_CI[1,], upper = R0_CI[2,], 
                        stringsAsFactors = F, check.names = F)
  full_sir = list("full_sir" = full_sir, 
                  "p_values" = pval_paird, 
                  "prediction_last" = pred_df_last, 
                  "prediction_full" = pred_df_full)
  rm(seg_sir, R0_mcmc)
  return(full_sir)
}


