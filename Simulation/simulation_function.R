get_info_simulation = function(result, digits = 3) {
  # Get info for ecrem
  coefs_ecrem = rbindlist(lapply(result,FUN = function(x){data.table(t((x$coefs_real_final)))}))
  ave_ecrem = apply(coefs_ecrem,2, mean )
  coefs_std = t((t(coefs_ecrem) - result[[1]]$coefs_real_gt_1)^2)
  rmse_ecrem = sqrt(1/nrow(coefs_std)*colSums(coefs_std))
  alpha = rbindlist(lapply(result,FUN = function(x){
    x$coefs_upper_tmp = x$coefs_real_final + 2*sqrt(diag(x$vcov))
    x$coefs_lower_tmp = x$coefs_real_final - 2*sqrt(diag(x$vcov))
    data.table(t((x$coefs_real_gt_1 <x$coefs_upper_tmp) &  (x$coefs_real_gt_1 >x$coefs_lower_tmp)))}))
  cp_ecrem = colSums(alpha)/nrow(alpha)
  # Get info for rem
  coefs_rem = rbindlist(lapply(result,FUN = function(x){data.table(t((x$coefs_unchange)))}))
  ave_rem = apply(coefs_rem,2, mean )
  coefs_std = t((t(coefs_rem) - result[[1]]$coefs_real_gt_1)^2)
  rmse_rem = sqrt(1/nrow(coefs_std)*colSums(coefs_std))
  if(is.null(result[[1]]$upper_ci_unchanged)){
    alpha = rbindlist(lapply(result,FUN = function(x){data.table(t((x$coefs_real_gt_1 <x$upper_ci_gt) &  (x$coefs_real_gt_1 >x$lower_ci_gt)))}))
  } else {
    alpha = rbindlist(lapply(result,FUN = function(x){data.table(t((x$coefs_real_gt_1 <x$upper_ci_unchanged) &  (x$coefs_real_gt_1 >x$lower_ci_unchanged)))}))
  }
  cp_rem = colSums(alpha)/nrow(alpha)
  res = cbind(coefs = result[[1]]$coefs_real_gt_1,ave_ecrem,rmse_ecrem,cp_ecrem,ave_rem,rmse_rem,cp_rem)
  rownames(res) = str_to_title(names(result[[1]]$coefs_real_gt_1))
  tmp = rownames(res)
  tmp = gsub("\\_",replacement = " ",tmp)
  rownames(res) = tmp
  res = round(res,digits = digits)
  class(res) = "result_simulation"
  return(res)
}

mcem_estimation_new = function(cl = NULL, K = 50,M = 100, beta_real_gt,beta_fake_gt,
                           comb_event_stream = comb_event_stream,
                           n_actors = n_actors, eps = 0.1,
                           starting_vals_real = NULL,
                           starting_vals_fake = NULL, exo_cov, 
                           exo_cat, seed, build, build_time = NULL){
  
  set.seed(seed)
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt!=0]
  betas_fake = matrix(nrow = K+1,ncol = length(which_covs_fake))
  betas_real = matrix(nrow = K+1,ncol = length(which_covs_real))
  # Set starting parameters
  if(is.null(starting_vals_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    test_unchange = count_rem_undirected_list(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    starting_vals_real = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                           y = test_unchange$status,
                           offset = test_unchange$offset, family = poisson(),method = 3)$coefficients
  }
  if(is.null(starting_vals_fake)){
    # For the fake parameters, we sample a random starting point
    starting_vals_fake = rnorm(n =  length(which_covs_fake))
  }
  
  betas_fake[1,] =  starting_vals_fake
  betas_real[1,] = starting_vals_real
  
  for(k in 1:K) {
    
    
    if(is.null(cl)){
      result = lapply(1:(M), e_step_compl_list, 
                         comb_event_stream = comb_event_stream, 
                         n_actors = n_actors, 
                         beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                         which_covs_real = which_covs_real, 
                         which_covs_fake = which_covs_fake,seed = k + seed, 
                         exo_cov = exo_cov, exo_cat = exo_cat)
    } else {
      result = parLapply(cl, 1:(M), e_step_compl_list, 
                         comb_event_stream = comb_event_stream, 
                         n_actors = n_actors, 
                         beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                         which_covs_real = which_covs_real, 
                         which_covs_fake = which_covs_fake,seed = k + seed, 
                         exo_cov = exo_cov, exo_cat = exo_cat)
    }
    
    # result = parLapply(cl, 1:(M), e_step_compl_new, 
    #                    comb_event_stream = comb_event_stream, 
    #                    n_actors = n_actors, 
    #                    beta_real =  beta_real_gt[beta_real_gt != 0],beta_fake = beta_fake_gt[beta_fake_gt != 0],
    #                    which_covs_real = which_covs_real, 
    #                    which_covs_fake = which_covs_fake,seed = k + seed, 
    #                    exo_cov = exo_cov, exo_cat = exo_cat)
    
    # result = parLapply(cl, 1:(M), get_draw_both,
    #                    comb_event_stream = comb_event_stream,
    #                    n_actors = n_actors,
    #                    beta_real =mod_real_gt$coefficients, beta_fake = mod_fake_gt$coefficients,
    #                    which_covs_real= which_covs_real,
    #                    which_covs_fake = which_covs_fake,seed = k + seed - 123,
    #                    exo_cov = exo_cov, exo_cat = exo_cat)
    # result = do.call(what = "cbind", args = as.list(result))
    # 
    # #
    # comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
    #                                   n_actors = n_actors,
    #                                   beta_fake = beta_real,
    #                                   beta_real = beta_fake,
    #                                   which_covs_fake = which_covs_real,
    #                                   which_covs_real=which_covs_fake ,
    #                                   seed = seed*m,
    #                                   exo_cov = exo_cov,
    #                                   exo_cat = exo_cov)
    # 
    # comb_event_stream$average = rowSums(result)/ncol(result)
    # comb_event_stream$index = 1:length(comb_event_stream$side_a)
    # ggplot(data = comb_event_stream, aes( y = average, x = factor(ind_real)))+
    #   geom_boxplot()

    
    
    # mean(comb_event_stream$average[comb_event_stream$ind_real == 0])
    # 
    result = rbindlist(result)
    # build_time = comb_event_stream$times[100]
    # Should there be some build-up time?
    if(build) {
      result$include =  (result$time >= build_time) +  (result$time_end > build_time)
      result$time[result$include == 1] = build_time
      result = result[include!= 0]
    }
    
    
    count_real = result[add_var == 1]
    count_fake = result[add_var == 0]
    
    # mod_real =   fastglm(x =   as.matrix(count_real[,.(intercept,degree_sum,degree_abs,triangle,repetition,cov_cont,cov_cat)]), 
    #                      y = count_real$status, 
    #                      offset = count_real$offset, family = poisson(),method = 3)
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  +
    #                  cov_cat + offset(log(time_end-time)), data = count_real[time_end>  comb_event_stream$times[200]],family = poisson)
    
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    
    mod_fake =   fastglm(x =   as.matrix(count_fake[,which_covs_fake, with = F]),
                         y = count_fake$status,
                         offset = count_fake$offset, family = poisson(),method = 3)
    
    betas_fake[k+1,] = mod_fake$coefficients
    betas_real[k+1,] = mod_real$coefficients
    cat("Iteration ", k," completed\n")
    cat("Change =", sum(abs(betas_real[k+1,] - betas_real[k,] ))/(sum(abs( betas_real[k,] ))),"\n")
    cat("Paramerers real = ", mod_real$coefficients, "\n")
    cat("Paramerers fake = ", mod_fake$coefficients,"\n")
    
    if(sum(abs(betas_real[k+1,] - betas_real[k,] ))/(sum(abs( betas_real[k,] )))< eps) {
      
      score_per_simulation = lapply(1:(M), score_function, 
                                             data = count_real, model = mod_real, which_covs_real = which_covs_real)
      if(is.null(cl)){
        hessian_per_simulation = lapply(1:(M), hessian_function_new, 
                                        data = count_real, model = mod_real, which_covs_real = which_covs_real)
      } else {
        hessian_per_simulation = parLapply(cl, 1:(M), hessian_function_new, 
                                           data = count_real, model = mod_real, which_covs_real = which_covs_real)
      }
      
      
      mean = 1/M*Reduce(x = score_per_simulation, "+")
      new_info = -1/M*Reduce(x = hessian_per_simulation, "+") +  1/M*Reduce(x = lapply(score_per_simulation, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
      # new_hessian = 1/M*Reduce(x = hessian_per_simulation, "+") + 
      #   1/M*Reduce(x =  lapply(score_per_simulation, FUN = function(x){x%*% t(x)}), "+") +
      #   1/M*Reduce(x = score_per_simulation, "+")%*% t(Reduce(x = score_per_simulation, "+"))
      # 
      vcov =  solve(new_info)
      break
    }
  }
  return(list(  betas_fake= betas_fake,betas_real = betas_real, vcov = vcov ))
}


mcem_estimation = function(cl = NULL, K = 50,M = 100, beta_real_gt,beta_fake_gt,
                               comb_event_stream = comb_event_stream,
                               n_actors = n_actors, eps = 0.1,
                               starting_vals_real = NULL,
                               starting_vals_fake = NULL, exo_cov, 
                               exo_cat, seed, build, build_time = NULL){
  
  set.seed(seed)
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt!=0]
  betas_fake = matrix(nrow = K+1,ncol = length(which_covs_fake))
  betas_real = matrix(nrow = K+1,ncol = length(which_covs_real))
  # Set starting parameters
  if(is.null(starting_vals_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    test_unchange = count_rem_undirected_new(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    starting_vals_real = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                                 y = test_unchange$status,
                                 offset = test_unchange$offset, family = poisson(),method = 3)$coefficients
  }
  if(is.null(starting_vals_fake)){
    # For the fake parameters, we sample a random starting point
    starting_vals_fake = rnorm(n =  length(which_covs_fake))
  }
  
  betas_fake[1,] =  starting_vals_fake
  betas_real[1,] = starting_vals_real
  
  for(k in 1:K) {
    
    
    if(is.null(cl)){
      result = lapply(1:(M), e_step_compl_new, 
                      comb_event_stream = comb_event_stream, 
                      n_actors = n_actors, 
                      beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                      which_covs_real = which_covs_real, 
                      which_covs_fake = which_covs_fake,seed = k + seed, 
                      exo_cov = exo_cov, exo_cat = exo_cat)
    } else {
      result = parLapply(cl, 1:(M), e_step_compl_new, 
                         comb_event_stream = comb_event_stream, 
                         n_actors = n_actors, 
                         beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                         which_covs_real = which_covs_real, 
                         which_covs_fake = which_covs_fake,seed = k + seed, 
                         exo_cov = exo_cov, exo_cat = exo_cat)
    }
    
    # result = parLapply(cl, 1:(M), e_step_compl_new, 
    #                    comb_event_stream = comb_event_stream, 
    #                    n_actors = n_actors, 
    #                    beta_real =  beta_real_gt[beta_real_gt != 0],beta_fake = beta_fake_gt[beta_fake_gt != 0],
    #                    which_covs_real = which_covs_real, 
    #                    which_covs_fake = which_covs_fake,seed = k + seed, 
    #                    exo_cov = exo_cov, exo_cat = exo_cat)
    
    # result = parLapply(cl, 1:(M), get_draw_both,
    #                    comb_event_stream = comb_event_stream,
    #                    n_actors = n_actors,
    #                    beta_real =mod_real_gt$coefficients, beta_fake = mod_fake_gt$coefficients,
    #                    which_covs_real= which_covs_real,
    #                    which_covs_fake = which_covs_fake,seed = k + seed - 123,
    #                    exo_cov = exo_cov, exo_cat = exo_cat)
    # result = do.call(what = "cbind", args = as.list(result))
    # 
    # #
    # comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
    #                                   n_actors = n_actors,
    #                                   beta_fake = beta_real,
    #                                   beta_real = beta_fake,
    #                                   which_covs_fake = which_covs_real,
    #                                   which_covs_real=which_covs_fake ,
    #                                   seed = seed*m,
    #                                   exo_cov = exo_cov,
    #                                   exo_cat = exo_cov)
    # 
    # comb_event_stream$average = rowSums(result)/ncol(result)
    # comb_event_stream$index = 1:length(comb_event_stream$side_a)
    # ggplot(data = comb_event_stream, aes( y = average, x = factor(ind_real)))+
    #   geom_boxplot()
    
    
    
    # mean(comb_event_stream$average[comb_event_stream$ind_real == 0])
    # 
    result = rbindlist(result)
    # build_time = comb_event_stream$times[100]
    # Should there be some build-up time?
    if(build) {
      result$include =  (result$time >= build_time) +  (result$time_end > build_time)
      result$time[result$include == 1] = build_time
      result = result[include!= 0]
    }
    
    
    count_real = result[add_var == 1]
    count_fake = result[add_var == 0]
    
    # mod_real =   fastglm(x =   as.matrix(count_real[,.(intercept,degree_sum,degree_abs,triangle,repetition,cov_cont,cov_cat)]), 
    #                      y = count_real$status, 
    #                      offset = count_real$offset, family = poisson(),method = 3)
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  +
    #                  cov_cat + offset(log(time_end-time)), data = count_real[time_end>  comb_event_stream$times[200]],family = poisson)
    
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    
    mod_fake =   fastglm(x =   as.matrix(count_fake[,which_covs_fake, with = F]),
                         y = count_fake$status,
                         offset = count_fake$offset, family = poisson(),method = 3)
    
    betas_fake[k+1,] = mod_fake$coefficients
    betas_real[k+1,] = mod_real$coefficients
    cat("Iteration ", k," completed\n")
    cat("Change =", sum(abs(betas_real[k+1,] - betas_real[k,] ))/(sum(abs( betas_real[k,] ))),"\n")
    cat("Paramerers real = ", mod_real$coefficients, "\n")
    cat("Paramerers fake = ", mod_fake$coefficients,"\n")
    
    if(sum(abs(betas_real[k+1,] - betas_real[k,] ))/(sum(abs( betas_real[k,] )))< eps) {
      
      score_per_simulation = lapply(1:(M), score_function, 
                                    data = count_real, model = mod_real, which_covs_real = which_covs_real)
      if(is.null(cl)){
        hessian_per_simulation = lapply(1:(M), hessian_function_new, 
                                        data = count_real, model = mod_real, which_covs_real = which_covs_real)
      } else {
        hessian_per_simulation = parLapply(cl, 1:(M), hessian_function_new, 
                                           data = count_real, model = mod_real, which_covs_real = which_covs_real)
      }
      
      
      mean = 1/M*Reduce(x = score_per_simulation, "+")
      new_info = -1/M*Reduce(x = hessian_per_simulation, "+") +  1/M*Reduce(x = lapply(score_per_simulation, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
      # new_hessian = 1/M*Reduce(x = hessian_per_simulation, "+") + 
      #   1/M*Reduce(x =  lapply(score_per_simulation, FUN = function(x){x%*% t(x)}), "+") +
      #   1/M*Reduce(x = score_per_simulation, "+")%*% t(Reduce(x = score_per_simulation, "+"))
      # 
      vcov =  solve(new_info)
      break
    }
  }
  return(list(  betas_fake= betas_fake,betas_real = betas_real, vcov = vcov ))
}


mcem_estimation_old = function(cl, K = 50,M = 100, beta_real_gt,beta_fake_gt,
                           comb_event_stream = comb_event_stream,
                           n_actors = n_actors, eps = 0.1,
                           starting_vals_real = NULL,
                           starting_vals_fake = NULL, exo_cov, exo_cat, seed){
  
  set.seed(seed)
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt!=0]
  betas_fake = matrix(nrow = K+1,ncol = length(which_covs_fake))
  betas_real = matrix(nrow = K+1,ncol = length(which_covs_real))
  # Set starting parameters
  if(is.null(starting_vals_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    test_unchange = count_rem_undirected(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    starting_vals_real = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                                 y = test_unchange$status,
                                 offset = test_unchange$offset, family = poisson(),method = 3)$coefficients
  }
  if(is.null(starting_vals_fake)){
    # For the fake parameters, we sample a random starting point
    starting_vals_fake = rnorm(n =  length(which_covs_fake))
  }
  
  betas_fake[1,] =  starting_vals_fake
  betas_real[1,] = starting_vals_real
  
  for(k in 1:K) {
    result = parLapply(cl, 1:(M), e_ste_compl, 
                       comb_event_stream = comb_event_stream, 
                       n_actors = n_actors, 
                       beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                       which_covs_real = which_covs_real, 
                       which_covs_fake = which_covs_fake,seed = k + seed, 
                       exo_cov = exo_cov, exo_cat = exo_cat)
    
    result = rbindlist(result)
    
    count_real = result[add_var == 1]
    count_fake = result[add_var == 0]
    
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    mod_fake =   fastglm(x =   as.matrix(count_fake[,which_covs_fake, with = F]),
                         y = count_fake$status,
                         offset = count_fake$offset, family = poisson(),method = 3)
    betas_fake[k+1,] = mod_fake$coefficients
    betas_real[k+1,] = mod_real$coefficients
    cat("Iteration ", k," completed\n")
    cat("Change =", sum(abs(betas_real[k+1,]  - betas_real[k,] )),"\n")
    cat("Paramerers real = ", mod_real$coefficients, "\n")
    cat("Paramerers fake = ", mod_fake$coefficients,"\n")
    
    if(sum(abs(betas_real[k+1,] - betas_real[k,] ))< eps) {
      break
    }
  }
  return(list(  betas_fake= betas_fake,betas_real = betas_real ))
}


e_ste_compl = function(m,comb_event_stream, n_actors, beta_real, beta_fake, 
                       which_covs_real , 
                       which_covs_fake ,seed = k, 
                       exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, 
            beta_real = beta_real, beta_fake = beta_fake, 
            which_covs_real = which_covs_real, which_covs_fake=which_covs_fake ,seed = seed*m,exo_cov = exo_cov, exo_cat = exo_cov)
  
  return(rbind(count_rem_undirected(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 1, ind = m), 
               count_rem_undirected(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 0, ind = m)))
}

e_step_compl_new = function(m,comb_event_stream,
                            n_actors, beta_real, beta_fake, 
                       which_covs_real , 
                       which_covs_fake ,seed = k, 
                       exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                    n_actors = n_actors, 
                                    beta_real = beta_real,
                                    beta_fake = beta_fake, 
                                    which_covs_real = which_covs_real,
                                    which_covs_fake=which_covs_fake ,
                                    seed = seed*m,
                                    exo_cov = exo_cov, 
                                    exo_cat = exo_cat)
  
  return(rbind(count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 1, ind = m), 
               count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 0, ind = m)))
}


e_step_compl_list = function(m,comb_event_stream,
                            n_actors, beta_real, beta_fake, 
                            which_covs_real , 
                            which_covs_fake ,seed = k, 
                            exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step_list(comb_event_stream = comb_event_stream,
                                    n_actors = n_actors, 
                                    beta_real = beta_real,
                                    beta_fake = beta_fake, 
                                    which_covs_real = which_covs_real,
                                    which_covs_fake=which_covs_fake ,
                                    seed = seed*m,
                                    exo_cov = exo_cov, 
                                    exo_cat = exo_cat)
  
  return(rbind(count_rem_undirected_list(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 1, ind = m), 
               count_rem_undirected_list(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 0, ind = m)))
}


e_step_compl_both = function(m,comb_event_stream,
                            n_actors, beta_real, beta_fake, 
                            which_covs_real , 
                            which_covs_fake ,seed = k, 
                            exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step_both(comb_event_stream = comb_event_stream,
                                    n_actors = n_actors, 
                                    beta_real = beta_real,
                                    beta_fake = beta_fake, 
                                    which_covs_real = which_covs_real,
                                    which_covs_fake=which_covs_fake ,
                                    seed = seed*m,
                                    exo_cov = exo_cov, 
                                    exo_cat = exo_cov)
  
  return(rbind(count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 1, ind = m), 
               count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat, add_var = 0, ind = m)))
}
get_draw = function(m,comb_event_stream,
                            n_actors, beta_real, beta_fake, 
                            which_covs_real , 
                            which_covs_fake ,seed = k, 
                            exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                    n_actors = n_actors, 
                                    beta_real = beta_real,
                                    beta_fake = beta_fake, 
                                    which_covs_real = which_covs_real,
                                    which_covs_fake=which_covs_fake ,
                                    seed = seed*m,
                                    exo_cov = exo_cov, 
                                    exo_cat = exo_cov)
  
  return(comb_event_stream$ind)
}
get_draw_both = function(m,comb_event_stream,
                    n_actors, beta_real, beta_fake, 
                    which_covs_real , 
                    which_covs_fake ,seed = k, 
                    exo_cov, exo_cat) {
  # set.seed(seed*m)
  # comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, beta_real = beta_real, beta_fake = beta_fake)
  comb_event_stream = st_e_step_both(comb_event_stream = comb_event_stream,
                                    n_actors = n_actors, 
                                    beta_real = beta_real,
                                    beta_fake = beta_fake, 
                                    which_covs_real = which_covs_real,
                                    which_covs_fake=which_covs_fake ,
                                    seed = seed*m,
                                    exo_cov = exo_cov, 
                                    exo_cat = exo_cov)
  
  return(comb_event_stream$ind)
}

sem_estimation = function(K = 50, beta_real_gt,beta_fake_gt,comb_event_stream = comb_event_stream,n_actors = n_actors, 
                          starting_vals_real = NULL,
                          starting_vals_fake = NULL,
                          exo_cov,exo_cat, seed){
  set.seed(seed)
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt!=0]
  betas_fake = matrix(nrow = K+1,ncol = length(which_covs_fake))
  betas_real = matrix(nrow = K+1,ncol = length(which_covs_real))
  # Set starting parameters
  if(is.null(starting_vals_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    test_unchange = count_rem_undirected(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    starting_vals_real = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                                 y = test_unchange$status,
                                 offset = test_unchange$offset, family = poisson(),method = 3)$coefficients
  }
  if(is.null(starting_vals_fake)){
    # For the fake parameters, we sample a random starting point
    starting_vals_fake = rnorm(n =  length(which_covs_fake))
  }
  # Set the first parameters
  betas_fake[1,] =  starting_vals_fake
  betas_real[1,] = starting_vals_real
  
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step(comb_event_stream = comb_event_stream,n_actors = n_actors, 
                                  beta_real = betas_real[k,], beta_fake = betas_fake[k,], 
                                  which_covs_real = which_covs_real, 
                                  which_covs_fake = which_covs_fake,seed = seed+k,exo_cov = exo_cov,exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    mod_fake =   fastglm(x =   as.matrix(count_fake[,which_covs_fake, with = F]),
                         y = count_fake$status,
                         offset = count_fake$offset, family = poisson(),method = 3)
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = mod_fake$coefficients
    betas_real[k+1,] = mod_real$coefficients
    cat("Iteration ", k," completed\n")
  }
  return(list(betas_real = betas_real,betas_fake = betas_fake))
}


sem_estimation_new = function(K = 50, beta_real_gt,beta_fake_gt,comb_event_stream = comb_event_stream,n_actors = n_actors, 
                          starting_vals_real = NULL,
                          starting_vals_fake = NULL,
                          exo_cov,exo_cat, seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt!=0]
  betas_fake = matrix(nrow = K+1,ncol = length(which_covs_fake))
  betas_real = matrix(nrow = K+1,ncol = length(which_covs_real))
  # Set starting parameters
  if(is.null(starting_vals_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    test_unchange = count_rem_undirected_new(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    starting_vals_real = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                                 y = test_unchange$status,
                                 offset = test_unchange$offset, family = poisson(),method = 3)$coefficients
  }
  if(is.null(starting_vals_fake)){
    # For the fake parameters, we sample a random starting point
    starting_vals_fake = rnorm(n =  length(which_covs_fake))
  }
  # Set the first parameters
  betas_fake[1,] =  starting_vals_fake
  betas_real[1,] = starting_vals_real
  
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = betas_real[k,],
                                      beta_fake = betas_fake[k,], 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    mod_fake =   fastglm(x =   as.matrix(count_fake[,which_covs_fake, with = F]),
                         y = count_fake$status,
                         offset = count_fake$offset, family = poisson(),method = 3)
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = mod_fake$coefficients
    betas_real[k+1,] = mod_real$coefficients
    cat("Iteration ", k," completed\n")
  }
  return(list(betas_real = betas_real,betas_fake = betas_fake))
}


mi_estimation = function(K = 50, mod_real,mod_fake,formula_real, formula_fake,
                         comb_event_stream = comb_event_stream,
                         n_actors = n_actors, 
                         exo_cov,exo_cat,
                         seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
  which_covs_real =names(coef(mod_real))
  which_covs_real[1] = "Intercept"
  which_covs_fake ="Intercept"
  betas_fake = matrix(nrow = K+1,ncol = length(coef(mod_fake)))
  betas_real = matrix(nrow = K+1,ncol = length(coef(mod_real)))
  # Set the first parameters
  betas_fake[1,] =  coef(mod_fake)
  betas_real[1,] = coef(mod_real)
  
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)

    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    betas_real[k+1,] = rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_real$coefficients =    betas_real[k+1,]
    mod_fake$coefficients = betas_fake[k+1,]
    cat("Iteration ", k," completed\n")
  }
  return(list(betas_real = betas_real,betas_fake = betas_fake))
}


mi_estimation_with_se_alt = function(K = 50, K_2 = 10, mod_real,mod_fake,formula_real, formula_fake,
                         comb_event_stream = comb_event_stream,pi_start,
                         n_actors = n_actors, 
                         exo_cov,exo_cat,
                         seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
  which_covs_real =names(coef(mod_real))
  which_covs_real[1] = "Intercept"
  which_covs_fake ="Intercept"
  betas_fake = matrix(nrow = K+1,ncol = length(coef(mod_fake)))
  betas_real = matrix(nrow = K+1,ncol = length(coef(mod_real)))
  # Set the first parameters
  betas_fake[1,] =  coef(mod_fake)
  betas_real[1,] = coef(mod_real)
  list_model_real = list()
  list_model_fake = list()
  pi  = c()
  pi[1] = pi_start
  # Warum up ----
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step_alt(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,pi = pi[1],
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    pi[k+1] = mean(comb_event_stream$ind)
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)
    
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    betas_real[k+1,] = rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_real$coefficients =    betas_real[k+1,]
    mod_fake$coefficients = betas_fake[k+1,]
    cat("Iteration ", k," completed\n")
  }
  pi_new  = c()
  pi_new[1] = pi[k]
  # Standard error and estiation ----
  for(k in 1:K_2) {
    # Simulation Step
    comb_event_stream = st_e_step_alt(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,pi = pi_start,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    pi_new[k+1] = mean(comb_event_stream$ind)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)
    
    list_model_real[[k]] = mod_real
    list_model_fake[[k]] = mod_fake

    mod_real$coefficients =   rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_fake$coefficients =  rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    cat("Iteration ", k," completed\n")
  }
  
  est_coefs_fake = lapply(list_model_fake,function(x) return(coef(x)))
  mean_fake = 1/length(list_model_real)*Reduce("+",est_coefs_fake)
  
  est_coefs = lapply(list_model_real,function(x) return(coef(x)))
  est_var = lapply(list_model_real,function(x) return(vcov(x)))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/K_2*Reduce("+",est_coefs)
  B_bar =1/(K_2-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/K_2*Reduce("+",est_var)
  final_var_alt = V_bar + (1+1/K_2)*B_bar

  return(list(betas_real = betas_real,betas_fake = betas_fake,
              vcov = final_var_alt, final_coefs_real = mean, final_coefs_fake = mean_fake, B = B_bar, W = V_bar))
}

mi_estimation_with_se = function(K = 50, K_2 = 10, mod_real,mod_fake,formula_real, formula_fake,
                                 comb_event_stream = comb_event_stream,
                                 n_actors = n_actors, 
                                 exo_cov,exo_cat,
                                 seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
  which_covs_real =names(coef(mod_real))
  which_covs_real[1] = "Intercept"
  which_covs_fake ="Intercept"
  betas_fake = matrix(nrow = K+1,ncol = length(coef(mod_fake)))
  betas_real = matrix(nrow = K+1,ncol = length(coef(mod_real)))
  # Set the first parameters
  betas_fake[1,] =  coef(mod_fake)
  betas_real[1,] = coef(mod_real)
  list_model_real = list()
  list_model_fake = list()
  # Warum up ----
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)
    
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    betas_real[k+1,] = rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_real$coefficients =    betas_real[k+1,]
    mod_fake$coefficients = betas_fake[k+1,]
    cat("Iteration ", k," completed\n")
  }
  percentage_real = c()
  # Standard error and estiation ----
  for(k in 1:K_2) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    percentage_real[k] = mean(comb_event_stream$ind)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)
    
    list_model_real[[k]] = mod_real
    list_model_fake[[k]] = mod_fake
    
    mod_real$coefficients =   rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_fake$coefficients =  rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    cat("Iteration ", k," completed\n")
  }
  
  est_coefs_fake = lapply(list_model_fake,function(x) return(coef(x)))
  mean_fake = 1/length(list_model_real)*Reduce("+",est_coefs_fake)
  
  est_coefs = lapply(list_model_real,function(x) return(coef(x)))
  est_var = lapply(list_model_real,function(x) return(vcov(x)))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/K_2*Reduce("+",est_coefs)
  B_bar =1/(K_2-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/K_2*Reduce("+",est_var)
  final_var_alt = V_bar + (1+1/K_2)*B_bar
  
  return(list(betas_real = betas_real,betas_fake = betas_fake,
              vcov = final_var_alt, final_coefs_real = mean, final_coefs_fake = mean_fake, B = B_bar, W = V_bar, 
              percentage_real = mean(percentage_real)))
}

mi_estimation_se = function(K = 50, mod_real_start,mod_fake_start,formula_real, formula_fake,
                         comb_event_stream = comb_event_stream,
                         n_actors = n_actors, 
                         exo_cov,exo_cat,
                         seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
  mod_fake = mod_fake_start
  mod_real = mod_real_start
  which_covs_real =names(coef(mod_real))
  which_covs_real[1] = "Intercept"
  which_covs_fake ="Intercept"
  betas_fake = matrix(nrow = K+1,ncol = length(coef(mod_fake)))
  betas_real = matrix(nrow = K+1,ncol = length(coef(mod_real)))
  # Set the first parameters

  betas_fake[1,] =  coef(mod_fake)
  betas_real[1,] = coef(mod_real)
  model_real_list = list()
  model_fake_list = list()
  for(k in 1:K) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real = mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    # comb_event_stream[which(comb_event_stream$ind == comb_event_stream$ind_real)]
    # Given the simulated data, we carry out the complete-case analysis
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    # Should there be some build-up time?
    if(build) {
      count_real$include =  (count_real$time >= build_time) +  (count_real$time_end > build_time)
      count_real$time[count_real$include == 1] = build_time
      count_real = count_real[include!= 0]
      
      count_fake$include =  (count_fake$time >= build_time) +  (count_fake$time_end > build_time)
      count_fake$time[count_fake$include == 1] = build_time
      count_fake = count_fake[include!= 0]
      
    }
    
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    # if there are no fake events stop the loop
    if(nrow(count_fake) == 0){
      return(mod_real)
    }
    
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake)
    
    # mod_real = glm(status ~ degree_sum +degree_abs+ triangle +repetition +cov_cont  + 
    #                  cov_cat + offset(log(time_end-time)), data = count_real,family = poisson)
    # mod_fake = glm(status ~ 1 + offset(log(time_end-time)), data = count_fake,family = poisson)
    # Update the covariates 
    # beta_fake[1] = mod_fake$coefficients
    betas_fake[k+1,] = rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    betas_real[k+1,] = rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    mod_real$coefficients =    betas_real[k+1,]
    mod_fake$coefficients = betas_fake[k+1,]
    model_real_list[[k]] = mod_real
    model_fake_list[[k]] = mod_fake
    cat("Iteration ", k," completed\n")
  }
  
  est_coefs_fake = lapply(model_fake_list,function(x) return(coef(x)))
  mean_fake = 1/length(model_real_list)*Reduce("+",est_coefs_fake)
  
  est_coefs = lapply(model_real_list,function(x) return(coef(x)))
  est_var = lapply(model_real_list,function(x) return(vcov(x)))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/length(model_real_list)*Reduce("+",est_coefs)
  B_bar =1/(length(model_real_list)-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/length(model_real_list)*Reduce("+",est_var)
  B = B_bar
  W = V_bar
  final_var_alt = V_bar + B_bar
  
  return(list(vcov = final_var_alt, final_coefs_real = mean, final_coefs_fake = mean_fake))
}



simulation_complete = function(x, n_actors = 40,K = 50,exo_cov, exo_cat,
                               n = 500, beta_real_gt, beta_fake_gt) {
  set.seed(x)
  # exo_cov = rnorm(n = n_actors)
  # exo_cat = sample(size  = n_actors,x = c(1,2,3),replace = T)
  # 
  test_undirected_real = simulation_rem_undirected_new(beta = beta_real_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  test_undirected_fake = simulation_rem_undirected_new(beta = beta_fake_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  
  # Check out the Ground Truth
  test_count_undirected_gt = count_rem_undirected_new(event_data = test_undirected_real$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake =names(beta_fake_gt)[beta_fake_gt!=0]
  
  mod_real_gt =   fastglm(x =   as.matrix(test_count_undirected_gt[,which_covs_real, with = F]),
                          y = test_count_undirected_gt$status,
                          offset = test_count_undirected_gt$offset, family = poisson(),method = 3)
  
  
  end_real = max(test_undirected_real$data$times)
  test_undirected_fake$data = test_undirected_fake$data[times< end_real] 
  test_count_fake_gt = count_rem_undirected_new(event_data = test_undirected_fake$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  if(nrow(test_count_fake_gt)>0){
    mod_fake_gt =   fastglm(x =   as.matrix(test_count_fake_gt[,which_covs_fake, with = F]),
                            y = test_count_fake_gt$status,
                            offset = test_count_fake_gt$offset, family = poisson(),method = 3)
    
  } else {
    mod_fake_gt =  mod_real_gt
  }
  
  
  test_undirected_real$data$ind_real = 1
  test_undirected_fake$data$ind_real = 0
  comb_event_stream = rbind(test_undirected_real$data, test_undirected_fake$data)
  comb_event_stream = comb_event_stream[order(times)]
  comb_event_stream$ind = F
  comb_event_stream$weight = 1
  test_unchange = count_rem_undirected_new(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  mod_unchange = fastglm(x =   as.matrix(test_unchange[,which_covs_real, with = F]),
                         y = test_unchange$status,
                         offset = test_unchange$offset, family = poisson(),method = 3)
  
  rmse_naiv_1 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- beta_real_gt[beta_real_gt!= 0])^2))
  rmse_naiv_2 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- mod_real_gt$coefficients)^2))
  
  result_sem = sem_estimation_new(K = K,beta_real_gt = beta_real_gt,starting_vals_fake = 0,
                                  beta_fake_gt = beta_fake_gt[1], 
                                  comb_event_stream =comb_event_stream,
                                  n_actors =n_actors,exo_cov =exo_cov,
                                  exo_cat = exo_cat,seed = 123, 
                                  build = F,build_time = comb_event_stream$times[100])
  if(class(result_sem) == "fastglm"){
    coefs_real_final =result_sem$coefficients
    coefs_fake_final = -1000
    
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    se_tmp = list(vcov =result_sem$se)
    
    res = list(rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*se_tmp$vcov, 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*mod_fake_gt$se, 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*mod_fake_gt$se)
  } else {
    coefs_real_final = apply(result_sem$betas_real[(K-20):K,],MARGIN = 2,mean)
    coefs_fake_final = mean(result_sem$betas_fake[(K-20):K])
    
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    se_tmp = mi_se_estimation(M = 20,beta_real_gt = coefs_real_final, 
                              which_covs_real = which_covs_real, 
                              which_covs_fake = which_covs_fake[1], 
                              beta_fake_gt = coefs_fake_final,
                              comb_event_stream = comb_event_stream, 
                              n_actors =n_actors,exo_cov =exo_cov,
                              exo_cat = exo_cat,seed = 123)
    
    res = list(rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*se_tmp$vcov, 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*mod_real_gt$se, 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*mod_real_gt$se)
  }
  
  
  return(res)
}


simulation_complete_final = function(x, n_actors = 40,K = 50,K_2 = 20, exo_cov, exo_cat,
                               n = 500, beta_real_gt, beta_fake_gt) {
  set.seed(x)
  # exo_cov = rnorm(n = n_actors)
  # exo_cat = sample(size  = n_actors,x = c(1,2,3),replace = T)
  # 
  test_undirected_real = simulation_rem_undirected_new(beta = beta_real_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  test_undirected_fake = simulation_rem_undirected_new(beta = beta_fake_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  
  # Check out the Ground Truth
  test_count_undirected_gt = count_rem_undirected_new(event_data = test_undirected_real$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake =names(beta_fake_gt)[beta_fake_gt!=0]
  
  formula_real = formula(paste("status ~  offset(offset)+ ",paste(which_covs_real[-1],collapse = " + ")))
  # formula_fake = formula("status ~ cov_cont+ offset(offset)" )
  formula_fake = formula("status ~ offset(offset)" )
  
  mod_real_gt =   glm(formula_real,family = poisson(), data = test_count_undirected_gt)
  
  end_real = max(test_undirected_real$data$times)
  test_undirected_fake$data = test_undirected_fake$data[times< end_real] 
  test_count_fake_gt = count_rem_undirected_new(event_data = test_undirected_fake$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  if(nrow(test_count_fake_gt)>0){
    mod_fake_gt = glm(formula_fake,family = poisson(), data = test_count_fake_gt)
  } else {
    mod_fake_gt =  mod_real_gt
  }
  
  
  test_undirected_real$data$ind_real = 1
  test_undirected_fake$data$ind_real = 0
  comb_event_stream = rbind(test_undirected_real$data, test_undirected_fake$data)
  comb_event_stream = comb_event_stream[order(times)]
  comb_event_stream$ind = F
  comb_event_stream$weight = 1
  test_unchange = count_rem_undirected_new(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  mod_unchange = glm(formula_real,family = poisson(), data = test_unchange)
  rmse_naiv_1 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- beta_real_gt[beta_real_gt!= 0])^2))
  rmse_naiv_2 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- mod_real_gt$coefficients)^2))
  # Set the starting parameter of the measurement error model to 1 
  mod_fake_unchange = mod_fake_gt
  mod_fake_unchange$coefficients= c(-1)

  result_sem = mi_estimation_with_se(K = K,K_2 = K_2,mod_real = mod_unchange,
                             formula_real = formula_real,
                             formula_fake =formula_fake,
                             mod_fake = mod_fake_unchange,
                             comb_event_stream =comb_event_stream,
                             n_actors =n_actors,exo_cov =exo_cov,
                             exo_cat = exo_cat,seed = 23,
                             build = F,build_time = comb_event_stream$times[100])
  # result_alt = da_estimation_clustered_mgcv(K = 20, K_2 = 5,rand_prop_real = 0.5,
  #                                           beta_real_gt = beta_real_gt[-c(6,7)]!=0,
  #                                           separable = F,
  #                                           beta_fake_gt = beta_fake_gt[-c(6,7)]!=0,
  #                                           comb_event_stream = comb_event_stream,
  #                                           beg_formula_real = "status ~  offset(offset) +",
  #                                           beg_formula_fake= "status ~ offset(offset) +",
  #                                           n_actors = n_actors,
  #                                           seed = 1123,eps = 0.01,
  #                                           change_all = "all",
  #                                           exo_cat  = list("exo_cat" = exo_cat),
  #                                           exo_cov = list("exo_cov" =exo_cov),
  #                                           build = F)
  # # # debugonce(mi_estimation_with_se_alt)
  # result_sem = mi_estimation_with_se_alt(K = K,K_2 = K_2,mod_real = mod_unchange,pi_start = 0.5,
  #                                    formula_real = formula_real,
  #                                    formula_fake =formula_fake,
  #                                    mod_fake = mod_fake_unchange,
  #                                    comb_event_stream =comb_event_stream,
  #                                    n_actors =n_actors,exo_cov =exo_cov,
  #                                    exo_cat = exo_cat,seed = 23,
  #                                    build = F,build_time = comb_event_stream$times[100])
  # 
  # beta_fake_gt[-c(6,7)]!=0
 
  
  #browser()
  if("glm" %in% class(result_sem)){
    coefs_real_final =result_sem$coefficients
    coefs_fake_final = -1000
    
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    se_tmp = list(vcov = sqrt(diag(vcov(result_sem))))

    res = list(seed = x, 
               rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               vcov = vcov(mod_real_gt), 
               vcov_gt = vcov(mod_real_gt), 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*se_tmp$vcov, 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*se_tmp$vcov,
               percentage_real_gt = 0, 
               percentage_real_ours =0)
  } else {
    coefs_real_final = result_sem$final_coefs_real
    coefs_fake_final = result_sem$final_coefs_fake
   
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    res = list(seed = x, 
               rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               vcov = result_sem$vcov, 
               vcov_gt = vcov(mod_real_gt), 
               W = result_sem$W, 
               B = result_sem$B, 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*sqrt(diag(result_sem$vcov)), 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*sqrt(diag(result_sem$vcov)), 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_real_gt))), 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_real_gt))), 
               upper_ci_unchanged = mod_unchange$coefficients + qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_unchange))), 
               lower_ci_unchanged = mod_unchange$coefficients - qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_unchange))), 
               percentage_real_gt = nrow(test_undirected_real$data)/(nrow(test_undirected_real$data) + nrow(test_undirected_fake$data)), 
               percentage_real_ours = result_sem$percentage_real)
  
    }
  
  
  return(res)
}
simulation_complete_glm = function(x, n_actors = 40,K = 50,K_2 = 20, exo_cov, exo_cat,
                                   n = 500, beta_real_gt, beta_fake_gt) {
  set.seed(x)
  # exo_cov = rnorm(n = n_actors)
  # exo_cat = sample(size  = n_actors,x = c(1,2,3),replace = T)
  # 
  test_undirected_real = simulation_rem_undirected_new(beta = beta_real_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  test_undirected_fake = simulation_rem_undirected_new(beta = beta_fake_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  
  # Check out the Ground Truth
  test_count_undirected_gt = count_rem_undirected_new(event_data = test_undirected_real$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  which_covs_real =names(beta_real_gt)[beta_real_gt!=0]
  which_covs_fake =names(beta_fake_gt)[beta_fake_gt!=0]
  
  formula_real = formula(paste("status ~  offset(offset)+ ",paste(which_covs_real[-1],collapse = " + ")))
  formula_fake = formula("status ~  offset(offset)" )
  
  mod_real_gt =   glm(formula_real,family = poisson(), data = test_count_undirected_gt)
  
  end_real = max(test_undirected_real$data$times)
  test_undirected_fake$data = test_undirected_fake$data[times< end_real] 
  test_count_fake_gt = count_rem_undirected_new(event_data = test_undirected_fake$data,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  
  if(nrow(test_count_fake_gt)>0){
    mod_fake_gt = glm(formula_fake,family = poisson(), data = test_count_fake_gt)
  } else {
    mod_fake_gt =  mod_real_gt
  }
  
  
  test_undirected_real$data$ind_real = 1
  test_undirected_fake$data$ind_real = 0
  comb_event_stream = rbind(test_undirected_real$data, test_undirected_fake$data)
  comb_event_stream = comb_event_stream[order(times)]
  comb_event_stream$ind = F
  comb_event_stream$weight = 1
  test_unchange = count_rem_undirected_new(event_data = comb_event_stream,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
  mod_unchange = glm(formula_real,family = poisson(), data = test_unchange)
  rmse_naiv_1 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- beta_real_gt[beta_real_gt!= 0])^2))
  rmse_naiv_2 = sqrt(1/length(mod_unchange$coefficients)*sum((mod_unchange$coefficients- mod_real_gt$coefficients)^2))
  # Set the starting parameter of the measurement error model to 1 
  mod_fake_unchange = mod_fake_gt
  mod_fake_unchange$coefficients= -1
  result_sem = mi_estimation(K = K,mod_real = mod_unchange,
                             formula_real = formula_real, 
                             formula_fake =formula_fake, 
                             mod_fake = mod_fake_unchange, 
                             comb_event_stream =comb_event_stream,
                             n_actors =n_actors,exo_cov =exo_cov,
                             exo_cat = exo_cat,seed = 123, 
                             build = F,build_time = comb_event_stream$times[100])
  
  
  if("glm" %in% class(result_sem)){
    coefs_real_final =result_sem$coefficients
    coefs_fake_final = -1000
    
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    se_tmp = list(vcov = sqrt(diag(vcov(result_sem))))
    
    res = list(rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*se_tmp$vcov, 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*se_tmp$vcov, 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*se_tmp$vcov)
  } else {
    coefs_real_final = apply(result_sem$betas_real[(K-20):K,],MARGIN = 2,mean)
    coefs_fake_final = mean(result_sem$betas_fake[(K-20):K])
    
    rmse_ours_1 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 beta_real_gt[beta_real_gt!= 0])^2))
    rmse_ours_2 =sqrt(1/length(mod_unchange$coefficients)*sum((coefs_real_final-
                                                                 mod_real_gt$coefficients)^2))
    se_tmp = mi_se_estimation_glm(M = K_2,mod_real = mod_unchange,
                                  formula_real = formula_real, 
                                  formula_fake =formula_fake,
                                  which_covs_real =which_covs_real,
                                  which_covs_fake =which_covs_fake[1],  
                                  mod_fake = mod_fake_unchange, 
                                  comb_event_stream = comb_event_stream, 
                                  n_actors =n_actors,exo_cov =exo_cov,
                                  exo_cat = exo_cat,seed = 123)
    
    res = list(rmse_naiv_1 = rmse_naiv_1, 
               rmse_naiv_2 = rmse_naiv_2, 
               rmse_ours_1 = rmse_ours_1, 
               rmse_ours_2 = rmse_ours_2, 
               coefs_real_final = coefs_real_final, 
               coefs_fake_final =coefs_fake_final, 
               coefs_real_gt_1 =beta_real_gt[beta_real_gt!= 0], 
               coefs_fake_gt_1 =beta_fake_gt[beta_fake_gt!= 0],
               coefs_real_gt_2 = mod_real_gt$coefficients, 
               coefs_fake_gt_2 =mod_fake_gt$coefficients, 
               coefs_unchange= mod_unchange$coefficients, 
               upper_ci = coefs_real_final + qnorm(p = 1- 0.025)*sqrt(diag(se_tmp$vcov)), 
               lower_ci = coefs_real_final - qnorm(p = 1- 0.025)*sqrt(diag(se_tmp$vcov)), 
               upper_ci_gt = mod_real_gt$coefficients + qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_real_gt))), 
               lower_ci_gt = mod_real_gt$coefficients - qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_real_gt))), 
               upper_ci_unchanged = mod_unchange$coefficients + qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_unchange))), 
               lower_ci_unchanged = mod_unchange$coefficients - qnorm(p = 1- 0.025)*sqrt(diag(vcov(mod_unchange))))
  }
  
  
  return(res)
}



# Function that only does the simulation itself to see what amount of false-positives there are 
simulation_check_number = function(x, n_actors = 40,K = 50,K_2 = 20, exo_cov, exo_cat,
                                   n = 500, beta_real_gt, beta_fake_gt) {
  set.seed(x)
  # exo_cov = rnorm(n = n_actors)
  # exo_cat = sample(size  = n_actors,x = c(1,2,3),replace = T)
  # 
  test_undirected_real = simulation_rem_undirected_new(beta = beta_real_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)
  test_undirected_fake = simulation_rem_undirected_new(beta = beta_fake_gt,n_actors = n_actors,n = n,exo_cov = exo_cov,exo_cat = exo_cat,seed = x)

  end_real = max(test_undirected_real$data$times)
  test_undirected_fake$data = test_undirected_fake$data[times< end_real] 
  
  test_undirected_real$data$ind_real = 1
  test_undirected_fake$data$ind_real = 0
  comb_event_stream = rbind(test_undirected_real$data, test_undirected_fake$data)
  return(table(comb_event_stream$ind_real)/length(comb_event_stream$ind_real))
}


mi_se_estimation = function(M = 50, beta_real_gt,beta_fake_gt,comb_event_stream = comb_event_stream,n_actors = n_actors, 
                              exo_cov,exo_cat,which_covs_real, which_covs_fake, seed,  build = T,build_time = comb_event_stream$times[100]){
  set.seed(seed)
 
  model_list = list()
  for(k in 1:M) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real =beta_real_gt,
                                      beta_fake = beta_fake_gt, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake,
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    
    mod_real =   fastglm(x =   as.matrix(count_real[,which_covs_real, with = F]),
                         y = count_real$status,
                         offset = count_real$offset, family = poisson(),method = 3)
    model_list[[k]] = mod_real
    cat("Iteration ", k," completed\n")
  }
  est_coefs = lapply(model_list,function(x) return(coef(x)))
  est_var = lapply(model_list,function(x) return(x$se))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/length(model_list)*Reduce("+",est_coefs)
  B_bar =1/(length(model_list)-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/length(model_list)*Reduce("+",est_var)
  B = B_bar
  W = V_bar
  final_var_alt = V_bar + B_bar
  return(list(vcov =diag(final_var_alt)))
}


mi_se_estimation_glm = function(M = 50,mod_real,
                                formula_real = formula_real, 
                                formula_fake =formula_fake, 
                                mod_fake, 
                                comb_event_stream = comb_event_stream,n_actors = n_actors, 
                                exo_cov,exo_cat,
                                which_covs_real, 
                                which_covs_fake, 
                                seed,  build = T,
                                build_time = comb_event_stream$times[100]){
  set.seed(seed)
  
  model_list = list()
  for(k in 1:M) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real =mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake[1],
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    model_list[[k]] = mod_real
    cat("Iteration ", k," completed\n")
  }
  est_coefs = lapply(model_list,function(x) return(coef(x)))
  est_var = lapply(model_list,function(x) return(vcov(x)))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/length(model_list)*Reduce("+",est_coefs)
  B_bar =1/(length(model_list)-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/length(model_list)*Reduce("+",est_var)
  B = B_bar
  W = V_bar
  final_var_alt = V_bar + B_bar
  return(list(vcov =final_var_alt))
}



mi_se_estimation_glm = function(M = 50,mod_real,
                                formula_real = formula_real, 
                                formula_fake =formula_fake, 
                                mod_fake, 
                                comb_event_stream = comb_event_stream,n_actors = n_actors, 
                                exo_cov,exo_cat,
                                which_covs_real, 
                                which_covs_fake, 
                                seed,  build = T,
                                build_time = comb_event_stream$times[100]){
  set.seed(seed)
  
  model_list = list()
  for(k in 1:M) {
    # Simulation Step
    comb_event_stream = st_e_step_new(comb_event_stream = comb_event_stream,
                                      n_actors = n_actors, 
                                      beta_real =mod_real$coefficients,
                                      beta_fake = mod_fake$coefficients, 
                                      which_covs_real = which_covs_real, 
                                      which_covs_fake = which_covs_fake[1],
                                      seed = seed+k,
                                      exo_cov = exo_cov,
                                      exo_cat = exo_cat)
    
    count_real = count_rem_undirected_new(event_data = comb_event_stream[ind == 1],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    count_fake = count_rem_undirected_new(event_data = comb_event_stream[ind == 0],
                                          exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat)
    mod_real =   glm(formula_real,family = poisson(), data = count_real)
    model_list[[k]] = mod_real
    cat("Iteration ", k," completed\n")
  }
  est_coefs = lapply(model_list,function(x) return(coef(x)))
  est_var = lapply(model_list,function(x) return(vcov(x)))
  # B_bar = apply(est_coefs,MARGIN = 2,var)
  mean = 1/length(model_list)*Reduce("+",est_coefs)
  B_bar =1/(length(model_list)-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
  V_bar = 1/length(model_list)*Reduce("+",est_var)
  B = B_bar
  W = V_bar
  final_var_alt = V_bar + B_bar
  return(list(vcov =final_var_alt))
}



update_covs_directed_clustered = function(data_present,which_hist,froms, tos, n_actors){
  ind_wheres = findrows(a = froms,b = tos,x = n_actors)
  
  for(i in 1:length(froms)) {
    
    if(data_present$repetition[ind_where[i]] == 0) {
      #Check which(data_present$from == from)
      tmp_out_sender_ind = findfrom(a = from,x = n_actors)
      tmp_in_sender_ind = findfrom(a = to,x = n_actors)
      tmp_out_receiver_ind =  findto(b = from,x = n_actors)
      tmp_in_receiver_ind =findto(b = to,x = n_actors)
      
      data_present$out_sender[tmp_out_sender_ind] = data_present$out_sender[tmp_out_sender_ind] + 1
      
      # Outdegree Receiver
      # Increment 1 the outdegree receiver for all events that may take place with i being the receiver
      
      #Check which(data_present$to == from)
      data_present$out_receiver[tmp_out_receiver_ind] = data_present$out_receiver[tmp_out_receiver_ind] +1
      
      # Indegree Sender
      # Increment 1 the indegree sender for all events that may take place with j being the sender
      #Check which(data_present$from == to)
      data_present$in_sender[tmp_in_sender_ind] = data_present$in_sender[tmp_in_sender_ind] +1
      
      # Indegree Receiver
      # Increment 1 the Indegree receiver for all events that may take place with j being the receiver
      #Check which(data_present$to == to)
      data_present$in_receiver[tmp_in_receiver_ind] = data_present$in_receiver[tmp_in_receiver_ind] +1
      
      # Transitivity
      # Increment 1 the transitivity if h was in sender j and may be in sender i 
      # or h was in receiver i and will be in receiver j 
      
      # What actors did to already send to? 
      # which((data_present$from == to) & (data_present$history == TRUE))
      ind_sender_to_j = intersect(tmp_in_sender_ind, which_hist)
      
      # Whats the position of from to those h 
      tmp = data_present$to[ind_sender_to_j]
      tmp = tmp[tmp!= from]
      ind_transitivity_1 = findrows_1(bs = tmp,a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$to[which((data_present$from == to) & (data_present$history == TRUE))] )
      
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(tmp_out_receiver_ind, which_hist)
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_transitivity_2 = findrows_2(as =  tmp,b = to, x = n_actors)
      
      #which(data_present$to == to & data_present$from %in% data_present$from[which((data_present$to == from) & (data_present$history == TRUE))])
      ind_transitivity = unique(c(ind_transitivity_1,ind_transitivity_2))
      ind_transitivity = ind_transitivity[ind_transitivity<= n_actors * (n_actors-1)]
      
      data_present$transitivity[ind_transitivity] = data_present$transitivity[ind_transitivity] +1
      
      
      # Shared contactees
      # Increment 1 shared contactees if h was in receiver j in the past and h may be sender future of i 
      # or if h was in receiver of j and will be in receiver of i
      
      #What possible ties did already send to j in the past 
      ind_sender_to_j = intersect(tmp_in_receiver_ind, which_hist)
      ind_sender_to_j = ind_sender_to_j[ind_sender_to_j!= ind_where]
      # which((data_present$to == to) & (data_present$history == TRUE))
      # Which ties may i send to of these 
      
      ind_shared_contactees_1 = findrows_1(bs =  data_present$from[ind_sender_to_j],
                                           a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$from[which((data_present$to == to) & (data_present$history == TRUE))])
      # Which of these may send to i 
      
      ind_shared_contactees_2 = findrows_2(as =  data_present$from[ind_sender_to_j],
                                           b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$from[ind_sender_to_j])
      
      ind_shared_contactees = unique(c(ind_shared_contactees_1,ind_shared_contactees_2))
      ind_shared_contactees = ind_shared_contactees[ind_shared_contactees<= n_actors * (n_actors-1)]
      data_present$shared_contactees[ind_shared_contactees] = data_present$shared_contactees[ind_shared_contactees] +1
      # Triangle Closure 
      # Increment 1 the triangle closure if h was in receiver i and may be in sender j 
      # or h was in sender j and will be in receiver i 
      #What possible ties were already sent from j 
      ind_sender_from_j =   intersect(which_hist,tmp_in_sender_ind)
      ind_sender_from_j = ind_sender_from_j[ind_sender_from_j!= ind_where]
      
      # which((data_present$from == to) & (data_present$history == TRUE))
      # Which of these ties have i as the receiver? 
      tmp = data_present$to[ind_sender_from_j]
      tmp = tmp[tmp!= from]
      ind_triangle_1 =  findrows_2(as = tmp,
                                   b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$to[ind_sender_from_j])
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(which_hist,tmp_out_receiver_ind)
      ind_sender_to_i = ind_sender_to_i[ind_sender_to_i!= ind_where]
      
      # which((data_present$to == from) & (data_present$history == TRUE))
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_triangle_2 = findrows_1(bs =  tmp,
                                  a = to, x = n_actors)
      # which(data_present$from == to & data_present$to %in% data_present$from[ind_sender_to_i])
      ind_triangle = unique(c(ind_triangle_1,ind_triangle_2))
      ind_triangle = ind_triangle[ind_triangle<= n_actors * (n_actors-1)]
      
      data_present$triangle[ind_triangle] = data_present$triangle[ind_triangle] +1
      
      # Shared contacters
      # Increment 1 shared contacters if h was in sender i and will be in sender j 
      # or h was in sender i and will be in receiver j 
      #What possible ties were already sent from i 
      ind_sender_from_i = intersect(which_hist,tmp_out_sender_ind)
      ind_sender_from_i = ind_sender_from_i[ind_sender_from_i!= ind_where]
      
      # which((data_present$from == from) & (data_present$history == TRUE))
      # Which ties may i send to
      
      
      ind_shared_contacters_1 = findrows_1(bs =  data_present$to[ind_sender_from_i],
                                           a = to, x = n_actors)
      
      # which(data_present$from == to & data_present$to %in% data_present$to[ind_sender_from_i])
      
      ind_shared_contacters_2 = findrows_2(as = data_present$to[ind_sender_from_i],
                                           b = to, x = n_actors)
      
      # which(data_present$to == to & data_present$from %in% data_present$to[ind_sender_from_i])
      ind_shared_contacters = unique(c(ind_shared_contacters_1,ind_shared_contacters_2))
      ind_shared_contacters = ind_shared_contacters[ind_shared_contacters<= n_actors * (n_actors-1)]
      data_present$shared_contacters[ind_shared_contacters] = data_present$shared_contacters[ind_shared_contacters] +1 
      # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
     
      
    } else {
      # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
      
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
    }
    
    which_hist <- c(which_hist,ind_wheres)
    
  }
  return(list(data_present = data_present, which_hist = which_hist))
  
  # return(data_present)
}
  
update_covs_directed = function(data_present,which_hist,from, to, n_actors){
  ind_where = findrow(a = from,b = to,x = n_actors)
  
  if(data_present$repetition[ind_where] == 0) {
    #Check which(data_present$from == from)
    tmp_out_sender_ind = findfrom(a = from,x = n_actors)
    tmp_in_sender_ind = findfrom(a = to,x = n_actors)
    tmp_out_receiver_ind =  findto(b = from,x = n_actors)
    tmp_in_receiver_ind =findto(b = to,x = n_actors)
    
    data_present$out_sender[tmp_out_sender_ind] = data_present$out_sender[tmp_out_sender_ind] + 1
    
    # Outdegree Receiver
    # Increment 1 the outdegree receiver for all events that may take place with i being the receiver
    
    #Check which(data_present$to == from)
    data_present$out_receiver[tmp_out_receiver_ind] = data_present$out_receiver[tmp_out_receiver_ind] +1
    
    # Indegree Sender
    # Increment 1 the indegree sender for all events that may take place with j being the sender
    #Check which(data_present$from == to)
    data_present$in_sender[tmp_in_sender_ind] = data_present$in_sender[tmp_in_sender_ind] +1
    
    # Indegree Receiver
    # Increment 1 the Indegree receiver for all events that may take place with j being the receiver
    #Check which(data_present$to == to)
    data_present$in_receiver[tmp_in_receiver_ind] = data_present$in_receiver[tmp_in_receiver_ind] +1
    
    # Transitivity
    # Increment 1 the transitivity if h was in sender j and may be in sender i 
    # or h was in receiver i and will be in receiver j 
    
    # What actors did to already send to? 
    # which((data_present$from == to) & (data_present$history == TRUE))
    ind_sender_to_j = intersect(tmp_in_sender_ind, which_hist)
    
    # Whats the position of from to those h 
    tmp = data_present$to[ind_sender_to_j]
    tmp = tmp[tmp!= from]
    ind_transitivity_1 = findrows_1(bs = tmp,a = from, x = n_actors)
    
    # which(data_present$from == from & data_present$to %in% data_present$to[which((data_present$from == to) & (data_present$history == TRUE))] )
    
    # What possible ties did already send to i in the past
    ind_sender_to_i= intersect(tmp_out_receiver_ind, which_hist)
    tmp = data_present$from[ind_sender_to_i]
    tmp = tmp[tmp!= to]
    
    ind_transitivity_2 = findrows_2(as =  tmp,b = to, x = n_actors)
    
    #which(data_present$to == to & data_present$from %in% data_present$from[which((data_present$to == from) & (data_present$history == TRUE))])
    ind_transitivity = unique(c(ind_transitivity_1,ind_transitivity_2))
    ind_transitivity = ind_transitivity[ind_transitivity<= n_actors * (n_actors-1)]
    
    data_present$transitivity[ind_transitivity] = data_present$transitivity[ind_transitivity] +1
    
    
    # Shared contactees
    # Increment 1 shared contactees if h was in receiver j in the past and h may be sender future of i 
    # or if h was in receiver of j and will be in receiver of i
    
    #What possible ties did already send to j in the past 
    ind_sender_to_j = intersect(tmp_in_receiver_ind, which_hist)
    ind_sender_to_j = ind_sender_to_j[ind_sender_to_j!= ind_where]
    # which((data_present$to == to) & (data_present$history == TRUE))
    # Which ties may i send to of these 
    
    ind_shared_contactees_1 = findrows_1(bs =  data_present$from[ind_sender_to_j],
                                         a = from, x = n_actors)
    
    # which(data_present$from == from & data_present$to %in% data_present$from[which((data_present$to == to) & (data_present$history == TRUE))])
    # Which of these may send to i 
    
    ind_shared_contactees_2 = findrows_2(as =  data_present$from[ind_sender_to_j],
                                         b = from, x = n_actors)
    
    # which(data_present$to == from & data_present$from %in% data_present$from[ind_sender_to_j])
    
    ind_shared_contactees = unique(c(ind_shared_contactees_1,ind_shared_contactees_2))
    ind_shared_contactees = ind_shared_contactees[ind_shared_contactees<= n_actors * (n_actors-1)]
    data_present$shared_contactees[ind_shared_contactees] = data_present$shared_contactees[ind_shared_contactees] +1
    # Triangle Closure 
    # Increment 1 the triangle closure if h was in receiver i and may be in sender j 
    # or h was in sender j and will be in receiver i 
    #What possible ties were already sent from j 
    ind_sender_from_j =   intersect(which_hist,tmp_in_sender_ind)
    ind_sender_from_j = ind_sender_from_j[ind_sender_from_j!= ind_where]
    
    # which((data_present$from == to) & (data_present$history == TRUE))
    # Which of these ties have i as the receiver? 
    tmp = data_present$to[ind_sender_from_j]
    tmp = tmp[tmp!= from]
    ind_triangle_1 =  findrows_2(as = tmp,
                                 b = from, x = n_actors)
    
    # which(data_present$to == from & data_present$from %in% data_present$to[ind_sender_from_j])
    # What possible ties did already send to i in the past
    ind_sender_to_i= intersect(which_hist,tmp_out_receiver_ind)
    ind_sender_to_i = ind_sender_to_i[ind_sender_to_i!= ind_where]
    
    # which((data_present$to == from) & (data_present$history == TRUE))
    tmp = data_present$from[ind_sender_to_i]
    tmp = tmp[tmp!= to]
    
    ind_triangle_2 = findrows_1(bs =  tmp,
                                a = to, x = n_actors)
    # which(data_present$from == to & data_present$to %in% data_present$from[ind_sender_to_i])
    ind_triangle = unique(c(ind_triangle_1,ind_triangle_2))
    ind_triangle = ind_triangle[ind_triangle<= n_actors * (n_actors-1)]
    
    data_present$triangle[ind_triangle] = data_present$triangle[ind_triangle] +1
    
    # Shared contacters
    # Increment 1 shared contacters if h was in sender i and will be in sender j 
    # or h was in sender i and will be in receiver j 
    #What possible ties were already sent from i 
    ind_sender_from_i = intersect(which_hist,tmp_out_sender_ind)
    ind_sender_from_i = ind_sender_from_i[ind_sender_from_i!= ind_where]
    
    # which((data_present$from == from) & (data_present$history == TRUE))
    # Which ties may i send to
    
    
    ind_shared_contacters_1 = findrows_1(bs =  data_present$to[ind_sender_from_i],
                                         a = to, x = n_actors)
    
    # which(data_present$from == to & data_present$to %in% data_present$to[ind_sender_from_i])
    
    ind_shared_contacters_2 = findrows_2(as = data_present$to[ind_sender_from_i],
                                         b = to, x = n_actors)
    
    # which(data_present$to == to & data_present$from %in% data_present$to[ind_sender_from_i])
    ind_shared_contacters = unique(c(ind_shared_contacters_1,ind_shared_contacters_2))
    ind_shared_contacters = ind_shared_contacters[ind_shared_contacters<= n_actors * (n_actors-1)]
    data_present$shared_contacters[ind_shared_contacters] = data_present$shared_contacters[ind_shared_contacters] +1 
    # Reciprocity
    # Increment 1 the reciprocity for the inverse receiver/sender tie
    tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
    #Check which((data_present$to == from) &  (data_present$from == to))
    data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    which_hist <- c(which_hist,ind_where)
    
  } else {
    # Reciprocity
    # Increment 1 the reciprocity for the inverse receiver/sender tie
    tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
    #Check which((data_present$to == from) &  (data_present$from == to))
    data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
    
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
  }
  return(list(data_present = data_present, which_hist = which_hist))
  
  # return(data_present)
}


simulation_rem_new = function(beta, n_actors, n,exo_cov){
  time_beg = Sys.time()
  
  pb <- txtProgressBar(min = 0, max =n, style = 3)
  obs = data.table(from =numeric(length = n),from =numeric(length = n))
  # What are the possibles combinations one could observe? 
  possibe_comb = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  # possibe_comb_comp = matrix(combinations(x = n_actors),ncol = 2,byrow = T)
  # names(possibe_comb) = c("from","to")
  # # Which events were already observed?
  # possibe_comb$history = FALSE
  n_possible = nrow(possibe_comb)
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("from","to")
  data_present$intercept= 1
  data_present$out_sender= numeric(length = nrow(possibe_comb))
  data_present$in_sender= numeric(length = nrow(possibe_comb))
  data_present$out_receiver= numeric(length = nrow(possibe_comb))
  data_present$in_receiver= numeric(length = nrow(possibe_comb))
  data_present$reciprocity= numeric(length = nrow(possibe_comb))
  data_present$repetition= numeric(length = nrow(possibe_comb))
  data_present$transitivity= numeric(length = nrow(possibe_comb))
  data_present$shared_contactees= numeric(length = nrow(possibe_comb))
  data_present$triangle= numeric(length = nrow(possibe_comb))
  data_present$shared_contacters= numeric(length = nrow(possibe_comb))
  
  data_present$exo_from = exo_cov[data_present$from]
  data_present$exo_to = exo_cov[data_present$to]
  
  lambda =  numeric(length = nrow(possibe_comb))
  event_data = data.table("from" = numeric(length = n), "to"= numeric(length = n), "times" = numeric(length = n))
  which_hist = c()
  for(i in 1:n){
    setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda  = as.vector(exp(beta %*% t(data_present[,-c(1,2)])))
    sum_lambda = sum(lambda)
    prob_lambda = lambda/sum_lambda
    waiting_time = rexp(n = 1, rate = sum_lambda)
    event = sample(x = n_possible,size = 1,prob = prob_lambda)
    # Update the covariates in data_present accordingly 
    from = data_present$from[event]
    to = data_present$to[event]
    event_data[i,"from"] = from
    event_data[i,"to"] = to 
    event_data[i,"times"] = waiting_time
    
    #Step b: Indicate that the event has happenend at least once
    tmp =  update_covs_directed(data_present = data_present, 
                                which_hist = which_hist,from = from, 
                                to = to, n_actors =n_actors) 
    data_present = tmp$data_present
    which_hist =tmp$which_hist
 
  }
  event_data$times = cumsum(event_data$times)
  close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(list("data"= event_data, "end_stage" = data_present))
}


simulation_rem = function(beta, n_actors, n,exo_cov){
  time_beg = Sys.time()
  
  pb <- txtProgressBar(min = 0, max =n, style = 3)
  obs = data.table(from =numeric(length = n),from =numeric(length = n))
  # What are the possibles combinations one could observe? 
  possibe_comb = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  # possibe_comb_comp = matrix(combinations(x = n_actors),ncol = 2,byrow = T)
  # names(possibe_comb) = c("from","to")
  # # Which events were already observed?
  # possibe_comb$history = FALSE
  n_possible = nrow(possibe_comb)
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("from","to")
  data_present$intercept= 1
  data_present$out_sender= numeric(length = nrow(possibe_comb))
  data_present$in_sender= numeric(length = nrow(possibe_comb))
  data_present$out_receiver= numeric(length = nrow(possibe_comb))
  data_present$in_receiver= numeric(length = nrow(possibe_comb))
  data_present$reciprocity= numeric(length = nrow(possibe_comb))
  data_present$repetition= numeric(length = nrow(possibe_comb))
  data_present$transitivity= numeric(length = nrow(possibe_comb))
  data_present$shared_contactees= numeric(length = nrow(possibe_comb))
  data_present$triangle= numeric(length = nrow(possibe_comb))
  data_present$shared_contacters= numeric(length = nrow(possibe_comb))
  
  data_present$exo_from = exo_cov[data_present$from]
  data_present$exo_to = exo_cov[data_present$to]
  
  lambda =  numeric(length = nrow(possibe_comb))
  event_data = data.table("from" = numeric(length = n), "to"= numeric(length = n), "times" = numeric(length = n))
  which_hist = c()
  for(i in 1:n){
    setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda  = as.vector(exp(beta %*% t(data_present[,-c(1,2)])))
    sum_lambda = sum(lambda)
    prob_lambda = lambda/sum_lambda
    waiting_time = rexp(n = 1, rate = sum_lambda)
    event = sample(x = n_possible,size = 1,prob = prob_lambda)
    # Update the covariates in data_present accordingly 
    from = data_present$from[event]
    to = data_present$to[event]
    event_data[i,"from"] = from
    event_data[i,"to"] = to 
    event_data[i,"times"] = waiting_time
    
    #Step b: Indicate that the event has happenend at least once
    
    ind_where = findrow(a = from,b = to,x = n_actors)
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    if(data_present$repetition[ind_where] == 0) {
      #Check which(data_present$from == from)
      tmp_out_sender_ind = findfrom(a = from,x = n_actors)
      tmp_in_sender_ind = findfrom(a = to,x = n_actors)
      tmp_out_receiver_ind =  findto(b = from,x = n_actors)
      tmp_in_receiver_ind =findto(b = to,x = n_actors)
      
      data_present$out_sender[tmp_out_sender_ind] = data_present$out_sender[tmp_out_sender_ind] + 1
      
      # Outdegree Receiver
      # Increment 1 the outdegree receiver for all events that may take place with i being the receiver
      
      #Check which(data_present$to == from)
      data_present$out_receiver[tmp_out_receiver_ind] = data_present$out_receiver[tmp_out_receiver_ind] +1
      
      # Indegree Sender
      # Increment 1 the indegree sender for all events that may take place with j being the sender
      #Check which(data_present$from == to)
      data_present$in_sender[tmp_in_sender_ind] = data_present$in_sender[tmp_in_sender_ind] +1
      
      # Indegree Receiver
      # Increment 1 the Indegree receiver for all events that may take place with j being the receiver
      #Check which(data_present$to == to)
      data_present$in_receiver[tmp_in_receiver_ind] = data_present$in_receiver[tmp_in_receiver_ind] +1
      
      # Transitivity
      # Increment 1 the transitivity if h was in sender j and may be in sender i 
      # or h was in receiver i and will be in receiver j 
      
      # What actors did to already send to? 
      # which((data_present$from == to) & (data_present$history == TRUE))
      ind_sender_to_j = intersect(tmp_in_sender_ind, which_hist)
      
      # Whats the position of from to those h 
      tmp = data_present$to[ind_sender_to_j]
      tmp = tmp[tmp!= from]
      ind_transitivity_1 = findrows_1(bs = tmp,a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$to[which((data_present$from == to) & (data_present$history == TRUE))] )
      
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(tmp_out_receiver_ind, which_hist)
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_transitivity_2 = findrows_2(as =  tmp,b = to, x = n_actors)
      
      #which(data_present$to == to & data_present$from %in% data_present$from[which((data_present$to == from) & (data_present$history == TRUE))])
      ind_transitivity = unique(c(ind_transitivity_1,ind_transitivity_2))
      ind_transitivity = ind_transitivity[ind_transitivity<= n_actors * (n_actors-1)]
      
      data_present$transitivity[ind_transitivity] = data_present$transitivity[ind_transitivity] +1
      
      
      # Shared contactees
      # Increment 1 shared contactees if h was in receiver j in the past and h may be sender future of i 
      # or if h was in receiver of j and will be in receiver of i
      
      #What possible ties did already send to j in the past 
      ind_sender_to_j = intersect(tmp_in_receiver_ind, which_hist)
      ind_sender_to_j = ind_sender_to_j[ind_sender_to_j!= ind_where]
      # which((data_present$to == to) & (data_present$history == TRUE))
      # Which ties may i send to of these 
      
      ind_shared_contactees_1 = findrows_1(bs =  data_present$from[ind_sender_to_j],
                                           a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$from[which((data_present$to == to) & (data_present$history == TRUE))])
      # Which of these may send to i 
      
      ind_shared_contactees_2 = findrows_2(as =  data_present$from[ind_sender_to_j],
                                           b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$from[ind_sender_to_j])
      
      ind_shared_contactees = unique(c(ind_shared_contactees_1,ind_shared_contactees_2))
      ind_shared_contactees = ind_shared_contactees[ind_shared_contactees<= n_actors * (n_actors-1)]
      data_present$shared_contactees[ind_shared_contactees] = data_present$shared_contactees[ind_shared_contactees] +1
      # Triangle Closure 
      # Increment 1 the triangle closure if h was in receiver i and may be in sender j 
      # or h was in sender j and will be in receiver i 
      #What possible ties were already sent from j 
      ind_sender_from_j =   intersect(which_hist,tmp_in_sender_ind)
      ind_sender_from_j = ind_sender_from_j[ind_sender_from_j!= ind_where]
      
      # which((data_present$from == to) & (data_present$history == TRUE))
      # Which of these ties have i as the receiver? 
      tmp = data_present$to[ind_sender_from_j]
      tmp = tmp[tmp!= from]
      ind_triangle_1 =  findrows_2(as = tmp,
                                   b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$to[ind_sender_from_j])
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(which_hist,tmp_out_receiver_ind)
      ind_sender_to_i = ind_sender_to_i[ind_sender_to_i!= ind_where]
      
      # which((data_present$to == from) & (data_present$history == TRUE))
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_triangle_2 = findrows_1(bs =  tmp,
                                  a = to, x = n_actors)
      # which(data_present$from == to & data_present$to %in% data_present$from[ind_sender_to_i])
      ind_triangle = unique(c(ind_triangle_1,ind_triangle_2))
      ind_triangle = ind_triangle[ind_triangle<= n_actors * (n_actors-1)]
      
      data_present$triangle[ind_triangle] = data_present$triangle[ind_triangle] +1
      
      # Shared contacters
      # Increment 1 shared contacters if h was in sender i and will be in sender j 
      # or h was in sender i and will be in receiver j 
      #What possible ties were already sent from i 
      ind_sender_from_i = intersect(which_hist,tmp_out_sender_ind)
      ind_sender_from_i = ind_sender_from_i[ind_sender_from_i!= ind_where]
      
      # which((data_present$from == from) & (data_present$history == TRUE))
      # Which ties may i send to
      
      
      ind_shared_contacters_1 = findrows_1(bs =  data_present$to[ind_sender_from_i],
                                           a = to, x = n_actors)
      
      # which(data_present$from == to & data_present$to %in% data_present$to[ind_sender_from_i])
      
      ind_shared_contacters_2 = findrows_2(as = data_present$to[ind_sender_from_i],
                                           b = to, x = n_actors)
      
      # which(data_present$to == to & data_present$from %in% data_present$to[ind_sender_from_i])
      ind_shared_contacters = unique(c(ind_shared_contacters_1,ind_shared_contacters_2))
      ind_shared_contacters = ind_shared_contacters[ind_shared_contacters<= n_actors * (n_actors-1)]
      data_present$shared_contacters[ind_shared_contacters] = data_present$shared_contacters[ind_shared_contacters] +1 
      # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      which_hist = c(which_hist,ind_where)
      
    } else {
      # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
      
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    }
  }
  event_data$times = cumsum(event_data$times)
  close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(list("data"= event_data, "end_stage" = data_present))
}

findrow_undirected = function(a,b,x){
  return(c(0,cumsum(x-1:( max(a)-1)))[a] + (b-a))
  # return(sum(x-1:(a-1)) + (b-a))
}
find_side_a = function(a,x){
 
  if(a != x){
    # return(sum(x-1:(a-1)) : (sum(x-1:(a))-1 ) +1)
    return(c(0,cumsum(x-1:( a-1)))[a] + ((a+1):x-a))
  }
}
find_side_b = function(b,x){
  # return(c(0,cumsum(x-1:(x-2)))[1:(b-1)] + (b-1:(b-1)))
  if(b != 1){
    return(c(1,cumsum(x-1:(x-1)) +1)[1:(b-1)]  + (b-1):1- 1)
  }
}

findrows_undirected = function(as, bs, x) {
  return(c(0,cumsum(x-1:( max(as)-1)))[as] + (bs-as))
}

st_e_step_both = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist_real = c()
  which_hist_fake = c()
  data_present_real = data_present
  data_present_fake = data_present
  for(i in 1:n){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present_real[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present_fake[,matched_fake, with = F])))
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    
    # Sample the event (fake or true)
    prob = lambda_real[ind_where]/
      (lambda_real[ind_where] + lambda_fake[ind_where])
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      tmp =  update_covs_undirected_both(data_present =data_present_real,which_hist = which_hist_real,
                                         side_a =side_a, side_b = side_b, n_actors =n_actors,real =T) 
      data_present_real = tmp$data_present
      which_hist_real =tmp$which_hist
    } else {
      tmp = update_covs_undirected_both(data_present =data_present_fake,which_hist = which_hist_fake,
                                                      side_a =side_a, side_b = side_b, n_actors =n_actors, real = F) 
      data_present_fake = tmp$data_present
      which_hist_fake =tmp$which_hist
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 


st_e_step_old = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  for(i in 1:n){

    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    
    # Sample the event (fake or true)
    prob = lambda_real[ind_where]/
      (lambda_real[ind_where] + lambda_fake[ind_where])
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      tmp = update_covs_undirected(data_present =data_present,which_hist = which_hist,side_a =side_a, side_b = side_b, n_actors =n_actors) 
      which_hist = tmp$which_hist
      data_present = tmp$data_present
      }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 
score_function = function(tmp_ind, data, model, which_covs_real) {
  # the function returns x*y - exp(theta^Tx + log(offset))*x
  data = data[ind == tmp_ind]
  return(colSums(as.matrix(data[,which_covs_real, with = F]) *data$status - 
                   as.vector(exp(as.matrix(data[,which_covs_real, with = F]) %*% model$coefficients +data$offset))*
                   as.matrix(data[,which_covs_real, with = F])))
}




hessian_function_new = function(tmp_ind,data,model, which_covs_real, target, offset){
  data_tmp = data[ind == tmp_ind,which_covs_real, with = F]
  data.list <- split(data_tmp, seq(nrow(data_tmp)))
  data.list = lapply(data.list, function(x){t(matrix(x))})
  trying = bdiag(data.list)
  res = t(data_tmp)%*%trying
  p = length(which_covs_real)
  n = nrow(data_tmp)
  res_array = array(data = res, dim = c(p,p,n))
  res_array = aperm(res_array, c(3,2,1)) 
  
  res_array= res_array[,,] * - as.vector(exp(as.matrix(data[ind == tmp_ind,which_covs_real, with = F]) %*% model$coefficients +data[ind == tmp_ind]$offset))
  
  return(colSums(res_array, dims = 1))
}


hessian_function = function(tmp_ind,data,model, which_covs_real, target, offset){
  data = data[ind == tmp_ind]
  
  y = aaply(1:nrow(data),.margins = 1,
            .fun =  function(x){(as.matrix(data[,which_covs_real, with = F])[x,])%*% t(as.matrix(data[,which_covs_real, with = F])[x,])})
  y = y[,,] * - as.vector(exp(as.matrix(data[,which_covs_real, with = F]) %*% model$coefficients +data$offset))
  return(colSums(y, dims = 1))
}

st_e_step_new = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  comb_event_stream$between_time = c(comb_event_stream$times[1],diff(comb_event_stream$times))
  for(i in 1:n){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    
    # lambda_real  = as.vector(comb_event_stream$between_time[i]*exp(beta_real %*% t(data_present[,matched_real, with = F])))
    # lambda_fake  = as.vector(comb_event_stream$between_time[i]*exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    # 
    
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    
    # Sample the event (fake or true)
    prob = lambda_real[ind_where]/
      (lambda_real[ind_where] + lambda_fake[ind_where])
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      tmp = update_covs_undirected(data_present =data_present,which_hist = which_hist,side_a =side_a, side_b = side_b, n_actors =n_actors) 
      which_hist = tmp$which_hist
      data_present = tmp$data_present
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 

st_e_step_alt = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, pi, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  comb_event_stream$between_time = c(comb_event_stream$times[1],diff(comb_event_stream$times))
  last_event_time = 0
  for(i in 1:n){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 

    # lambda_real  = as.vector(comb_event_stream$between_time[i]*exp(beta_real %*% t(data_present[,matched_real, with = F])))
    # lambda_fake  = as.vector(comb_event_stream$between_time[i]*exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    # 
    
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    time_tmp = comb_event_stream$times[i]
    
    offset = log(time_tmp- last_event_time)
    last_event_time = time_tmp
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F]) +offset))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])+offset))
    
    
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    
    # Sample the event (fake or true)
    # prob = lambda_real[ind_where]/
    #   (lambda_real[ind_where] + lambda_fake[ind_where])
    prob = (pi*dpois(x = 1, lambda = lambda_real[ind_where]))/
      (pi*dpois(x = 1, lambda = lambda_real[ind_where]) + (1-pi)*dpois(x = 1, lambda = lambda_fake[ind_where]))
    
    
    
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      tmp = update_covs_undirected(data_present =data_present,which_hist = which_hist,side_a =side_a, side_b = side_b, n_actors =n_actors) 
      which_hist = tmp$which_hist
      data_present = tmp$data_present
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 


st_e_step_list = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  
  if(number_cont!=0){
    names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  comb_event_stream$between_time = c(comb_event_stream$times[1],diff(comb_event_stream$times))
  for(i in 1:n){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    
    # lambda_real  = as.vector(comb_event_stream$between_time[i]*exp(beta_real %*% t(data_present[,matched_real, with = F])))
    # lambda_fake  = as.vector(comb_event_stream$between_time[i]*exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    # 
    
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    
    # Sample the event (fake or true)
    prob = lambda_real[ind_where]/
      (lambda_real[ind_where] + lambda_fake[ind_where])
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      tmp = update_covs_undirected(data_present =data_present,which_hist = which_hist,side_a =side_a, side_b = side_b, n_actors =n_actors) 
      which_hist = tmp$which_hist
      data_present = tmp$data_present
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 


st_e_step = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_cov = exo_cov,exo_cat = exo_cat) {
  time_beg = Sys.time()
  set.seed(seed)
  n = nrow(comb_event_stream)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  for(i in 1:n){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    # cat(i,"\n")
    # cat(mean(lambda_1),"\n")
    # cat(mean(lambda_2),"\n")
    side_a = comb_event_stream$side_a[i]
    side_b = comb_event_stream$side_b[i]
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    # Sample the event (fake or true)
    prob = lambda_real[ind_where]/
      (lambda_real[ind_where] + lambda_fake[ind_where])
    # ind_true = rbinom(n = 1, size = 1, 
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where]),
    #                             lambda_2[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    # ind_true = rbinom(n = 1, size = 1,
    #                   prob =  c(lambda_1[ind_where]/(lambda_1[ind_where] + lambda_2[ind_where])))
    
    comb_event_stream$ind[i] = rbinom(n = 1, size = 1,prob = prob)
    if(comb_event_stream$ind[i]){
      if(data_present$repetition[ind_where] == 0) {
        #Check which(data_present$from == from)
        tmp_deg_side_a_1 = find_side_a(a = side_a,x = n_actors)
        tmp_deg_side_a_2 =find_side_b(b = side_a,x = n_actors)
        tmp_deg_side_b_1 = find_side_a(a = side_b,x = n_actors)
        tmp_deg_side_b_2 =find_side_b(b = side_b,x = n_actors)
        
        data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
        data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
        data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
        data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
        
        # Transitivity
        # Increment 1 the transitivity if h was in sender j and may be in sender i 
        # or h was in receiver i and will be in receiver j 
        
        # What events already happened with side_a?
        ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist)
        ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist)
        # with whom were those events? 
        triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
        
        # What events already happened with side_b?
        ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist)
        ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist)
        # with whom were those events? 
        triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
        
        if(length(triangle_tmp_1)!=0) {
          side_a_tmp = triangle_tmp_1
          side_b_tmp = rep(side_b, length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        if(length(triangle_tmp_2)!=0) {
          side_a_tmp = triangle_tmp_2
          side_b_tmp = rep(side_a, length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        
        data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
        which_hist = c(which_hist,ind_where)
        data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
        data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
      } else {
        
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      }
      data_present$repetition = data_present$repetition>0
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(comb_event_stream)
} 


update_covs_undirected_both = function(data_present,which_hist,side_a, side_b, n_actors, change = F, real){
  ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
  #Step c: Search all edges that will occurr in the future and will also be affected 
  # Assuming that the actual event is i-> j 
  # 
  # Where will changes occurr? 
  # Outdegree Sender
  
  if(data_present$repetition[ind_where] == 0) {
    #Check which(data_present$from == from)
    tmp_deg_side_a_1 = find_side_a(a = side_a,x = n_actors)
    tmp_deg_side_a_2 =find_side_b(b = side_a,x = n_actors)
    tmp_deg_side_b_1 = find_side_a(a = side_b,x = n_actors)
    tmp_deg_side_b_2 =find_side_b(b = side_b,x = n_actors)
    
    data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
    data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
    data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
    data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
    
    # Transitivity
    # Increment 1 the transitivity if h was in sender j and may be in sender i 
    # or h was in receiver i and will be in receiver j 
    
    # What events already happened with side_a?
    ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist)
    ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist)
    # with whom were those events? 
    triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
    
    # What events already happened with side_b?
    ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist)
    ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist)
    # with whom were those events? 
    triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
    
    if(length(triangle_tmp_1)!=0) {
      side_a_tmp = triangle_tmp_1
      side_b_tmp = rep(side_b, length(side_a_tmp))
      ind_change = side_a_tmp>side_b_tmp
      tmp = side_a_tmp[ind_change]
      side_a_tmp[ind_change] = side_b_tmp[ind_change]
      side_b_tmp[ind_change] = tmp
      triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
    }
    
    if(length(triangle_tmp_2)!=0) {
      side_a_tmp = triangle_tmp_2
      side_b_tmp = rep(side_a, length(side_a_tmp))
      ind_change = side_a_tmp>side_b_tmp
      tmp = side_a_tmp[ind_change]
      side_a_tmp[ind_change] = side_b_tmp[ind_change]
      side_b_tmp[ind_change] = tmp
      triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
    }
    
    data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    which_hist <- c(which_hist,ind_where)
  
    data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
    data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
    
    if(change) {
      change_ind = c(tmp_deg_side_a_1,
                     tmp_deg_side_a_2,
                     tmp_deg_side_b_1,
                     tmp_deg_side_b_2, 
                     triangle_tmp_1, 
                     triangle_tmp_2, 
                     ind_where)
      data_present$change[change_ind] = 1
    }
    
  } else {
    
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    if(change) {
      data_present$change[ind_where] = 1
    }
    
  }
  

  return(list(data_present = data_present, which_hist = which_hist))
  
  # return(data_present)
}


update_covs_undirected = function(data_present,which_hist,side_a, side_b, n_actors, change = F){
  ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
  #Step c: Search all edges that will occurr in the future and will also be affected 
  # Assuming that the actual event is i-> j 
  # 
  # Where will changes occurr? 
  # Outdegree Sender
  
  if(data_present$repetition[ind_where] == 0) {
    #Check which(data_present$from == from)
    tmp_deg_side_a_1 = find_side_a(a = side_a,x = n_actors)
    tmp_deg_side_a_2 =find_side_b(b = side_a,x = n_actors)
    tmp_deg_side_b_1 = find_side_a(a = side_b,x = n_actors)
    tmp_deg_side_b_2 =find_side_b(b = side_b,x = n_actors)
    
    data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
    data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
    data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
    data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
    
    # Transitivity
    # Increment 1 the transitivity if h was in sender j and may be in sender i 
    # or h was in receiver i and will be in receiver j 
    
    # What events already happened with side_a?
    ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist)
    ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist)
    # with whom were those events? 
    triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
    
    # What events already happened with side_b?
    ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist)
    ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist)
    # with whom were those events? 
    triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
    
    if(length(triangle_tmp_1)!=0) {
      side_a_tmp = triangle_tmp_1
      side_b_tmp = rep(side_b, length(side_a_tmp))
      ind_change = side_a_tmp>side_b_tmp
      tmp = side_a_tmp[ind_change]
      side_a_tmp[ind_change] = side_b_tmp[ind_change]
      side_b_tmp[ind_change] = tmp
      triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
    }
    
    if(length(triangle_tmp_2)!=0) {
      side_a_tmp = triangle_tmp_2
      side_b_tmp = rep(side_a, length(side_a_tmp))
      ind_change = side_a_tmp>side_b_tmp
      tmp = side_a_tmp[ind_change]
      side_a_tmp[ind_change] = side_b_tmp[ind_change]
      side_b_tmp[ind_change] = tmp
      triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
    }
    
    data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    which_hist = c(which_hist,ind_where)
    data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
    data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
    
    if(change) {
      change_ind = c(tmp_deg_side_a_1,
                     tmp_deg_side_a_2,
                     tmp_deg_side_b_1,
                     tmp_deg_side_b_2, 
                     triangle_tmp_1, 
                     triangle_tmp_2, 
                     ind_where)
      data_present$change[change_ind] = 1
    }
    
  } else {
    
    # Repetition 
    # Increment the repetition by 1 for the tie that was observed  
    data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    if(change) {
      data_present$change[ind_where] = 1
    }
    
  }
  
  
  
  return(list(data_present = data_present, which_hist = which_hist))
}


simulation_rem_undirected_new = function(beta, n_actors, n,exo_cov, exo_cat, seed = 123){
  time_beg = Sys.time()
  set.seed(seed)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  obs = data.table(from =numeric(length = n),from =numeric(length = n))
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  lambda =  numeric(length = nrow(data_present))
  event_data = data.table("side_a" = numeric(length = n), "side_b"= numeric(length = n), "times" = numeric(length = n))
  which_hist = c()
  
  match_covs = match(names(beta),names(data_present) )
  
  
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda  = as.vector(exp(beta %*% t(data_present[,match_covs, with = F])))
    sum_lambda = sum(lambda)
    prob_lambda = lambda/sum_lambda
    waiting_time = rexp(n = 1, rate = sum_lambda)
    event = sample(x = n_possible,size = 1,prob = prob_lambda)
    # Update the covariates in data_present accordingly 
    side_a = data_present$side_a[event]
    side_b = data_present$side_b[event]
    event_data[i,"side_a"] = side_a
    event_data[i,"side_b"] = side_b 
    event_data[i,"times"] = waiting_time
    
    tmp = update_covs_undirected(data_present =data_present,which_hist =  which_hist,
                                 side_a =side_a, side_b =  side_b, n_actors =n_actors )
    data_present = tmp$data_present
    which_hist= tmp$which_hist
  }
  event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(list("data"= event_data, "end_stage" = data_present))
}

simulation_rem_undirected_list = function(beta, n_actors, n,exo_cov, exo_cat, seed = 123){
  time_beg = Sys.time()
  set.seed(seed)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  obs = data.table(from =numeric(length = n),from =numeric(length = n))
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  
  if(number_cont!=0){
    names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  match_covs = match(names(beta),names(data_present) )
  lambda =  numeric(length = nrow(data_present))
  event_data = data.table("side_a" = numeric(length = n), "side_b"= numeric(length = n), "times" = numeric(length = n))
  which_hist = c()
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda  = as.vector(exp(beta %*% t(data_present[,match_covs, with = F])))
    sum_lambda = sum(lambda)
    prob_lambda = lambda/sum_lambda
    waiting_time = rexp(n = 1, rate = sum_lambda)
    event = sample(x = n_possible,size = 1,prob = prob_lambda)
    # Update the covariates in data_present accordingly 
    side_a = data_present$side_a[event]
    side_b = data_present$side_b[event]
    event_data[i,"side_a"] = side_a
    event_data[i,"side_b"] = side_b 
    event_data[i,"times"] = waiting_time
    
    tmp = update_covs_undirected(data_present =data_present,which_hist =  which_hist,
                                 side_a =side_a, side_b =  side_b, n_actors =n_actors )
    data_present = tmp$data_present
    which_hist= tmp$which_hist
  }
  event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(list("data"= event_data, "end_stage" = data_present))
}


simulation_rem_undirected = function(beta, n_actors, n,exo_cov, exo_cat, seed = 123){
  time_beg = Sys.time()
  set.seed(seed)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  obs = data.table(from =numeric(length = n),from =numeric(length = n))
  # What are the possibles combinations one could observe? 
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  lambda =  numeric(length = nrow(data_present))
  event_data = data.table("side_a" = numeric(length = n), "side_b"= numeric(length = n), "times" = numeric(length = n))
  which_hist = c()
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda  = as.vector(exp(beta %*% t(data_present[,-c(1,2,10,11)])))
    sum_lambda = sum(lambda)
    prob_lambda = lambda/sum_lambda
    waiting_time = rexp(n = 1, rate = sum_lambda)
    event = sample(x = n_possible,size = 1,prob = prob_lambda)
    # Update the covariates in data_present accordingly 
    side_a = data_present$side_a[event]
    side_b = data_present$side_b[event]
    event_data[i,"side_a"] = side_a
    event_data[i,"side_b"] = side_b 
    event_data[i,"times"] = waiting_time
    
    #Step b: Indicate that the event has happenend at least once

    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    if(data_present$repetition[ind_where] == 0) {
      #Check which(data_present$from == from)
      tmp_deg_side_a_1 = find_side_a(a = side_a,x = n_actors)
      tmp_deg_side_a_2 =find_side_b(b = side_a,x = n_actors)
      tmp_deg_side_b_1 = find_side_a(a = side_b,x = n_actors)
      tmp_deg_side_b_2 =find_side_b(b = side_b,x = n_actors)
    
      data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
      data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
      data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
      data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
      
      # Transitivity
      # Increment 1 the transitivity if h was in sender j and may be in sender i 
      # or h was in receiver i and will be in receiver j 

      # What events already happened with side_a?
      ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist)
      ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist)
      # with whom were those events? 
      triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
      
      # What events already happened with side_b?
      ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist)
      ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist)
      # with whom were those events? 
      triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
      
      if(length(triangle_tmp_1)!=0) {
        side_a_tmp = triangle_tmp_1
        side_b_tmp = rep(side_b, length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }
      
      if(length(triangle_tmp_2)!=0) {
        side_a_tmp = triangle_tmp_2
        side_b_tmp = rep(side_a, length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }

      data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
     
       data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      which_hist = c(which_hist,ind_where)
      data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
      data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
    } else {
     
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
    }
  }
  event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(list("data"= event_data, "end_stage" = data_present))
}


# event_data = test_undirected$data
count_rem_undirected = function(event_data,exo_cov,exo_cat, n_actors, add_var = 1,ind = 1){
  time_beg = Sys.time()
  n = nrow(event_data)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  data_present$change = 1
  data_present$status = data_present$time = 0
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)

  data_ges = vector(mode = "list", length = nrow(event_data))
  data_ges[[1]] = data_present
  
  time_tmp = 0
  which_hist = c()
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    side_a = event_data$side_a[i] 
    side_b = event_data$side_b[i]  
    waiting_time = event_data$times[i] 
    
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    data_present$status[ind_where] = 1
    data_present$time = waiting_time
    # data_present$time[ind_where] = waiting_time
    # data_present$time[-ind_where] = time_tmp
    # time_tmp = waiting_time
    
    
    
    data_present$change = 0
    data_present$change[ind_where] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    if(data_present$repetition[ind_where] == 0) {
      #Check which(data_present$from == from)

      
      tmp_deg_side_a_1 = find_side_a(a = side_a,x = n_actors)
      tmp_deg_side_a_2 =find_side_b(b = side_a,x = n_actors)
      tmp_deg_side_b_1 = find_side_a(a = side_b,x = n_actors)
      tmp_deg_side_b_2 =find_side_b(b = side_b,x = n_actors)
    
      
      data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
      data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
      data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
      data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
      
      # Transitivity
      # Increment 1 the transitivity if h was in sender j and may be in sender i 
      # or h was in receiver i and will be in receiver j 
      # What events already happened with side_a?
      ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist)
      ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist)
      # with whom were those events? 
      triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
      
      # What events already happened with side_b?
      ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist)
      ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist)
      # with whom were those events? 
      triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
      
      if(length(triangle_tmp_1)!=0) {
        side_a_tmp = triangle_tmp_1
        side_b_tmp = rep(side_b, length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }
      
      if(length(triangle_tmp_2)!=0) {
        side_a_tmp = triangle_tmp_2
        side_b_tmp = rep(side_a, length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }
      
      data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
      

      
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
      data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
      
      change_ind = c(tmp_deg_side_a_1,
                     tmp_deg_side_a_2,
                     tmp_deg_side_b_1,
                     tmp_deg_side_b_2, 
                     triangle_tmp_1, 
                     triangle_tmp_2, 
                     ind_where)
      
      data_present$change[change_ind] = 1
      
      # Update the history 
      which_hist = c(which_hist,ind_where)
      
    } else {
      
      # Increment the repetition by 1 for the tie that was observed  
      
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      data_present$change[ind_where] = 1
    }
    
    data_ges[[i+1]] = data_present[change == 1]
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$side_a, concat_data$side_b,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = max(event_data$times)
  tmp = concat_data[,.(time_end =  c(time[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != time]
  # close(pb)
  diff = Sys.time() -time_beg
  # Add info that we need to say if the counted events are true (1) or false (0) in add_var, to indicate the simulation run (ind), and to calculate the offset 
  concat_data$add_var = add_var
  concat_data$ind = ind
  concat_data$offset = log(concat_data$time_end-concat_data$time)
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}

count_rem_undirected_list = function(event_data,exo_cov,exo_cat, n_actors, add_var = 1,ind = 1){
  time_beg = Sys.time()
  n = nrow(event_data)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  # Stop if there is no event data
  if(nrow(event_data) == 0){
    return(event_data)
  }
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  
  if(number_cont!=0){
    names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
      }
   
  }
  
  if(number_cat!=0){
    names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  # data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  # data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  data_present$change = 1
  data_present$status = data_present$time = 0
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  data_ges = vector(mode = "list", length = nrow(event_data))
  data_ges[[1]] = data_present
  
  time_tmp = 0
  which_hist = c()
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    side_a = event_data$side_a[i] 
    side_b = event_data$side_b[i]  
    waiting_time = event_data$times[i] 
    
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    data_present$status[ind_where] = 1
    data_present$time = waiting_time
    # data_present$time[ind_where] = waiting_time
    # data_present$time[-ind_where] = time_tmp
    # time_tmp = waiting_time
    
    
    
    data_present$change = 0
    data_present$change[ind_where] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    tmp = update_covs_undirected(data_present = data_present, which_hist = which_hist,
                           side_a = side_a, side_b = side_b, n_actors =n_actors, change = T)
    
    data_present = tmp$data_present
    which_hist = tmp$which_hist
    
    data_ges[[i+1]] = data_present[change == 1]
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$side_a, concat_data$side_b,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = max(event_data$times)
  tmp = concat_data[,.(time_end =  c(time[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != time]
  # close(pb)
  diff = Sys.time() -time_beg
  # Add info that we need to say if the counted events are true (1) or false (0) in add_var, to indicate the simulation run (ind), and to calculate the offset 
  concat_data$add_var = add_var
  concat_data$ind = ind
  concat_data$offset = log(concat_data$time_end-concat_data$time)
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}


count_rem_undirected_new = function(event_data,exo_cov,exo_cat, n_actors, add_var = 1,ind = 1){
  time_beg = Sys.time()
  n = nrow(event_data)
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  # Stop if there is no event data
  if(nrow(event_data) == 0){
    return(event_data)
  }
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$cov_cont = exo_cov[data_present$side_a] + exo_cov[data_present$side_b]
  data_present$cov_cat = as.numeric(exo_cat[data_present$side_a] == exo_cat[data_present$side_b])
  
  # data_present$cov_cont_side_a = exo_cov[data_present$from]
  # data_present$cov_cont_side_b = exo_cov[data_present$to]
  # data_present$cov_cat_side_a = exo_cat[data_present$from]
  # data_present$cov_cat_side_b = exo_cat[data_present$to]  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  data_present$change = 1
  data_present$status = data_present$time = 0
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  
  data_ges = vector(mode = "list", length = nrow(event_data))
  data_ges[[1]] = data_present
  
  time_tmp = 0
  which_hist = c()
  for(i in 1:n){
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    side_a = event_data$side_a[i] 
    side_b = event_data$side_b[i]  
    waiting_time = event_data$times[i] 
    
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_where = findrow_undirected(a = side_a,b = side_b,x = n_actors)
    data_present$status[ind_where] = 1
    data_present$time = waiting_time
    # data_present$time[ind_where] = waiting_time
    # data_present$time[-ind_where] = time_tmp
    # time_tmp = waiting_time
    
    
    
    data_present$change = 0
    data_present$change[ind_where] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    tmp = update_covs_undirected(data_present = data_present, which_hist = which_hist,
                                 side_a = side_a, side_b = side_b, n_actors =n_actors, change = T)
    
    data_present = tmp$data_present
    which_hist = tmp$which_hist
    
    data_ges[[i+1]] = data_present[change == 1]
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$side_a, concat_data$side_b,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = max(event_data$times)
  tmp = concat_data[,.(time_end =  c(time[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != time]
  # close(pb)
  diff = Sys.time() -time_beg
  # Add info that we need to say if the counted events are true (1) or false (0) in add_var, to indicate the simulation run (ind), and to calculate the offset 
  concat_data$add_var = add_var
  concat_data$ind = ind
  concat_data$offset = log(concat_data$time_end-concat_data$time)
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}


count_rem = function(event_data,exo_cov, n_actors){
  time_beg = Sys.time()
  n = nrow(event_data)
  pb <- txtProgressBar(min = 0, max =n, style = 3)

  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  n_possible = nrow(data_present)
  
  names(data_present) = c("from","to")
  data_present$intercept= 1
  data_present$out_sender= numeric(length = n_possible)
  data_present$in_sender= numeric(length = n_possible)
  data_present$out_receiver= numeric(length = n_possible)
  data_present$in_receiver= numeric(length =n_possible)
  data_present$reciprocity= numeric(length = n_possible)
  data_present$repetition= numeric(length = n_possible)
  data_present$transitivity= numeric(length = n_possible)
  data_present$shared_contactees= numeric(length = n_possible)
  data_present$triangle= numeric(length = n_possible)
  data_present$shared_contacters= numeric(length = n_possible)
  data_present$exo_from = exo_cov[data_present$from]
  data_present$exo_to = exo_cov[data_present$to]
  data_present$change = 1
  data_present$status = data_present$time = 0
  data_ges = vector(mode = "list", length = nrow(event_data))
  data_ges[[1]] = data_present
  
  time_tmp = 0
  which_hist = c()
  for(i in 1:n){
    setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    from = event_data$from[i] 
    to = event_data$to[i]  
    waiting_time = event_data$times[i] 
    
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_where = findrow(a = from,b = to,x = n_actors)
    data_present$status[ind_where] = 1
    data_present$time = waiting_time
    # data_present$time[ind_where] = waiting_time
    # data_present$time[-ind_where] = time_tmp
    # time_tmp = waiting_time
    
   
    
    data_present$change = 0
    data_present$change[ind_where] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    if(data_present$repetition[ind_where] == 0) {
      #Check which(data_present$from == from)
      tmp_out_sender_ind = findfrom(a = from,x = n_actors)
      tmp_in_sender_ind = findfrom(a = to,x = n_actors)
      tmp_out_receiver_ind =  findto(b = from,x = n_actors)
      tmp_in_receiver_ind =findto(b = to,x = n_actors)
   
      data_present$out_sender[tmp_out_sender_ind] = data_present$out_sender[tmp_out_sender_ind] + 1
      
      # Outdegree Receiver
      # Increment 1 the outdegree receiver for all events that may take place with i being the receiver
      
      #Check which(data_present$to == from)
      data_present$out_receiver[tmp_out_receiver_ind] = data_present$out_receiver[tmp_out_receiver_ind] +1
      
      # Indegree Sender
      # Increment 1 the indegree sender for all events that may take place with j being the sender
      #Check which(data_present$from == to)
      data_present$in_sender[tmp_in_sender_ind] = data_present$in_sender[tmp_in_sender_ind] +1
      
      # Indegree Receiver
      # Increment 1 the Indegree receiver for all events that may take place with j being the receiver
      #Check which(data_present$to == to)
      data_present$in_receiver[tmp_in_receiver_ind] = data_present$in_receiver[tmp_in_receiver_ind] +1
      
      # Transitivity
      # Increment 1 the transitivity if h was in sender j and may be in sender i 
      # or h was in receiver i and will be in receiver j 
      
      # What actors did to already send to? 
      # which((data_present$from == to) & (data_present$history == TRUE))
      ind_sender_to_j = intersect(tmp_in_sender_ind, which_hist)
      
      # Whats the position of from to those h 
      tmp = data_present$to[ind_sender_to_j]
      tmp = tmp[tmp!= from]
      ind_transitivity_1 = findrows_1(bs = tmp,a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$to[which((data_present$from == to) & (data_present$history == TRUE))] )
      
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(tmp_out_receiver_ind, which_hist)
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_transitivity_2 = findrows_2(as =  tmp,b = to, x = n_actors)
      
      #which(data_present$to == to & data_present$from %in% data_present$from[which((data_present$to == from) & (data_present$history == TRUE))])
      ind_transitivity = unique(c(ind_transitivity_1,ind_transitivity_2))
      ind_transitivity = ind_transitivity[ind_transitivity<= n_actors * (n_actors-1)]
      
      data_present$transitivity[ind_transitivity] = data_present$transitivity[ind_transitivity] +1
      
      

      # Shared contactees
      # Increment 1 shared contactees if h was in receiver j in the past and h may be sender future of i 
      # or if h was in receiver of j and will be in receiver of i
      
      #What possible ties did already send to j in the past 
      ind_sender_to_j = intersect(tmp_in_receiver_ind, which_hist)
      ind_sender_to_j = ind_sender_to_j[ind_sender_to_j!= ind_where]
      # which((data_present$to == to) & (data_present$history == TRUE))
      # Which ties may i send to of these 
      
      ind_shared_contactees_1 = findrows_1(bs =  data_present$from[ind_sender_to_j],
                                           a = from, x = n_actors)
      
      # which(data_present$from == from & data_present$to %in% data_present$from[which((data_present$to == to) & (data_present$history == TRUE))])
      # Which of these may send to i 
      
      ind_shared_contactees_2 = findrows_2(as =  data_present$from[ind_sender_to_j],
                                           b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$from[ind_sender_to_j])
      
      ind_shared_contactees = unique(c(ind_shared_contactees_1,ind_shared_contactees_2))
      ind_shared_contactees = ind_shared_contactees[ind_shared_contactees<= n_actors * (n_actors-1)]
      data_present$shared_contactees[ind_shared_contactees] = data_present$shared_contactees[ind_shared_contactees] +1
      

      
      # Triangle Closure 
      # Increment 1 the triangle closure if h was in receiver i and may be in sender j 
      # or h was in sender j and will be in receiver i 
      #What possible ties were already sent from j 
      ind_sender_from_j =   intersect(which_hist,tmp_in_sender_ind)
      ind_sender_from_j = ind_sender_from_j[ind_sender_from_j!= ind_where]
      
      # which((data_present$from == to) & (data_present$history == TRUE))
      # Which of these ties have i as the receiver? 
      tmp = data_present$to[ind_sender_from_j]
      tmp = tmp[tmp!= from]
      ind_triangle_1 =  findrows_2(as = tmp,
                                   b = from, x = n_actors)
      
      # which(data_present$to == from & data_present$from %in% data_present$to[ind_sender_from_j])
      # What possible ties did already send to i in the past
      ind_sender_to_i= intersect(which_hist,tmp_out_receiver_ind)
      ind_sender_to_i = ind_sender_to_i[ind_sender_to_i!= ind_where]
      
      # which((data_present$to == from) & (data_present$history == TRUE))
      tmp = data_present$from[ind_sender_to_i]
      tmp = tmp[tmp!= to]
      
      ind_triangle_2 = findrows_1(bs =  tmp,
                                  a = to, x = n_actors)
      # which(data_present$from == to & data_present$to %in% data_present$from[ind_sender_to_i])
      ind_triangle = unique(c(ind_triangle_1,ind_triangle_2))
      ind_triangle = ind_triangle[ind_triangle<= n_actors * (n_actors-1)]
      
      data_present$triangle[ind_triangle] = data_present$triangle[ind_triangle] +1
      

      
      
      # Shared contacters
      # Increment 1 shared contacters if h was in sender i and will be in sender j 
      # or h was in sender i and will be in receiver j 
      #What possible ties were already sent from i 
      ind_sender_from_i = intersect(which_hist,tmp_out_sender_ind)
      ind_sender_from_i = ind_sender_from_i[ind_sender_from_i!= ind_where]
      
      # which((data_present$from == from) & (data_present$history == TRUE))
      # Which ties may i send to
      
      
      ind_shared_contacters_1 = findrows_1(bs =  data_present$to[ind_sender_from_i],
                                           a = to, x = n_actors)
      
      # which(data_present$from == to & data_present$to %in% data_present$to[ind_sender_from_i])
      
      ind_shared_contacters_2 = findrows_2(as = data_present$to[ind_sender_from_i],
                                           b = to, x = n_actors)
      
      # which(data_present$to == to & data_present$from %in% data_present$to[ind_sender_from_i])
      ind_shared_contacters = unique(c(ind_shared_contacters_1,ind_shared_contacters_2))
      ind_shared_contacters = ind_shared_contacters[ind_shared_contacters<= n_actors * (n_actors-1)]
      data_present$shared_contacters[ind_shared_contacters] = data_present$shared_contacters[ind_shared_contacters] +1 
     
      
      
       # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      change_ind = c(tmp_out_sender_ind,
                     tmp_in_sender_ind,
                     tmp_out_receiver_ind,
                     tmp_in_receiver_ind, 
                     ind_transitivity, 
                     ind_shared_contactees, 
                     ind_triangle, 
                     ind_shared_contacters, 
                     tmp_reciproc_ind, 
                     ind_where)
      
      data_present$change[change_ind] = 1
      
      # Update the history 
      which_hist = c(which_hist,ind_where)
      
    } else {
      # Reciprocity
      # Increment 1 the reciprocity for the inverse receiver/sender tie
      tmp_reciproc_ind = findrow(a = to,b = from,x = n_actors)
      #Check which((data_present$to == from) &  (data_present$from == to))
      data_present$reciprocity[tmp_reciproc_ind] = data_present$reciprocity[tmp_reciproc_ind] +1
     
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_where] = data_present$repetition[ind_where]+1
      data_present$change[c(tmp_reciproc_ind,ind_where)] = 1
    }
    data_ges[[i+1]] = data_present[change == 1]
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$from, concat_data$to,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = max(event_data$times)
  tmp = concat_data[,.(time_end =  c(time[-1], max_time), 
                 status = c(status[-1],0)), 
              by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != time]
  close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}



count_rem_new = function(event_data,exo_cov, n_actors){
  time_beg = Sys.time()
  n = nrow(event_data)
  pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  # This is a dataframe that changes and gives us the state of the covariates
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  n_possible = nrow(data_present)
  
  names(data_present) = c("from","to")
  data_present$intercept= 1
  data_present$out_sender= numeric(length = n_possible)
  data_present$in_sender= numeric(length = n_possible)
  data_present$out_receiver= numeric(length = n_possible)
  data_present$in_receiver= numeric(length =n_possible)
  data_present$reciprocity= numeric(length = n_possible)
  data_present$repetition= numeric(length = n_possible)
  data_present$transitivity= numeric(length = n_possible)
  data_present$shared_contactees= numeric(length = n_possible)
  data_present$triangle= numeric(length = n_possible)
  data_present$shared_contacters= numeric(length = n_possible)
  data_present$exo_from = exo_cov[data_present$from]
  data_present$exo_to = exo_cov[data_present$to]
  data_present$change = 1
  data_present$status = data_present$time = 0
  data_ges = vector(mode = "list", length = nrow(event_data))
  data_ges[[1]] = data_present
  
  time_tmp = 0
  which_hist = c()
  for(i in 1:n){
    setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    from = event_data$from[i] 
    to = event_data$to[i]  
    waiting_time = event_data$times[i] 
    
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_where = findrow(a = from,b = to,x = n_actors)
    data_present$status[ind_where] = 1
    data_present$time = waiting_time
    # data_present$time[ind_where] = waiting_time
    # data_present$time[-ind_where] = time_tmp
    # time_tmp = waiting_time
    
    
    
    data_present$change = 0
    data_present$change[ind_where] = 1
    
    tmp = update_covs_directed(data_present = data_present, which_hist = which_hist,from = from ,to = to, n_actors =n_actors)
    
    data_present = tmp$data_present
    which_hist = tmp$which_hist
    data_ges[[i+1]] = data_present[change == 1]
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$from, concat_data$to,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = max(event_data$times)
  tmp = concat_data[,.(time_end =  c(time[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != time]
  close(pb)
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}




