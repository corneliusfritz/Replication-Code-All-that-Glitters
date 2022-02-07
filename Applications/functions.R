
summary.da_ecrem_result = function(model) {
  names <- colnames(model$betas_real)
  which_exclude = grepl(pattern = "time",x = names)
  names = names[!which_exclude]
  names = gsub("repetition_binFALSE:", "Onset ",x = names)
  names = gsub("repetition_binTRUE:", "Continuation ",x = names)
  names = gsub("\\(Intercept\\)", "Intercept",x = names)
  # names = gsub("repetition_binTRUE", "Continuation Intercept",x = names)
  names = gsub("cov_", "z_",x = names)
  
  coefs_real <- model$beta_final[!which_exclude]
  coefs_real[which(names == "Continuation Intercept")] = coefs_real[which(names == "Continuation Intercept")] + coefs_real[which(names == "Onset Intercept")]
  
  se_real <-  sqrt(diag(model$vcov_real)[!which_exclude])
  se_real[which(names == "Continuation Intercept")] =
    sqrt(diag(model$vcov_real)[!which_exclude][which(names == "Continuation Intercept")]+
           diag(model$vcov_real)[!which_exclude][which(names == "Onset Intercept")] + 
           2*model$vcov_real[!which_exclude,!which_exclude][which(names == "Continuation Intercept"),which(names == "Onset Intercept") ])
  
  alpha = 0.05
  lower_real = coefs_real - se_real*qnorm(1-0.5*alpha)
  upper_real = coefs_real + se_real*qnorm(1-0.5*alpha)
  which_exclude = grepl(pattern = "time",x = names(model$betas_fake))
  coefs_fake <- model$betas_fake
  coefs_fake <-  apply(coefs_fake,  2,mean )
  coefs_fake = coefs_fake[grepl(pattern = "s\\(time\\)",x = names(coefs_fake))]
  names = gsub("Intercept", "AIntercept",x = names)
  
  order_tmp = order(names)
  names = gsub("AIntercept","Intercept",x = names)
  names = names[order_tmp]
  coefs_real = coefs_real[order_tmp]
  se_real = se_real[order_tmp]
  lower_real = lower_real[order_tmp]
  upper_real = upper_real[order_tmp]
  names = gsub("Continuation ",replacement = ' ',x = names)
  names = gsub("Onset ","",x = names)
  names = gsub("z_dyad_","",x = names)
  names = gsub("z_cat_","",x = names)
  names= gsub("_"," ",x = names)
  names = str_to_title(names)
  
  
  return(data.table(coef.names = names,
                    coef = coefs_real,
                    se = se_real,
                    ci.low = lower_real,
                    ci.up = upper_real
  )
  )
}
fun_mi_ecrem_result = function(model,alpha = 0.05, digits = 3) {
  names <- names(model$beta_final)
  which_exclude = grepl(pattern = "time",x = names)
  names = names[!which_exclude]
  
  # names = gsub("repetition_binFALSE:", "Onset ",x = names)
  # names = gsub("repetition_binTRUE:", "Continuation ",x = names)
  # names = gsub("\\(Intercept\\)", "Onset Intercept",x = names)
  # names = gsub("repetition_binTRUE", "Continuation Intercept",x = names)
  names = gsub("cov_", "z_",x = names)
  
  coefs_real <- model$beta_final[!which_exclude]
  # coefs_real[which(names == "Continuation Intercept")] = coefs_real[which(names == "Continuation Intercept")] + coefs_real[which(names == "Onset Intercept")]
  
  se_real <-  sqrt(diag(model$vcov_real)[!which_exclude])
  # se_real[which(names == "Continuation Intercept")] =
  #   sqrt(diag(model$vcov)[!which_exclude][which(names == "Continuation Intercept")]+
  #          diag(model$vcov)[!which_exclude][which(names == "Onset Intercept")] + 
  #          2*model$vcov[!which_exclude,!which_exclude][which(names == "Continuation Intercept"),which(names == "Onset Intercept") ])
  # 
  
  lower_real = coefs_real - se_real*qnorm(1-0.5*alpha)
  upper_real = coefs_real + se_real*qnorm(1-0.5*alpha)
  which_exclude = grepl(pattern = "time",x = names(model$beta_fake_final))
  coefs_fake <- model$beta_fake_final
  # coefs_fake <-  apply(coefs_fake,  2,mean )
  # coefs_fake = coefs_fake[grepl(pattern = "s\\(time\\)",x = names(coefs_fake))]
  names = gsub("Intercept", "AIntercept",x = names)
  
  order_tmp = order(names)
  names = gsub("AIntercept","Intercept",x = names)
  names = names[order_tmp]
  coefs_real = coefs_real[order_tmp]
  se_real = se_real[order_tmp]
  lower_real = lower_real[order_tmp]
  upper_real = upper_real[order_tmp]
  names = gsub("Continuation ",replacement = ' ',x = names)
  names = gsub("Onset ","",x = names)
  names = gsub("z_dyad_","",x = names)
  names = gsub("z_cat_","",x = names)
  names= gsub("_"," ",x = names)
  names = str_to_title(names)
  
  tr <- data.table(coef.names = names,
                   coef = coefs_real,
                   se = se_real, 
                   t_vals = coefs_real/se_real,
                   ci.low = lower_real,
                   ci.up = upper_real
                   # gof =c(NA,coefs_fake),
                   # gof.decimal = c(NA,T),
                   # gof.names = c("Measurement Error",names(coefs_fake))
                   
                   
  )
  
  tr$CI = trimws(paste0("[",round(tr$ci.low, digits = digits), ",",round(tr$ci.up, digits = digits),"]"))
  tr$ci.up = tr$ci.low = NULL
  tr$coef = special_round(tr$coef, digits = digits)
  tr$se = special_round(tr$se, digits = digits)
  tr$t_vals = special_round(tr$t_vals, digits = digits)
  return(tr)
}

fun_rem_result = function(model, digits = 3, alpha = 0.05) {
  names <- names(model$coefficients)
  which_exclude = grepl(pattern = "time",x = names)
  names = names[!which_exclude]
  # names = gsub("repetition_binFALSE:", "Onset ",x = names)
  # names = gsub("repetition_binTRUE:", "Continuation ",x = names)
  # names = gsub("\\(Intercept\\)", "Onset Intercept",x = names)
  # names = gsub("repetition_binTRUE", "Continuation Intercept",x = names)
  names = gsub("cov_", "z_",x = names)
  
  coefs_real <- model$coefficients[!which_exclude]
  # coefs_real[which(names == "Continuation Intercept")] = coefs_real[which(names == "Continuation Intercept")] + coefs_real[which(names == "Onset Intercept")]
  
  se_real <-  sqrt(diag(model$Ve)[!which_exclude])
  # se_real[which(names == "Continuation Intercept")] =
  #   sqrt(diag(model$Ve)[!which_exclude][which(names == "Continuation Intercept")]+
  #          diag(model$Ve)[!which_exclude][which(names == "Onset Intercept")] + 
  #          2*model$Ve[!which_exclude,!which_exclude][which(names == "Continuation Intercept"),which(names == "Onset Intercept") ])
  # 
  lower_real = coefs_real - se_real*qnorm(1-0.5*alpha)
  upper_real = coefs_real + se_real*qnorm(1-0.5*alpha)
  
  names = gsub("Intercept", "AIntercept",x = names)
  
  order_tmp = order(names)
  names = gsub("AIntercept","Intercept",x = names)
  names = names[order_tmp]
  coefs_real = coefs_real[order_tmp]
  se_real = se_real[order_tmp]
  lower_real = lower_real[order_tmp]
  upper_real = upper_real[order_tmp]
  names = gsub("Continuation "," ",x = names)
  names = gsub("Onset ","",x = names)
  names = gsub("z_dyad_","",x = names)
  names = gsub("z_cat_","",x = names)
  names= gsub("_"," ",x = names)
  names = stringr::str_to_title(names)
  
  tr <- data.table(coef.names = names,
                   coef = coefs_real,
                   se = se_real, 
                   t_vals = coefs_real/se_real,
                   ci.low = lower_real,
                   ci.up = upper_real)
  
  tr$CI = trimws(paste0("[",round(tr$ci.low, digits = digits), ",",round(tr$ci.up, digits = digits),"]"))
  tr$ci.up = tr$ci.low = NULL
  tr$coef = special_round(tr$coef, digits = digits)
  tr$se = special_round(tr$se, digits = digits)
  tr$t_vals = special_round(tr$t_vals, digits = digits)
  
  return(tr)
}




special_round = function(x, digits = 3){
  x = round(x, digits = digits)
  x[x == 0] = "< 0.001"
  return(x)
}

mecm_iteration = function(M,cl = NULL,  
                          build = F, pi, 
                          end_time,
                          beg_time,
                          comb_event_stream, 
                          formula_real, 
                          formula_fake,
                          n_actors, 
                          change_all,
                          mod_real, 
                          mod_fake, 
                          which_covs_real = which_covs_real, 
                          which_covs_fake = which_covs_fake,seed, ind_important = F,
                          exo_cov, exo_cat, exo_dyad,
                          exo_cov_fake, exo_cat_fake, exo_dyad_fake,
                          bam = bam , ...){
  if(is.null(cl)){
    result = lapply(1:(M), e_step_compl_clustered_mgcv, 
                    end_time = end_time, beg_time = beg_time, 
                    comb_event_stream = comb_event_stream, 
                    n_actors = n_actors, pi = pi, 
                    change_all = change_all,
                    mod_real = mod_real, mod_fake =mod_fake, 
                    which_covs_real = which_covs_real, 
                    which_covs_fake = which_covs_fake,seed = seed, ind_important = F,
                    exo_cov = exo_cov, exo_cat = exo_cat, exo_dyad = exo_dyad,
                    exo_cov_fake = exo_cov_fake, exo_cat_fake = exo_cat_fake, exo_dyad_fake = exo_dyad_fake, 
                    bam = bam)
  } else {
    result = parLapply(cl, 1:(M), e_step_compl_clustered_mgcv, 
                       end_time = end_time, beg_time = beg_time, 
                       comb_event_stream = comb_event_stream, 
                       n_actors = n_actors,  pi = pi, 
                       change_all = change_all,
                       mod_real = mod_real, mod_fake =mod_fake, 
                       which_covs_real = which_covs_real, 
                       which_covs_fake = which_covs_fake,seed = seed, ind_important = F,
                       exo_cov = exo_cov, exo_cat = exo_cat, exo_dyad = exo_dyad,
                       exo_cov_fake = exo_cov_fake, exo_cat_fake = exo_cat_fake, exo_dyad_fake = exo_dyad_fake, 
                       bam = bam)
  }
  
  sampled_events = result[[1]]$sampled_events
  pi_new = mean(sampled_events$number_real)
  # build_time = comb_event_stream$times[100]
  # Should there be some build-up time?
  # if(build) {
  #   result$include =  (result$time >= build_time) +  (result$time_end > build_time)
  #   result$time[result$include == 1] = build_time
  #   result = result[include!= 0]
  # }
  # if(M>1){
  #   result = result[,.(weight = sum(weight)),by = c(names(result)[-length(names(result))])]
  # }
  # 
  count_real =  result[[1]]$count_real
  # concat_data = concat_data[,.(weight = .N),by = c(names(concat_data)[c(3:20,25)])]
  
  count_fake = result[[1]]$count_fake
  count_real$repetition_bin = factor(count_real$repetition_bin ) 
  attr(formula_real,'.Environment') = environment()
  attr(formula_fake,'.Environment') = environment()
  
  if(bam[1]){
    mod_real =   bam(formula_real,family = poisson(), data = count_real,
                     weights = count_real$weight, discrete = T, nthreads = 20,
                     samfrac = 0.1, ...)
  } else {
    mod_real =   glm(formula_real,family = poisson(), 
                     data = count_real,weights = count_real$weight, ...)
  }
  
  if(bam[2]){
    mod_fake =  bam(formula_fake,family = poisson(), data = count_fake, 
                    discrete = T, nthreads = 20,weights = count_fake$weight, ...)
  } else {
    mod_fake =   glm(formula_fake,family = poisson(), data = count_fake,weights = count_fake$weight, ...)
  }
  gc(full = T)
  
  # if(max(abs(c(coef(mod_real),coef(mod_fake))))>100){
  #   browser()
  # }
  return(list(mod_real = mod_real, 
              mod_fake = mod_fake, 
              sampled_events = sampled_events, 
              pi_new = pi_new))
  
}


da_estimation_clustered_mgcv = function(K = 10, K_2 = 40, beta_real_gt,beta_fake_gt, full_bayes = F, 
                                        comb_event_stream = comb_event_stream_clustered, separable = T, 
                                        n_actors = n_actors, eps = 0.1,
                                        change_all = "lazy",
                                        rand_prop_real = 0.8,
                                        mod_real = NULL,
                                        mod_fake = NULL,
                                        beg_time = NULL, 
                                        end_time = NULL,
                                        exo_cov = list(), 
                                        exo_cat = list(),
                                        exo_dyad= list(), 
                                        exo_cov_fake = list(),
                                        exo_dyad_fake = list(), 
                                        exo_cat_fake  = list(),
                                        beg_formula_real = "status ~  repetition_bin+ offset(offset)+  s(time, by =repetition_bin,  bs = 'ps')+", 
                                        beg_formula_fake= "status ~  s(time, bs = 'ps')  + offset(offset)", 
                                        seed, build, build_time = NULL,only_beg = F, ...){
  # Setup----
  set.seed(seed)
  call = sys.calls()[[sys.nframe()]]
  stopifnot(change_all %in% c("all", "week", "lazy"))
  if(is.null(beg_time)){
    beg_time = min(comb_event_stream$times)
  }
  if(is.null(end_time)){
    end_time = max(comb_event_stream$times)
  }
  which_covs_real =names(beta_real_gt)[beta_real_gt==T]
  # add the exogenous covariates 
  if(length(exo_cat) != 0){
    which_covs_real = c(which_covs_real,paste("cov_cat",names(exo_cat),sep = "_"))
  }
  if(length(exo_dyad) != 0){
    which_covs_real = c(which_covs_real,paste("cov_dyad",names(exo_dyad),sep = "_"))
  }
  if(length(exo_cov) != 0){
    which_covs_real = c(which_covs_real,paste("cov_cont",names(exo_cov),sep = "_"))
  }
  
  which_covs_fake = names(beta_fake_gt)[beta_fake_gt==T]  # add the exogenous covariates 
  if(length(exo_cat_fake) != 0){
    which_covs_fake = c(which_covs_fake,paste("cov_cat",names(exo_cat_fake),sep = "_"))
  }
  if(length(exo_dyad_fake) != 0){
    which_covs_fake = c(which_covs_fake,paste("cov_dyad",names(exo_dyad_fake),sep = "_"))
  }
  if(length(exo_cov_fake) != 0){
    which_covs_fake = c(which_covs_fake,paste("cov_cont",names(exo_cov_fake),sep = "_"))
  }
  
  if(separable) {
    formula_real = formula(paste(beg_formula_real,paste(which_covs_real[-1],
                                                        ":repetition_bin",collapse = " + ")))
    # formula_real = formula(paste("status ~  offset(offset)+ repetition_bin + s(time,  bs = 'ps')+",paste(which_covs_real[-1],
    #                                                                                                    ":repetition_bin",collapse = " + ")))
    # formula_real = formula(paste("status ~  offset(offset)+ repetition_bin + time+",paste(which_covs_real[-1],
    #                                                                                                      ":repetition_bin",collapse = " + ")))
    # formula_fake = formula(paste("status ~  1  + offset(offset)",paste(which_covs_fake[-1],collapse = " + ")))
    # formula_fake = formula(paste("status ~  time  + offset(offset)",paste(which_covs_fake[-1],collapse = " + ")))
    formula_fake = formula(paste(beg_formula_fake,paste(which_covs_fake[-1],collapse = " + ")))
    
    # formula_fake = formula(paste("status ~  time  + offset(offset)",paste(which_covs_fake[-1],collapse = " + ")))
  } else {
    if(full_bayes) {
      formula_real =     formula(paste(beg_formula_real,paste0("s(",which_covs_real[-1],
                                                      ", bs = 're')",collapse = " + ")))
      if(length(which_covs_fake)>1){
        formula_fake =  paste(beg_formula_fake,paste("s(",which_covs_fake[-1],
                                                     ", bs = 're')",collapse = " + "))
      } 
      else {
        formula_fake = formula(paste(beg_formula_fake,paste(which_covs_fake[-1],collapse = " + ")))
      }
      } 
    else {
      formula_real = formula(paste(beg_formula_real,paste(which_covs_real[-1],collapse = " + ")))
      formula_fake = formula(paste(beg_formula_fake,paste(which_covs_fake[-1],collapse = " + ")))
    }
    
  }
  
  bam = c(F,F)
  if(length(grep(pattern = "s\\(", x = formula_real))>0){
    bam[1] = T
  }
  if(length(grep(pattern = "s\\(", x = formula_fake))>0){
    bam[2] = T
  }
  
  # If only a REM should be estimated, we assume that all events are real and follow through with estimation
  if(only_beg){
    # Since the window is set really high the result is equivalent to  count_rem_undirected_clustered
    test_unchange = count_rem_undirected_windowed(event_data = comb_event_stream,window = 365, 
                                                   exo_cov = exo_cov,n_actors = n_actors,
                                                   exo_cat = exo_cat, exo_dyad = exo_dyad,
                                                   change_all = change_all, end = end_time,
                                                   start =  beg_time)
    if(build) {
      test_unchange$include =  (test_unchange$times >= build_time) +  (test_unchange$times_end > build_time)
      test_unchange$times[test_unchange$include == 1] = build_time
      test_unchange = test_unchange[include!= 0]
    }
    test_unchange$repetition_bin = factor(test_unchange$repetition_bin )
    if(bam[1]){
      mod_real =bam(formula_real,  family = poisson(),
                    weights = test_unchange$weight, data = test_unchange, discrete = T, nthreads = 20, ...)
   
    } else {
      mod_real =glm(formula_real,  family = poisson(),weights = test_unchange$weight, data = test_unchange, ...)
    }
    return(mod_real)
  }
  
  # Set starting parameters
  if(is.null(mod_real)){
    # For the real parameters we use the set of values that assumes all events are real 
    comb_event_stream$ind_real = F
    comb_event_stream$ind_real[sample(x = 1:length(comb_event_stream$ind_real),
                                      size = floor(length(comb_event_stream$ind_real)*rand_prop_real))] = T
    test_unchange_real = count_rem_undirected_clustered(event_data = comb_event_stream[ind_real == T],
                                                        exo_cov = exo_cov,n_actors = n_actors,
                                                        exo_cat = exo_cat, exo_dyad = exo_dyad,
                                                        change_all = change_all, end = end_time,
                                                        start =  beg_time)
    if(build) {
      test_unchange_real$include =  (test_unchange_real$time >= build_time) +  (test_unchange_real$time_end > build_time)
      test_unchange_real$time[test_unchange_real$include == 1] = build_time
      test_unchange_real = test_unchange_real[include!= 0]
    }
    test_unchange_real$repetition_bin = factor(test_unchange_real$repetition_bin ) 
    if(bam[1]){
      mod_real =bam(formula_real,  family = poisson(),
                    weights = test_unchange_real$weight, data = test_unchange_real, discrete = T, nthreads = 20, ...)
      
    } else {
      mod_real =glm(formula_real,  family = poisson(),weights = test_unchange_real$weight, data = test_unchange_real, ...)
    }
    test_unchange_fake = count_rem_undirected_clustered(event_data = comb_event_stream[ind_real == F],
                                                        exo_cov = exo_cov_fake,n_actors = n_actors,
                                                        exo_cat = exo_cat_fake, exo_dyad = exo_dyad_fake,
                                                        change_all = change_all, end = end_time,
                                                        start =  beg_time)
    if(build) {
      test_unchange_fake$include =  (test_unchange_fake$time >= build_time) +  (test_unchange_fake$time_end > build_time)
      test_unchange_fake$time[test_unchange_fake$include == 1] = build_time
      test_unchange_fake = test_unchange_fake[include!= 0]
    }
    test_unchange_fake$repetition_bin = factor(test_unchange_fake$repetition_bin ) 
    
    # browser()
    if(bam[2]){
      mod_fake =bam(formula_fake,  family = poisson(),
                    weights = test_unchange_fake$weight, data = test_unchange_fake, discrete = T, nthreads = 20, ...)
      
    } else {
      mod_fake =glm(formula_fake,  family = poisson(),weights = test_unchange_fake$weight, data = test_unchange_fake, ...)
    }
    rm(test_unchange_real)
    rm(test_unchange_fake)
    gc()
  }
  betas_fake = matrix(nrow = K+1 ,ncol = length(mod_fake$coefficients))
  betas_real = matrix(nrow = K+1,ncol = length(mod_real$coefficients))
  colnames(betas_fake) = names(mod_fake$coefficients)
  colnames(betas_real) = names(mod_real$coefficients)
  
  betas_real[1,] = mod_real$coefficients
  betas_fake[1,] =  mod_fake$coefficients
  
  pis = numeric(length = K + K_2)
  pis[1] = rand_prop_real
  beg_models = list(mod_real =  list(Vp = mod_real$Vp,
                                     Ve = mod_real$Ve , 
                                     Vc = mod_real$Vc,
                                     coefficients = mod_real$coefficients), 
                    mod_fake =   list(Vp = mod_fake$Vp,
                                      Ve = mod_fake$Ve , 
                                      Vc = mod_fake$Vc,
                                      coefficients = mod_real$coefficients))
  names_mod_real = names(coef(mod_real))
  names_mod_fake = names(coef(mod_fake))
  
  # Warm-Up Phase----
  for(k in 1:K) {
    # if(k == 14) {
    #   browser()
    # }
    iteration_tmp = mecm_iteration(comb_event_stream = comb_event_stream, 
                                   formula_real = formula_real, 
                                   formula_fake = formula_fake, 
                                   n_actors = n_actors, M = 1,
                                   change_all = change_all,
                                   mod_real = mod_real, 
                                   mod_fake =mod_fake, 
                                   end_time = end_time,
                                   beg_time =  beg_time,
                                   pi  =  pis[1],
                                   which_covs_real = which_covs_real, 
                                   which_covs_fake = which_covs_fake,
                                   seed = k + seed,
                                   ind_important = F,
                                   exo_cov = exo_cov, 
                                   exo_cat = exo_cat, 
                                   exo_dyad = exo_dyad,
                                   exo_cov_fake = exo_cov_fake, 
                                   exo_cat_fake = exo_cat_fake, 
                                   exo_dyad_fake = exo_dyad_fake,
                                   bam = bam, ...)
    
    mod_real = iteration_tmp$mod_real
    mod_fake = iteration_tmp$mod_fake
    pis[k +1] =  iteration_tmp$pi_new
    if(bam[1]){
      mod_real$coefficients = rmvn(n = 1,mu = mod_real$coefficients,V = mod_real$Vp)
      names(mod_real$coefficients ) = names_mod_real
    } else {
      mod_real$coefficients  =   rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
      names(mod_real$coefficients ) = names_mod_real
    }
    
    if(bam[2]){
      mod_fake$coefficients = rmvn(n = 1,mu = mod_fake$coefficients,V = mod_fake$Vp)
      names(mod_fake$coefficients ) = names_mod_fake
    } else {
      mod_fake$coefficients  =   rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
      names(mod_fake$coefficients ) = names_mod_fake
    }
    betas_real[k+1,] = iteration_tmp$mod_real$coefficients
    betas_fake[k+1,]= iteration_tmp$mod_fake$coefficients
    cat("Iteration ", k," completed\n")
    
    
    cat("Pi =",    pis[k +1]  ,"\n")
    cat("Paramerers real = ", mod_real$coefficients, "\n")
    cat("Paramerers fake = ", mod_fake$coefficients,"\n")
  }
  # Calculate Standard Errors----
  model_real_list = list()
  model_fake_list = list()
  event_list = list()
  pis_final = c()
  for(k in 1:K_2) {
    iteration_tmp = mecm_iteration(comb_event_stream = comb_event_stream, 
                                   formula_real = formula_real, 
                                   formula_fake = formula_fake, 
                                   n_actors = n_actors, M = 1,
                                   change_all = change_all,
                                   mod_real = mod_real, 
                                   mod_fake =mod_fake, 
                                   end_time = end_time,
                                   beg_time =  beg_time,
                                   pi  =  pis[1],
                                   which_covs_real = which_covs_real, 
                                   which_covs_fake = which_covs_fake,
                                   seed = k + seed,
                                   ind_important = F,
                                   exo_cov = exo_cov, 
                                   exo_cat = exo_cat, 
                                   exo_dyad = exo_dyad,
                                   exo_cov_fake = exo_cov_fake, 
                                   exo_cat_fake = exo_cat_fake, 
                                   exo_dyad_fake = exo_dyad_fake,
                                   bam = bam, ...)
    mod_real = iteration_tmp$mod_real
    mod_fake = iteration_tmp$mod_fake
    model_real_list[[k]] = mod_real
    model_fake_list[[k]] = mod_fake
    event_list[[k]] = iteration_tmp$sampled_events
    pis_final[k] =  iteration_tmp$pi_new
    
    if(bam[1]){
      mod_real$coefficients = rmvn(n = 1,mu = mod_real$coefficients,V = mod_real$Vp)
    } else {
      mod_real$coefficients  =   rmvn(n = 1,mu = mod_real$coefficients,V = vcov(mod_real))
    }
    
    if(bam[2]){
      mod_fake$coefficients = rmvn(n = 1,mu = mod_fake$coefficients,V = mod_fake$Vp)
    } else {
      mod_fake$coefficients  =   rmvn(n = 1,mu = mod_fake$coefficients,V = vcov(mod_fake))
    }
    
    cat("Iteration ", k," completed\n")
  }
  # browser()
  # Application of Rubins Rule -----
  if(bam[1]){
    # browser()
    est_coefs = lapply(model_real_list,function(x) return(coef(x)))
    est_var = lapply(model_real_list,function(x) return(x$Vp))
    # B_bar = apply(est_coefs,MARGIN = 2,var)
    mean = 1/K_2*Reduce("+",est_coefs)
    B_bar =1/(K_2-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
    V_bar = 1/K_2*Reduce("+",est_var)
    final_var_real = V_bar + (1+1/K_2)*B_bar
    # cut info not needed 
    model_real_list = lapply(model_real_list, function(x){
      list(Vp = x$Vp,
           Ve = x$Ve , 
           Vc = x$Vc,
           coefficients = x$coefficients)      
    })
  } else {
    est_coefs = lapply(model_real_list,function(x) return(coef(x)))
    est_var = lapply(model_real_list,function(x) return(vcov(x)))
    # B_bar = apply(est_coefs,MARGIN = 2,var)
    mean = 1/K_2*Reduce("+",est_coefs)
    B_bar =1/(K_2-1)*Reduce(x = lapply(est_coefs, FUN = function(x){(x-mean)%*% t(x-mean)}), "+")
    V_bar = 1/K_2*Reduce("+",est_var)
    final_var_real = V_bar + (1+1/K_2)*B_bar
  }
  
  
  # Everything for the fake models
  if(bam[2]){
    est_coefs_fake = lapply(model_fake_list,function(x) return(coef(x)))
    est_var_fake = lapply(model_fake_list,function(x) return(x$Vp))
    mean_fake = 1/K_2*Reduce("+",est_coefs_fake)
    B_fake =1/(K_2-1)*Reduce(x = lapply(est_coefs_fake, FUN = function(x){(x-mean_fake)%*% t(x-mean_fake)}), "+")
    V_fake = 1/K_2*Reduce("+",est_var_fake)
    final_var_fake = V_fake + (1+1/K_2)*B_fake
    model_fake_list = lapply(model_fake_list, function(x){
      list(Vp = x$Vp,
           Ve = x$Ve , 
           Vc = x$Vc,
           coefficients = x$coefficients)      
    })
  } else {
    est_coefs_fake = lapply(model_fake_list,function(x) return(coef(x)))
    est_var_fake = lapply(model_fake_list,function(x) return(vcov(x)))
    mean_fake = 1/K_2*Reduce("+",est_coefs_fake)
    B_fake =1/(K_2-1)*Reduce(x = lapply(est_coefs_fake, FUN = function(x){(x-mean_fake)%*% t(x-mean_fake)}), "+")
    V_fake = 1/K_2*Reduce("+",est_var_fake)
    final_var_fake = V_fake + (1+1/K_2)*B_fake
  }
  
  
  res = list(D = K_2, 
             vcov_real = final_var_real, B_real = B_bar, W_real = V_bar, beta_final = mean,
             vcov_fake = final_var_fake, B_fake = B_fake, W_fake = V_fake, beta_fake_final = mean_fake,
             betas_real = betas_real, betas_fake = betas_fake,
             model_real_list = model_real_list,
             model_fake_list = model_fake_list,
             mod_fake = mod_fake, 
             event_list = event_list, 
             beg_models = beg_models, 
             pis_final = pis_final, 
             formula_real = formula_real, 
             formula_fake = formula_fake,
             call = call)
  class(res) = "da_ecrem_result"
  return(res)
}

e_step_compl_clustered_mgcv = function(m,comb_event_stream, pi, beg_time, end_time,
                                       n_actors, mod_real, mod_fake, 
                                       which_covs_real , change_all,
                                       which_covs_fake ,seed , 
                                       exo_cov, exo_cat,exo_dyad,
                                       exo_cov_fake, exo_cat_fake,exo_dyad_fake,ind_important,bam ) {
  
  exo_cov = c(exo_cov, exo_cov_fake)
  exo_cat =c(exo_cat, exo_cat_fake)
  exo_dyad = c(exo_dyad, exo_dyad_fake)
  exo_cov =  exo_cov[unique(c(names(exo_cov)))]
  exo_cat =  exo_cat[unique(c(names(exo_cat)))]
  exo_dyad =  exo_dyad[unique(c(names(exo_dyad)))]
  
  comb_event_stream = st_e_step_clustered_mgcv(comb_event_stream = comb_event_stream,
                                               n_actors = n_actors, 
                                               pi = pi, 
                                               beg_time = beg_time, 
                                               end_time = end_time, 
                                               mod_real = mod_real,
                                               mod_fake = mod_fake, 
                                               which_covs_real = which_covs_real,
                                               which_covs_fake=which_covs_fake ,
                                               seed = seed*m,
                                               exo_cov = exo_cov, 
                                               exo_dyad = exo_dyad,
                                               exo_cat = exo_cat,  
                                               bam = bam)
  event_real = comb_event_stream[number_real > 0]
  event_real$weight = event_real$number_real
  
  event_fake = comb_event_stream[number_fake > 0]
  event_fake$weight = event_fake$number_fake
  
  res = list(count_real = count_rem_undirected_clustered(event_data = event_real,exo_cov = exo_cov,
                                                         n_actors = n_actors,exo_cat = exo_cat,change_all = change_all,
                                                         exo_dyad = exo_dyad, add_var = 1, ind = m,ind_important = ind_important, 
                                                         end = end_time, start = beg_time), 
             count_fake = count_rem_undirected_clustered(event_data = event_fake,exo_cov = exo_cov,
                                                         n_actors = n_actors,exo_cat = exo_cat,change_all = change_all,
                                                         exo_dyad =exo_dyad, add_var = 0, ind = m,ind_important = ind_important, 
                                                         end = end_time, start = beg_time), 
             sampled_events =comb_event_stream )
  return(res)
}


st_e_step_clustered_mgcv = function(comb_event_stream,n_actors, 
                                    mod_real, mod_fake,
                                    pi, 
                                    beg_time, 
                                    end_time, 
                                    which_covs_real, 
                                    which_covs_fake, 
                                    seed = 123,exo_dyad, 
                                    exo_cov = exo_cov,
                                    exo_cat = exo_cat,bam) {
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
  data_present$repetition_bin = logical(length = nrow(data_present))
  data_present$repetition_bin = factor(data_present$repetition_bin)
  data_present$triangle_scaled = numeric(length = nrow(data_present))
  data_present$degree_sum_scaled = numeric(length = nrow(data_present))
  data_present$degree_abs_scaled = numeric(length = nrow(data_present))
  data_present$offset = 0
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  number_dyad = length(exo_dyad)
  if(number_dyad!=0){
    
    if(length(names(exo_dyad)) == 0){
      names_dyad = paste("cov_dyad", seq(1,number_cat),sep = "_")
    } else {
      names_dyad = paste("cov_dyad",names(exo_dyad),sep = "_")
    }
    for(i in 1:number_dyad){
      data_present[ , names_dyad[i] := as.numeric(exo_dyad[[i]][cbind(data_present$side_a, data_present$side_b)])] 
    }
    
  }
  if(number_cont!=0){
    if(length(names(exo_cov)) == 0){
      
      
      
      names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    } else {
      names_cont = paste("cov_cont",names(exo_cov),sep = "_")
    }
    
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    
    if(length(names(exo_cat)) == 0){
      names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    } else {
      names_cat = paste("cov_cat",names(exo_cat),sep = "_")
    }
    
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  n_possible = nrow(data_present)
  data_present$time = 0
  # Which covariates are actually used for the real and false events?
  matched_real = match(which_covs_real,names(data_present))
  matched_fake = match(which_covs_fake,names(data_present))
  # lambda =  numeric(length = nrow(data_present))
  which_hist = c()
  when_hist = c()
  comb_event_stream$times_rank = rank(comb_event_stream$times)
  comb_event_stream$number_real = 0
  comb_event_stream$number_fake = 0
  mean_probs = numeric(length = length(unique(comb_event_stream$times_rank)))
  last_event_time_real = beg_time 
  last_event_time_fake = beg_time 
  last_event_time = beg_time 
  
  n = 1
  for(i in unique(comb_event_stream$times_rank)){
    tmp_ind = comb_event_stream$times_rank == i
    side_as = comb_event_stream[tmp_ind]$side_a
    side_bs = comb_event_stream[tmp_ind]$side_b
    weight = comb_event_stream[tmp_ind]$weight
    ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
    # tmp_time = comb_event_stream[times_rank == i]$times[1]
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    data_present$offset = log(comb_event_stream[tmp_ind]$times[1] - last_event_time)
    last_event_time = comb_event_stream[tmp_ind]$times[1]
    if(bam[1]){
      # data_present$offset = log(comb_event_stream[tmp_ind]$times[1] - last_event_time_real)
      data_present$times = last_event_time_real
      lambda_real  = predict.gam(mod_real, type = "response", newdata = data_present)
    } else {
      # data_present$offset = log(comb_event_stream[tmp_ind]$times[1] - last_event_time_real)
      data_present$times = last_event_time_real
      lambda_real  = predict(mod_real, type = "response", newdata = data_present)
      
    }
    if(bam[2]){
      # data_present$offset = log(comb_event_stream[tmp_ind]$times[1] - last_event_time_fake)
      data_present$times = last_event_time_fake
      lambda_fake  = predict.gam(mod_fake, type = "response", newdata = data_present)
    } else {
      # data_present$offset = log(comb_event_stream[tmp_ind]$times[1] - last_event_time_fake)
      data_present$times = last_event_time_fake
      lambda_fake  = predict(mod_fake, type = "response", newdata = data_present)
    }
    
    # lambda_fake  = predict(mod_fake, type = "response", newdata = data_present)
    
    # Sample the event (fake or true)
    probs = (lambda_real[ind_wheres])/
      (lambda_real[ind_wheres] + lambda_fake[ind_wheres])
    mean_probs[n] =  mean(probs)
    n = n +1
    comb_event_stream[tmp_ind]$number_real = rbinom(n = length(side_as), size = weight,prob = probs)
    comb_event_stream[tmp_ind]$number_fake = weight - comb_event_stream[tmp_ind]$number_real
    # 
    # if(i>1000){
    #   browser()
    # }
    
    if(sum(comb_event_stream[tmp_ind]$number_fake) > 0){
      last_event_time_fake = last_event_time
      
      }
    if(sum(comb_event_stream[tmp_ind]$number_real) > 0){
      last_event_time_real = last_event_time
      side_as_change = side_as[comb_event_stream$number_real[tmp_ind]>0]
      side_bs_change = side_bs[comb_event_stream$number_real[tmp_ind]>0]
      
      tmp = update_covs_undirected_clustered(data_present =data_present,time = data_present$time,
                                             which_hist = which_hist,when_hist = when_hist, side_as =side_as, 
                                             side_bs = side_bs, n_actors =n_actors, 
                                             change = T,change_all = T) 
      which_hist = tmp$which_hist
      data_present = tmp$data_present
      data_present$repetition_bin = factor(data_present$repetition_bin)
      
      # data_present$time[data_present$change] = comb_event_stream[times_rank == i]$times[1]
    }
  }
  
  
  # event_data$times = cumsum(event_data$times)
  # close(pb)
  diff = Sys.time() -time_beg
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  # return(list(comb_event_stream, mean_probs))
  return(comb_event_stream)
} 



st_e_step_clustered = function(comb_event_stream,n_actors, beta_real, beta_fake,which_covs_real, which_covs_fake, seed = 123,exo_dyad, exo_cov = exo_cov,exo_cat = exo_cat) {
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
  data_present$repetition_bin = logical(length = nrow(data_present))
  data_present$triangle_scaled = numeric(length = nrow(data_present))
  data_present$degree_sum_scaled = numeric(length = nrow(data_present))
  data_present$degree_abs_scaled = numeric(length = nrow(data_present))
  
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  number_dyad = length(exo_dyad)
  if(number_dyad!=0){
    
    if(length(names(exo_dyad)) == 0){
      names_dyad = paste("cov_dyad", seq(1,number_cat),sep = "_")
    } else {
      names_dyad = paste("cov_dyad",names(exo_dyad),sep = "_")
    }
    for(i in 1:number_dyad){
      data_present[ , names_dyad[i] := as.numeric(exo_dyad[[i]][cbind(data_present$side_a, data_present$side_b)])] 
    }
    
  }
  if(number_cont!=0){
    if(length(names(exo_cov)) == 0){
      
      
      
      names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    } else {
      names_cont = paste("cov_cont",names(exo_cov),sep = "_")
    }
    
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    
    if(length(names(exo_cat)) == 0){
      names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    } else {
      names_cat = paste("cov_cat",names(exo_cat),sep = "_")
    }
    
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
  comb_event_stream$times_rank = rank(comb_event_stream$times)
  comb_event_stream$number_real = 0
  comb_event_stream$number_fake = 0
  for(i in unique(comb_event_stream$times_rank)){
    
    # setTxtProgressBar(pb, i)
    #Step a: Who will have an event? 
    # The first two columns are not used since the are only from and to 
    lambda_real  = as.vector(exp(beta_real %*% t(data_present[,matched_real, with = F])))
    lambda_fake  = as.vector(exp(beta_fake %*% t(data_present[,matched_fake, with = F])))
    weight = comb_event_stream[times_rank == i]$weight
    side_as = comb_event_stream[times_rank == i]$side_a
    side_bs = comb_event_stream[times_rank == i]$side_b
    ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
    # Sample the event (fake or true)
    probs = lambda_real[ind_wheres]/
      (lambda_real[ind_wheres] + lambda_fake[ind_wheres])
    
    comb_event_stream[times_rank == i]$number_real = rbinom(n = length(side_as), size = weight,prob = probs)
    comb_event_stream[times_rank == i]$number_fake = weight - comb_event_stream[times_rank == i]$number_real
    
    if(i>1000){
      browser()
    }
    
    if(sum(comb_event_stream[times_rank == i]$number_real) > 0){
      side_as_change = side_as[comb_event_stream[times_rank == i]$number_real>0]
      side_bs_change = side_bs[comb_event_stream[times_rank == i]$number_real>0]
      
      tmp = update_covs_undirected_clustered(data_present =data_present,which_hist = which_hist,side_as =side_as, side_bs = side_bs, n_actors =n_actors) 
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


e_step_compl_clustered = function(m,comb_event_stream,
                                  n_actors, beta_real, beta_fake, 
                                  which_covs_real , 
                                  which_covs_fake ,seed = k, 
                                  exo_cov, exo_cat,exo_dyad) {
  comb_event_stream = st_e_step_clustered(comb_event_stream = comb_event_stream,
                                          n_actors = n_actors, 
                                          beta_real = beta_real,
                                          beta_fake = beta_fake, 
                                          which_covs_real = which_covs_real,
                                          which_covs_fake=which_covs_fake ,
                                          seed = seed*m,
                                          exo_cov = exo_cov, 
                                          exo_dyad = exo_dyad,
                                          exo_cat = exo_cat)
  event_real = comb_event_stream[number_real > 0]
  event_real$weight = event_real$number_real
  
  event_fake = comb_event_stream[number_fake > 0]
  event_fake$weight = event_fake$number_fake
  
  return(rbind(count_rem_undirected_clustered(event_data = event_real,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat,exo_dyad = exo_dyad, add_var = 1, ind = m), 
               count_rem_undirected_clustered(event_data = event_fake,exo_cov = exo_cov,n_actors = n_actors,exo_cat = exo_cat,exo_dyad =exo_dyad, add_var = 0, ind = m)))
}



count_rem_undirected_clustered = function(event_data,exo_cov,
                                          exo_cat,exo_dyad, 
                                          n_actors, 
                                          add_var = 1,ind = 1,
                                          change_all = "all", end, start, 
                                          ind_important = F, get_all = F){
  time_beg = Sys.time()
  n = nrow(event_data)
  end_diff = F
  start_diff = F
  max_time = max(event_data$times)
  
  if(end > max(event_data$times)){
    end_diff = T
    max_time = end
  }
  if(start < min(event_data$times)){
    start_diff = T
  }
  
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$repetition_bin = logical(length = nrow(data_present))
  data_present$triangle_scaled = numeric(length = nrow(data_present))
  data_present$degree_sum_scaled = numeric(length = nrow(data_present))
  data_present$degree_abs_scaled = numeric(length = nrow(data_present))
  
  
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  number_dyad = length(exo_dyad)
  if(number_dyad!=0){
    
    if(length(names(exo_dyad)) == 0){
      names_dyad = paste("cov_dyad", seq(1,number_cat),sep = "_")
    } else {
      names_dyad = paste("cov_dyad",names(exo_dyad),sep = "_")
    }
    for(i in 1:number_dyad){
      data_present[ , names_dyad[i] := as.numeric(exo_dyad[[i]][cbind(data_present$side_a, data_present$side_b)])] 
    }
    
  }
  if(number_cont!=0){
    if(length(names(exo_cov)) == 0){
      names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    } else {
      names_cont = paste("cov_cont",names(exo_cov),sep = "_")
    }
    
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := abs(exo_cov[[i]][data_present$side_a] - exo_cov[[i]][data_present$side_b])] 
      
      # data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    
    if(length(names(exo_cat)) == 0){
      names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    } else {
      names_cat = paste("cov_cat",names(exo_cat),sep = "_")
    }
    
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  data_present$change = 1
  data_present$status = data_present$time = 0
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  
  
  
  which_hist = c()
  when_hist = c()
  event_data$times_rank = rank(event_data$times)
  data_ges = vector(mode = "list", length = length(unique(event_data$times_rank )))
  if(start_diff){
    data_ges[[1]] = data_present
    data_ges[[1]]$times = start
    k = 2
  } else {
    k = 1
  }
  
  
  for(i in unique(event_data$times_rank )){
    # setTxtProgressBar(pb, i)
    #Step a: Who had an event? 
    event_ind = event_data$times_rank == i
    side_as = event_data$side_a[event_ind] 
    side_bs = event_data$side_b[event_ind] 
    waiting_time = event_data$times[event_ind] 
    weights = event_data$weight[event_ind] 
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
    data_present$status[ind_wheres] = weights
    data_present$times = waiting_time[1]
    
    
    data_present$change = 0
    data_present$change[ind_wheres] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    
    tmp = update_covs_undirected_clustered(data_present = data_present, which_hist = which_hist,when_hist = when_hist,
                                           side_as = side_as, side_bs = side_bs, time = waiting_time[1],n_actors =n_actors, change = T,change_all = change_all)
    if(sum(tmp$data_present$status) != length(side_as)){
      browser()
    }
    data_present = tmp$data_present
    which_hist = tmp$which_hist
    
    data_ges[[k]] = data_present[change == 1]
    k = k +1
  }
  
  if(end_diff){
    data_present = tmp$data_present
    data_present$status = 0
    data_present$times = end
    data_ges[[k +1]] = data_present
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$side_a, concat_data$side_b,sep = "_")
  concat_data = concat_data[order(from_to)]
  tmp = concat_data[,.(time_end =  c(times[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$time_end = tmp$time_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[time_end != times]
  # close(pb)
  # Add info that we need to say if the counted events are true (1) or false (0) in add_var, to indicate the simulation run (ind), and to calculate the offset 
  concat_data$add_var = add_var
  concat_data$ind = ind
  concat_data$offset = log(concat_data$time_end-concat_data$times)
  
  if(get_all){
    return(concat_data)
  }
  
  if(ind_important){
    not_include = match(c("side_a", "side_b", "degree_side_a", "degree_side_b", "from_to"),names(concat_data))
  } else {
    not_include = match(c("side_a", "side_b", "degree_side_a", "degree_side_b", "from_to","ind"),names(concat_data))
  }
  # concat_data$weight = 1
  concat_data = concat_data[,.(weight = .N),by = c(names(concat_data)[-not_include])]
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}


count_rem_undirected_windowed = function(event_data,window, exo_cov, 
                                         exo_cat,exo_dyad, 
                                         n_actors, 
                                         add_var = 1,ind = 1,
                                         change_all = "all", end, start, 
                                         ind_important = F, get_all = F){
  time_beg = Sys.time()
  
  # Double all events where one event is the real event and the other one the exo change event 
  event_data$exo = 2
  max_time = max(event_data$times)
  
  event_data_exo = event_data
  event_data_exo$times = event_data_exo$times + window
  event_data_exo$exo = 1
  event_data = rbind(event_data, event_data_exo)
  event_data = event_data[order(times, exo)]
  
  # Get all times where only exo. changes occur 
  exo_times = unique(event_data[exo ==1]$times)
  endo_times = unique(event_data[exo ==2]$times)
  # Delete the times where the exo and end time coincide
  exo_times = exo_times[!exo_times%in% endo_times]
  # Get the first time of endo change before each exo_times
  endo_times = endo_times[findInterval(exo_times, endo_times)]
  
  n = nrow(event_data)
  end_diff = F
  start_diff = F
  
  if(end > max_time){
    end_diff = T
    event_data = rbind(event_data, event_data[length(event_data$times)])
    event_data$times[length(event_data$times)]= end 
    event_data$exo[length(event_data$times)] = 4
    
  }
  if(start < min(event_data$times)){
    start_diff = T
    event_data = rbind(event_data[1],event_data)
    event_data$times[1]= start 
    event_data$exo[1] = 4
  }
  if(end > max_time){
    max_time = end
  }
  
  event_data = event_data[times<=max_time]
  
  # pb <- txtProgressBar(min = 0, max =n, style = 3)
  
  data_present = data.table(matrix(combinations(x = n_actors),ncol = 2,byrow = T))
  names(data_present) = c("side_a","side_b")
  data_present$Intercept= 1
  data_present$degree_sum= numeric(length = nrow(data_present))
  data_present$degree_abs= numeric(length = nrow(data_present))
  data_present$triangle = numeric(length = nrow(data_present))
  data_present$repetition = numeric(length = nrow(data_present))
  data_present$repetition_bin = logical(length = nrow(data_present))
  data_present$degree_side_a= numeric(length = nrow(data_present))
  data_present$degree_side_b= numeric(length = nrow(data_present))
  # Windowed effects
  data_present$degree_sum_windowed= numeric(length = nrow(data_present))
  data_present$degree_abs_windowed= numeric(length = nrow(data_present))
  data_present$triangle_windowed = numeric(length = nrow(data_present))
  data_present$repetition_windowed = numeric(length = nrow(data_present))
  data_present$repetition_bin_windowed = logical(length = nrow(data_present))
  data_present$degree_side_a_windowed= numeric(length = nrow(data_present))
  data_present$degree_side_b_windowed= numeric(length = nrow(data_present))
  # Scaled effects 
  data_present$triangle_scaled = numeric(length = nrow(data_present))
  data_present$degree_sum_scaled = numeric(length = nrow(data_present))
  data_present$degree_abs_scaled = numeric(length = nrow(data_present))
  # Added to indicate which covariates changed (change), what time the change took place (time), and which event actually occurred (status)
  data_present$change = 1
  data_present$status = data_present$time = 0
  number_cont = length(exo_cov)
  number_cat = length(exo_cat)
  number_dyad = length(exo_dyad)
  if(number_dyad!=0){
    
    if(length(names(exo_dyad)) == 0){
      names_dyad = paste("cov_dyad", seq(1,number_cat),sep = "_")
    } else {
      names_dyad = paste("cov_dyad",names(exo_dyad),sep = "_")
    }
    for(i in 1:number_dyad){
      data_present[ , names_dyad[i] := as.numeric(exo_dyad[[i]][cbind(data_present$side_a, data_present$side_b)])] 
    }
    
  }
  if(number_cont!=0){
    if(length(names(exo_cov)) == 0){
      names_cont = paste("cov_cont", seq(1,number_cont),sep = "_")
    } else {
      names_cont = paste("cov_cont",names(exo_cov),sep = "_")
    }
    
    for(i in 1:number_cont){
      data_present[ , names_cont[i] := exo_cov[[i]][data_present$side_a] + exo_cov[[i]][data_present$side_b]] 
    }
    
  }
  
  if(number_cat!=0){
    
    if(length(names(exo_cat)) == 0){
      names_cat = paste("cov_cat", seq(1,number_cat),sep = "_")
    } else {
      names_cat = paste("cov_cat",names(exo_cat),sep = "_")
    }
    
    for(i in 1:number_cat){
      data_present[ , names_cat[i] := as.numeric(exo_cat[[i]][data_present$side_a] == exo_cat[[i]][data_present$side_b])] 
    }
    
  }
  
  
  # We only have to look at half the posible events since we look at undirected events 
  data_present = data_present[side_a<side_b] 
  # Include a vectors or events that already occurred and when they did 
  which_hist = c()
  when_hist = c()
  # Since the timings may be clustered we take care of each temporal rank individually 
  event_data$times_rank = rank(event_data$times)
  data_ges = vector(mode = "list", length = length(unique(event_data$times_rank )))
  data_ges[[1]] = data_present
  last_event_time_real = start
  k = 1
  for(i in unique(event_data$times_rank )){
    # setTxtProgressBar(pb, i)
    #Step a: Who had an event? 
    tmp_where = event_data$times_rank == i
    side_as = event_data$side_a[tmp_where] 
    side_bs = event_data$side_b[tmp_where] 
    waiting_time = event_data$times[tmp_where] 
    weights = event_data$weight[tmp_where] 
    exo_ind = event_data$exo[tmp_where] 
    
    # if(i>4000){
    #   browser()
    # }
    #Step b: Indicate that the event has happenend at least once
    data_present$status = 0
    ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
    ind_wheres_event = ind_wheres[exo_ind == 2]
    weights_event = weights[exo_ind == 2]
    data_present$status[ind_wheres_event] = weights_event
    data_present$times = waiting_time[1]
    
    data_present$change = 0
    data_present$change[ind_wheres] = 1
    
    #Step c: Search all edges that will occurr in the future and will also be affected 
    # Assuming that the actual event is i-> j 
    # 
    # Where will changes occurr? 
    # Outdegree Sender
    tmp = update_covs_undirected_windowed(data_present = data_present, which_hist = which_hist,when_hist = when_hist,
                                          side_as = side_as, side_bs = side_bs,exos = exo_ind, window = window, 
                                          time = waiting_time[1],n_actors =n_actors, change = T,change_all = change_all)
    if(sum(tmp$data_present$status) != length(side_as[exo_ind == 2])){
      browser()
    }
    data_present = tmp$data_present
    which_hist = tmp$which_hist
    when_hist = tmp$when_hist
    
    data_ges[[k]] = data_present[change == 1]
    k = k +1
  }
  concat_data = rbindlist(data_ges)
  concat_data$from_to = paste(concat_data$side_a, concat_data$side_b,sep = "_")
  concat_data = concat_data[order(from_to)]
  max_time = end
  tmp = concat_data[,.(times_end =  c(times[-1], max_time), 
                       status = c(status[-1],0)), 
                    by = from_to]
  concat_data$times_end = tmp$times_end
  concat_data$status = tmp$status
  
  concat_data = concat_data[times_end != times]
  # close(pb)
  # Add info that we need to say if the counted events are true (1) or false (0) in add_var, to indicate the simulation run (ind), and to calculate the offset 
  concat_data$add_var = add_var
  concat_data$ind = ind
  concat_data$offset = log(concat_data$times_end-concat_data$times)
  exo_times
  endo_times
  
  concat_data[times %in% exo_times]$times = endo_times[match(concat_data[times %in% exo_times]$times,exo_times)]
  
  if(get_all){
    return(concat_data)
  }
  if(ind_important){
    not_include = match(c("side_a", "side_b", "degree_side_a", "degree_side_b", "from_to", "times_end"),names(concat_data))
  } else {
    not_include = match(c("side_a", "side_b", "degree_side_a", "degree_side_b", "from_to", "ind"),names(concat_data))
    # not_include = match(c("ind"),names(concat_data))
    # not_include = match(c("side_a", "side_b", "degree_side_a", "degree_side_b", "from_to", "times_end","ind"),names(concat_data))
    
  }
  # concat_data$weight = 1
  concat_data = concat_data[,.(weight = .N),by = c(names(concat_data)[-not_include])]
  diff = Sys.time() -time_beg
  
  cat("Needed time for computation: ", diff,attributes(diff)$units,"\n")
  return(concat_data)
}

update_covs_undirected_windowed = function(data_present,which_hist,when_hist, window, side_as, side_bs, exos, time, n_actors, change = F, change_all = "all"){
  ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
  data_present$change = 0
  
  # This is only the case if the event is the beginning or end of data collection (then we dont need to change anything)
  if(4 %in% exos){
    data_present$status = 0
    data_present$change = 1
    return(list(data_present = data_present, which_hist = which_hist, when_hist = when_hist))
  }
  
  for(i in 1:length(ind_wheres)){
    if(exos[i] == 2){
      if(data_present$repetition[ind_wheres[i]] == 0) {
        #Change all windowed and normal covs if the event did not already happen
        tmp_deg_side_a_1 = find_side_a(a = side_as[i],x = n_actors)
        tmp_deg_side_a_2 =find_side_b(b = side_as[i],x = n_actors)
        tmp_deg_side_b_1 = find_side_a(a = side_bs[i],x = n_actors)
        tmp_deg_side_b_2 =find_side_b(b = side_bs[i],x = n_actors)
        
        data_present$degree_side_a[tmp_deg_side_a_1] = data_present$degree_side_a[tmp_deg_side_a_1] + 1
        data_present$degree_side_b[tmp_deg_side_a_2] = data_present$degree_side_b[tmp_deg_side_a_2] + 1
        data_present$degree_side_a[tmp_deg_side_b_1] = data_present$degree_side_a[tmp_deg_side_b_1] + 1
        data_present$degree_side_b[tmp_deg_side_b_2] = data_present$degree_side_b[tmp_deg_side_b_2] + 1
        
        data_present$degree_side_a_windowed[tmp_deg_side_a_1] = data_present$degree_side_a_windowed[tmp_deg_side_a_1] + 1
        data_present$degree_side_b_windowed[tmp_deg_side_a_2] = data_present$degree_side_b_windowed[tmp_deg_side_a_2] + 1
        data_present$degree_side_a_windowed[tmp_deg_side_b_1] = data_present$degree_side_a_windowed[tmp_deg_side_b_1] + 1
        data_present$degree_side_b_windowed[tmp_deg_side_b_2] = data_present$degree_side_b_windowed[tmp_deg_side_b_2] + 1
        
        
        
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
          side_b_tmp = rep(side_bs[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        if(length(triangle_tmp_2)!=0) {
          side_a_tmp = triangle_tmp_2
          side_b_tmp = rep(side_as[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
        # Transitivity
        # Increment 1 the transitivity if h was in sender j and may be in sender i 
        # or h was in receiver i and will be in receiver j 
        which_hist_window = which_hist[when_hist>(time[1] - window)]
        # What events already happened with side_a?
        ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist_window)
        ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
        
        # What events already happened with side_b?
        ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist_window)
        ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
        
        if(length(triangle_tmp_1)!=0) {
          side_a_tmp = triangle_tmp_1
          side_b_tmp = rep(side_bs[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        if(length(triangle_tmp_2)!=0) {
          side_a_tmp = triangle_tmp_2
          side_b_tmp = rep(side_as[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] = 
          data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] +1
        
        
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]  + 1
        data_present$repetition_windowed[ind_wheres[i]] = data_present$repetition_windowed[ind_wheres[i]]  + 1
        data_present$repetition_bin =  data_present$repetition>0
        data_present$repetition_bin_windowed =  data_present$repetition_windowed>0
        # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
        
        data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
        data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
        data_present$degree_sum_windowed =   data_present$degree_side_a_windowed +   data_present$degree_side_b_windowed
        data_present$degree_abs_windowed =  abs( data_present$degree_side_a_windowed -   data_present$degree_side_b_windowed)
        
        data_present$triangle_scaled = data_present$triangle/(n_actors -2)
        data_present$degree_sum_scaled = data_present$degree_sum/(2*(n_actors-1) )
        data_present$degree_abs_scaled = data_present$degree_abs/(n_actors-1)
        
        if(change) {
          change_ind = c(tmp_deg_side_a_1,
                         tmp_deg_side_a_2,
                         tmp_deg_side_b_1,
                         tmp_deg_side_b_2, 
                         triangle_tmp_1, 
                         triangle_tmp_2, 
                         ind_wheres[i])
          data_present$change[change_ind] = 1
        }
        if(change_all == "all") {
          data_present$change = 1
        }
        if(change_all == "week") {
          if(time %% 7 == 0){
            data_present$change = 1  
          }
        }
        
      } else if(data_present$repetition_windowed[ind_wheres[i]] == 0) {
        # if(time>600){
        #   browser()
        #   
        # }
        #Change all windowed covs and the rep for the normal if the event happened already once but not lately
        tmp_deg_side_a_1 = find_side_a(a = side_as[i],x = n_actors)
        tmp_deg_side_a_2 =find_side_b(b = side_as[i],x = n_actors)
        tmp_deg_side_b_1 = find_side_a(a = side_bs[i],x = n_actors)
        tmp_deg_side_b_2 =find_side_b(b = side_bs[i],x = n_actors)
        
        
        data_present$degree_side_a_windowed[tmp_deg_side_a_1] = data_present$degree_side_a_windowed[tmp_deg_side_a_1] + 1
        data_present$degree_side_b_windowed[tmp_deg_side_a_2] = data_present$degree_side_b_windowed[tmp_deg_side_a_2] + 1
        data_present$degree_side_a_windowed[tmp_deg_side_b_1] = data_present$degree_side_a_windowed[tmp_deg_side_b_1] + 1
        data_present$degree_side_b_windowed[tmp_deg_side_b_2] = data_present$degree_side_b_windowed[tmp_deg_side_b_2] + 1
        
        
        # Transitivity
        # Increment 1 the transitivity if h was in sender j and may be in sender i 
        # or h was in receiver i and will be in receiver j 
        which_hist_window = which_hist[when_hist>(time[1] - window)]
        # What events already happened with side_a?
        ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist_window)
        ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
        
        # What events already happened with side_b?
        ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist_window)
        ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
        
        if(length(triangle_tmp_1)!=0) {
          side_a_tmp = triangle_tmp_1
          side_b_tmp = rep(side_bs[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        if(length(triangle_tmp_2)!=0) {
          side_a_tmp = triangle_tmp_2
          side_b_tmp = rep(side_as[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] = 
          data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] +1
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]  + 1
        data_present$repetition_windowed[ind_wheres[i]] = data_present$repetition_windowed[ind_wheres[i]]  + 1
        data_present$repetition_bin_windowed =  data_present$repetition_windowed>0
        # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
        
        data_present$degree_sum_windowed =   data_present$degree_side_a_windowed +   data_present$degree_side_b_windowed
        data_present$degree_abs_windowed =  abs( data_present$degree_side_a_windowed -   data_present$degree_side_b_windowed)
        
        
        
        if(change) {
          change_ind = c(tmp_deg_side_a_1,
                         tmp_deg_side_a_2,
                         tmp_deg_side_b_1,
                         tmp_deg_side_b_2, 
                         triangle_tmp_1, 
                         triangle_tmp_2, 
                         ind_wheres[i])
          data_present$change[change_ind] = 1
        }
        if(change_all == "all") {
          data_present$change = 1
        }
        if(change_all == "week") {
          if(time %% 7 == 0){
            data_present$change = 1  
          }
        }
        
      }
      # What to do if the windowed repition is larger than 2 -> only change the rep for both
      else if(data_present$repetition_windowed[ind_wheres[i]]>0){
        #Change only windowed normals covs  for rep
        
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
        data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]  + 1
        data_present$repetition_windowed[ind_wheres[i]] = data_present$repetition_windowed[ind_wheres[i]]  + 1
        if(change) {
          data_present$change[ind_wheres[i]] = 1
        }
        if(change_all == "all") {
          data_present$change = 1
        }
        if(change_all == "week") {
          if(time %% 7 == 0){
            data_present$change = 1  
          }
        }
      }
      
      
    } else{
      if(data_present$repetition_windowed[ind_wheres[i]] == 1) {
        #What to do if the event to be retracted was up until now only observed once -> retract the influence  
        tmp_deg_side_a_1 = find_side_a(a = side_as[i],x = n_actors)
        tmp_deg_side_a_2 =find_side_b(b = side_as[i],x = n_actors)
        tmp_deg_side_b_1 = find_side_a(a = side_bs[i],x = n_actors)
        tmp_deg_side_b_2 =find_side_b(b = side_bs[i],x = n_actors)
        
        data_present$degree_side_a_windowed[tmp_deg_side_a_1] = data_present$degree_side_a_windowed[tmp_deg_side_a_1] - 1
        data_present$degree_side_b_windowed[tmp_deg_side_a_2] = data_present$degree_side_b_windowed[tmp_deg_side_a_2] - 1
        data_present$degree_side_a_windowed[tmp_deg_side_b_1] = data_present$degree_side_a_windowed[tmp_deg_side_b_1] - 1
        data_present$degree_side_b_windowed[tmp_deg_side_b_2] = data_present$degree_side_b_windowed[tmp_deg_side_b_2] - 1
        
        
        
        # Transitivity
        # Increment 1 the transitivity if h was in sender j and may be in sender i 
        # or h was in receiver i and will be in receiver j 
        which_hist_window = which_hist[when_hist>(time[1] - window)]
        # What events already happened with side_a?
        ind_with_side_a_1 = intersect(c(tmp_deg_side_a_1), which_hist_window)
        ind_with_side_a_2 = intersect(c(tmp_deg_side_a_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_1 = c(data_present$side_b[ind_with_side_a_1],data_present$side_a[ind_with_side_a_2])
        
        # What events already happened with side_b?
        ind_with_side_b_1 = intersect(c(tmp_deg_side_b_1), which_hist_window)
        ind_with_side_b_2 = intersect(c(tmp_deg_side_b_2), which_hist_window)
        # with whom were those events? 
        triangle_tmp_2 = c(data_present$side_b[ind_with_side_b_1],data_present$side_a[ind_with_side_b_2])
        
        if(length(triangle_tmp_1)!=0) {
          side_a_tmp = triangle_tmp_1
          side_b_tmp = rep(side_bs[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        
        if(length(triangle_tmp_2)!=0) {
          side_a_tmp = triangle_tmp_2
          side_b_tmp = rep(side_as[i], length(side_a_tmp))
          ind_change = side_a_tmp>side_b_tmp
          tmp = side_a_tmp[ind_change]
          side_a_tmp[ind_change] = side_b_tmp[ind_change]
          side_b_tmp[ind_change] = tmp
          triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
        }
        data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] = 
          data_present$triangle_windowed[c(triangle_tmp_1,triangle_tmp_2)] - 1
        data_present$triangle_windowed[data_present$triangle_windowed<0] = 0
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        data_present$repetition_windowed[ind_wheres[i]] = data_present$repetition_windowed[ind_wheres[i]]  - 1
        data_present$repetition_bin_windowed =  data_present$repetition_windowed>0
        # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
        
        data_present$degree_sum_windowed =   data_present$degree_side_a_windowed +   data_present$degree_side_b_windowed
        data_present$degree_abs_windowed =  abs( data_present$degree_side_a_windowed -   data_present$degree_side_b_windowed)
        
        if(change) {
          change_ind = c(tmp_deg_side_a_1,
                         tmp_deg_side_a_2,
                         tmp_deg_side_b_1,
                         tmp_deg_side_b_2, 
                         triangle_tmp_1, 
                         triangle_tmp_2, 
                         ind_wheres[i])
          data_present$change[change_ind] = 1
        }
        if(change_all == "all") {
          data_present$change = 1
        }
        if(change_all == "week") {
          if(time %% 7 == 0){
            data_present$change = 1  
          }
        }
        
      } else {
        
        # Repetition 
        # Increment the repetition by 1 for the tie that was observed  
        # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
        data_present$repetition_windowed[ind_wheres[i]] = data_present$repetition_windowed[ind_wheres[i]]  - 1
        data_present$repetition_bin_windowed =  data_present$repetition_windowed>0
        if(change) {
          data_present$change[ind_wheres[i]] = 1
        }
        if(change_all == "all") {
          data_present$change = 1
        }
        if(change_all == "week") {
          if(time %% 7 == 0){
            data_present$change = 1  
          }
        }
      }
    }
  }
  
  which_hist = c(which_hist,ind_wheres[exos == 2])
  when_hist = c(when_hist,rep(time[1], times = length(ind_wheres[exos == 2])))
  
  return(list(data_present = data_present, which_hist = which_hist, when_hist = when_hist))
}
update_covs_undirected_clustered = function(data_present,which_hist,when_hist, side_as, side_bs, time, n_actors, change = F, change_all = "all"){
  ind_wheres = findrows_undirected(as = side_as,bs = side_bs,x = n_actors)
  #Step c: Search all edges that will occurr in the future and will also be affected 
  # Assuming that the actual event is i-> j 
  # 
  # Where will changes occurr? 
  # Outdegree Sender
  data_present$change = 0
  for(i in 1:length(ind_wheres)){
    if(data_present$repetition[ind_wheres[i]] == 0) {
      #Check which(data_present$from == from)
      tmp_deg_side_a_1 = find_side_a(a = side_as[i],x = n_actors)
      tmp_deg_side_a_2 =find_side_b(b = side_as[i],x = n_actors)
      tmp_deg_side_b_1 = find_side_a(a = side_bs[i],x = n_actors)
      tmp_deg_side_b_2 =find_side_b(b = side_bs[i],x = n_actors)
      
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
        side_b_tmp = rep(side_bs[i], length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_1 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }
      
      if(length(triangle_tmp_2)!=0) {
        side_a_tmp = triangle_tmp_2
        side_b_tmp = rep(side_as[i], length(side_a_tmp))
        ind_change = side_a_tmp>side_b_tmp
        tmp = side_a_tmp[ind_change]
        side_a_tmp[ind_change] = side_b_tmp[ind_change]
        side_b_tmp[ind_change] = tmp
        triangle_tmp_2 = findrows_undirected(as = side_a_tmp, bs = side_b_tmp,x = n_actors)
      }
      
      data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] = data_present$triangle[c(triangle_tmp_1,triangle_tmp_2)] +1
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]  + 1
      data_present$repetition_bin =  data_present$repetition>0
      # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
      
      data_present$degree_sum =   data_present$degree_side_a +   data_present$degree_side_b
      data_present$degree_abs =  abs( data_present$degree_side_a -   data_present$degree_side_b)
      
      data_present$triangle_scaled = data_present$triangle/(n_actors -2)
      data_present$degree_sum_scaled = data_present$degree_sum/(2*(n_actors-1) )
      data_present$degree_abs_scaled = data_present$degree_abs/(n_actors-1)
      
      if(change) {
        change_ind = c(tmp_deg_side_a_1,
                       tmp_deg_side_a_2,
                       tmp_deg_side_b_1,
                       tmp_deg_side_b_2, 
                       triangle_tmp_1, 
                       triangle_tmp_2, 
                       ind_wheres[i])
        data_present$change[change_ind] = 1
      }
      if(change_all == "all") {
        data_present$change = 1
      }
      if(change_all == "week") {
        if(time %% 7 == 0){
          data_present$change = 1  
        }
      }
      
    } else {
      
      # Repetition 
      # Increment the repetition by 1 for the tie that was observed  
      # data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]+1
      data_present$repetition[ind_wheres[i]] = data_present$repetition[ind_wheres[i]]  + 1
      data_present$repetition_bin =  data_present$repetition>0
      
      if(change) {
        data_present$change[ind_wheres[i]] = 1
      }
      if(change_all == "all") {
        data_present$change = 1
      }
      if(change_all == "week") {
        if(time %% 7 == 0){
          data_present$change = 1  
        }
      }
    }
  }
  which_hist = c(which_hist,ind_wheres)
  when_hist = c(when_hist,rep(time, times = length(ind_wheres)))
  
  
  return(list(data_present = data_present, which_hist = which_hist, when_hist = when_hist))
}
preprocess = function(event_data, add_cov_data, year_events, year_actors, threshold_fatalities = 0, threshold_participation =5, 
                      exclude_undefined_militia= T,  exclude_rioters = T, exclude_rebels = T) {
  # only look at data from Syria and specific year(s)
  acled_syria = event_data[country == "Syria" & year %in% year_events & fatalities >= threshold_fatalities]
  acled_syria = acled_syria[!event_type %in% c("Strategic developments")]
  # Cut the unidentified actors 
  acled_syria = acled_syria[!grep("Unidentified",acled_syria$actor1)]
  acled_syria = acled_syria[!grep("Unidentified",acled_syria$actor2)]
  tbl_actor = table(c(acled_syria$actor1,acled_syria$actor2))
  # Cut all actors that were active in less than 5 events 
  tmp_actors = names(tbl_actor)[tbl_actor>threshold_participation]
  if(exclude_undefined_militia){
    tmp_actors = tmp_actors[-grep("Militia",tmp_actors)]
  }
  if(exclude_rioters){
    tmp_actors = tmp_actors[-grep("Rioters",tmp_actors)]
  }
  if(exclude_rebels){
    tmp_actors = tmp_actors[!tmp_actors %in% c("Kurdish Ethnic Militia (Syria)", "Rioters (Syria)", "Islamist Militia (Syria)", "Islamist Rebels (Syria)","Opposition Rebels (Syria)")]
  }
  # And cut all empty actors (unkown actors)
  tmp_actors = tmp_actors[- (nchar(tmp_actors) == 0)]
  # Now limit the event stream to actors that are in the vector tmp_actors
  acled_syria = acled_syria[actor1 %in% tmp_actors & actor2 %in% tmp_actors]
  # Make the names more consistent
  acled_syria$actor1[grep(pattern = "Military Forces of Syria",x = acled_syria$actor1)] = "Military Forces of Syria"
  acled_syria$actor2[grep(pattern = "Military Forces of Syria",x = acled_syria$actor2)] = "Military Forces of Syria"
  acled_syria$actor1[grep(pattern = "Police Forces of Syria",x = acled_syria$actor1)] = "Police Forces of Syria"
  acled_syria$actor2[grep(pattern = "Police Forces of Syria",x = acled_syria$actor2)] = "Police Forces of Syria"
  acled_syria$actor1[grep(pattern = "Wa Harredh al Moa'mineen",x = acled_syria$actor1)] = "Wa Harredh al Moa'mineen"
  acled_syria$actor2[grep(pattern = "Wa Harredh al Moa'mineen",x = acled_syria$actor2)] = "Wa Harredh al Moa'mineen"
  acled_syria$actor1[grep(pattern = "Military Forces of Turkey",x = acled_syria$actor1)] = "Military Forces of Turkey"
  acled_syria$actor2[grep(pattern = "Military Forces of Turkey",x = acled_syria$actor2)] = "Military Forces of Turkey"
  acled_syria$actor1[grep(pattern = "Military Forces of Iran",x = acled_syria$actor1)] = "Military Forces of Iran"
  acled_syria$actor2[grep(pattern = "Military Forces of Iran",x = acled_syria$actor2)] = "Military Forces of Iran"
  acled_syria$actor1[grep(pattern = "Military Forces of Jordan",x = acled_syria$actor1)] = "Military Forces of Jordan"
  acled_syria$actor2[grep(pattern = "Military Forces of Jordan",x = acled_syria$actor2)] = "Military Forces of Jordan"
  acled_syria$actor1[grep(pattern = "Popular Mobilization Forces",x = acled_syria$actor1)] = "Military Forces of Iraq: Popular Mobilization Forces"
  acled_syria$actor2[grep(pattern = "Popular Mobilization Forces",x = acled_syria$actor2)] = "Military Forces of Iraq: Popular Mobilization Forces"
  acled_syria$actor1[grep(pattern = "Military Forces of Iraq ",x = acled_syria$actor1)] = "Military Forces of Iraq"
  acled_syria$actor2[grep(pattern = "Military Forces of Iraq ",x = acled_syria$actor2)] = "Military Forces of Iraq"
  acled_syria$actor1[grep(pattern = "YPG",x = acled_syria$actor1)] = "YPG"
  acled_syria$actor2[grep(pattern = "YPG",x = acled_syria$actor2)] = "YPG"
  acled_syria$actor1[grep(pattern = "Military Forces of Russia",x = acled_syria$actor1)] = "Military Forces of Russia"
  acled_syria$actor2[grep(pattern = "Military Forces of Russia",x = acled_syria$actor2)] = "Military Forces of Russia"
  # These are the same actors 
  acled_syria$actor1[grep(pattern = "PYD: Kurdish Democratic Union Party",x = acled_syria$actor1)] = "YPG"
  acled_syria$actor2[grep(pattern = "PYD: Kurdish Democratic Union Party",x = acled_syria$actor2)] = "YPG"
  acled_syria$actor1[grep(pattern = "Sham Legion",x = acled_syria$actor1)] = "Al Sham Corps"
  acled_syria$actor2[grep(pattern = "Sham Legion",x = acled_syria$actor2)] = "Al Sham Corps"
  
  # Delete Self-Loops
  acled_syria = acled_syria[actor1 != actor2]
  
  tbl_actor = unique(c(acled_syria$actor1,acled_syria$actor2))
  acled_syria$actor1_num = match(acled_syria$actor1,tbl_actor)
  # Check 
  # table(tbl_actor[acled_syria$actor1_num] == acled_syria$actor1)
  acled_syria$actor2_num = match(acled_syria$actor2,tbl_actor)
  
  acled_syria$actor1_num_old = acled_syria$actor1_num
  acled_syria$actor2_num_old = acled_syria$actor2_num
  # The events are not directed but also not ordered accordingly (side_a<side_b)
  tmp = acled_syria$actor1_num<acled_syria$actor2_num
  acled_syria$actor1_order[tmp] = acled_syria$actor1_num[tmp]
  acled_syria$actor1_order[!tmp] = acled_syria$actor2_num[!tmp]
  acled_syria$actor2_order[tmp] = acled_syria$actor2_num[tmp]
  acled_syria$actor2_order[!tmp] = acled_syria$actor1_num[!tmp]
  
  acled_syria$actor1_order_name[tmp] = acled_syria$actor1[tmp]
  acled_syria$actor1_order_name[!tmp] = acled_syria$actor2[!tmp]
  acled_syria$actor2_order_name[tmp] = acled_syria$actor2[tmp]
  acled_syria$actor2_order_name[!tmp] = acled_syria$actor1[!tmp]
  
  acled_syria$actor1 = acled_syria$actor1_order_name
  acled_syria$actor2 = acled_syria$actor2_order_name
  acled_syria$actor1_num = acled_syria$actor1_order
  acled_syria$actor2_num = acled_syria$actor2_order
  
  acled_syria$event_date = as.Date(acled_syria$event_date,format='%d %B %Y')
  # acled_syria = acled_syria[acled_syria$fatalities>5]
  # We should probably aggregate events between the same actors on one day (each drone strike is one event but should the complete strike be one event?)
  # acled_syria$tmp_id = paste(acled_syria$actor1_num,acled_syria$actor2_num,acled_syria$event_date,acled_syria$event_type,sep = "_")
  acled_syria$tmp_id = paste(acled_syria$actor1_num,acled_syria$actor2_num,acled_syria$event_date,sep = "_")
  
  # acled_syria$tmp_id[acled_syria$actor1_num>acled_syria$actor2_num] = paste(acled_syria$actor2_num,acled_syria$actor1_num,
  #                                                                           acled_syria$event_date,acled_syria$event_type,sep = "_")[acled_syria$actor1_num>acled_syria$actor2_num]
  # 
  acled_syria_cut = acled_syria[,.(actor1_num = actor1_num[1], 
                                   actor2_num = actor2_num[1], 
                                   actor1 = actor1[1], 
                                   actor2 = actor2[1],
                                   inter1 = inter1[1], 
                                   inter2 = inter2[1],
                                   weight = 1,
                                   type = event_type[1],
                                   time_precision = min(time_precision),
                                   event_date = event_date[1],
                                   fatalities = sum(fatalities), 
                                   latitude = mean(latitude), 
                                   longitude = mean(longitude)), 
                                by = tmp_id]
  
  acled_syria_cut$dyad = paste(acled_syria_cut$actor1, acled_syria_cut$actor2,sep = "_")
  acled_syria_cut$conflict_with_civilians_protest = factor(grepl("Civilians",acled_syria_cut$dyad),labels = c("rebel","civil"))
  acled_syria_cut$conflict_with_civilians_protest[grepl("Protest",acled_syria_cut$dyad)] = "civil"
  split_events = split(acled_syria_cut,acled_syria_cut$conflict_with_civilians  )
  event_stream = split_events$rebel
  
  # event_stream$time = 1:length(event_stream$event_date)
  event_stream$time = as.numeric(event_stream$event_date - min(event_stream$event_date) +1)
  event_stream = event_stream[order(time)]
  
  event_data_obs = data.table(side_a =event_stream$actor1, side_b = event_stream$actor2, times = event_stream$time, date = event_stream$event_date, weight = event_stream$weight)
  unique_actors = unique(c(event_data_obs$side_a, event_data_obs$side_b))
  n_actors = length(unique_actors)
  
  # actors_rebels = data.table(id = unique(c(split_events$rebel$actor1_num,split_events$rebel$actor2_num)))
  # actors_rebels$name = c(split_events$rebel$actor1,split_events$rebel$actor2)[match(actors_rebels$id,c(split_events$rebel$actor1_num,split_events$rebel$actor2_num))]
  # actors_rebels$type = c(split_events$rebel$actor1_type,split_events$rebel$actor2_type)[match(actors_rebels$id,c(split_events$rebel$actor1_num,split_events$rebel$actor2_num))]
  # actors_rebels = actors_rebels[order(id)]
  # actors_rebels$id = 1:length(actors_rebels$id)
  
  
  # Only look at actors that were active in 2018
  # active_actors =unique(c(event_data_obs$side_a[year(event_data_obs$date) == 2017], event_data_obs$side_b[year(event_data_obs$date) == 2017]))
  
  active_actors = unique(c(event_data_obs$side_a[year(event_data_obs$date) %in% c(year_actors)], event_data_obs$side_b[year(event_data_obs$date)  %in% c(year_actors)]))
  add_cov_data[is.na(add_cov_data)] = 0
  add_cov_data = add_cov_data[add_cov_data$noidea == 0]
  add_cov_data = add_cov_data[add_cov_data$label %in% active_actors]
  add_cov_data$sponsor = ""
  add_cov_data$sponsor[add_cov_data$sponsor_quatar == 1] = "1"
  add_cov_data$sponsor[add_cov_data$sponsor_saudi == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_saudi == 1],"2")
  add_cov_data$sponsor[add_cov_data$sponsor_turkey == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_turkey == 1],"3")
  add_cov_data$sponsor[add_cov_data$sponsor_US == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_US == 1], "4")
  add_cov_data$sponsor[add_cov_data$sponsor_jordan == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_jordan == 1],"5")
  add_cov_data$sponsor[add_cov_data$sponsor_russia == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_russia == 1],"6")
  add_cov_data$sponsor[add_cov_data$sponsor_iran == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_iran == 1], "7")
  add_cov_data$sponsor[add_cov_data$sponsor_EU == 1] = paste(add_cov_data$sponsor[add_cov_data$sponsor_EU== 1], "8")
  add_cov_data$sponsor = trimws(add_cov_data$sponsor)
  
  overlap = function(x,y){
    unlist(lapply(1:length(x), function(g){length(intersect(unlist(strsplit(x[g]," ")), unlist(strsplit(y[g]," "))))}))
  }
  
  x = add_cov_data$sponsor
  names(x) = add_cov_data$label
  common_sponsor = outer( x,x, overlap )
  
  # active_actors[!active_actors %in% add_cov_data$label]
  
  #Coding for the etho_var 
  # 0 Sunni-moderate, 
  # 1 Sunni-jihadist/salafist, 
  # 2 Shia,
  # 3 Alevi, 
  # 4 Kurds, 
  # 5 mixed/not salient/external gov.
  # Coding of actor types 
  # 1 = "State Forces"
  # 2 = "Rebel Groups"
  # 3 = "Political Militias"
  # 4 = "Identity Militias"
  # 5 = "Rioters"
  # 6 = "Protesters"
  # 7 = "Civilians"
  # 8 = "External/Other Forces"
  add_cov_data$id = 1:length(add_cov_data$label)
  # add_cov_data$state_force = F
  # add_cov_data$state_force[grep(pattern = "Military Force",x = add_cov_data$label)] = T
  # add_cov_data$state_force[grep(pattern = "Police Forces of Syria",x = add_cov_data$label)] = T
  add_cov_data$state_force = add_cov_data$type == "State Forces"
  add_cov_data$state_force[grep(pattern = "Military Force",x = add_cov_data$label)] = T
  add_cov_data$state_force[grep(pattern = "Police Forces of Syria",x = add_cov_data$label)] = T
  
  # add_cov_data$state_force[grep(pattern = "Military Forces of Iraq: Popular Mobilization Forces",x = add_cov_data$label)]
  
  equal = function(x,y, z){
    unlist(lapply(1:length(x), function(g){x[g] + y[g] == z}))
  }
  
  x = add_cov_data$state_force
  names(x) = add_cov_data$label
  both_state_force = outer( x,x, equal, z = 2 )
  both_rebel = outer( x,x, equal, z = 0)
  
  
  exo_cat = list("ethn_relig" = add_cov_data$ethn_relig)
  
  exo_dyad = list("common_sponsor_1" = common_sponsor>0, 
                  "common_sponsor_2" = common_sponsor>1, 
                  "common_sponsor_3" = common_sponsor>2, 
                  "common_sponsor_number" = common_sponsor, 
                  "both_state_force" = both_state_force, 
                  "both_rebel" = both_rebel)
  
  event_data_obs$side_a = match(event_data_obs$side_a,add_cov_data$label)
  event_data_obs$side_b = match(event_data_obs$side_b,add_cov_data$label)
  event_data_obs = event_data_obs[!is.na(side_b)]
  event_data_obs = event_data_obs[!is.na(side_a)]
  
  event_data_obs$side_a_old = event_data_obs$side_a
  event_data_obs$side_b_old = event_data_obs$side_b
  event_data_obs$side_a_order = event_data_obs$side_a
  event_data_obs$side_b_order = event_data_obs$side_b
  # The events are not directed but also not ordered accordingly (side_a<side_b)
  tmp = event_data_obs$side_a<event_data_obs$side_b
  event_data_obs$side_a_order[tmp] = event_data_obs$side_a[tmp]
  event_data_obs$side_a_order[!tmp] = event_data_obs$side_b[!tmp]
  event_data_obs$side_b_order[tmp] = event_data_obs$side_b[tmp]
  event_data_obs$side_b_order[!tmp] = event_data_obs$side_a[!tmp]
  
  
  event_data_obs$side_a = event_data_obs$side_a_order
  event_data_obs$side_b = event_data_obs$side_b_order
  
  event_data_obs$label_a =add_cov_data$label[event_data_obs$side_a]
  event_data_obs$label_b = add_cov_data$label[event_data_obs$side_b]
  event_data_obs$eth_a =add_cov_data$ethn_relig[event_data_obs$side_a]
  event_data_obs$eth_b = add_cov_data$ethn_relig[event_data_obs$side_b]
  add_cov_data$active = lapply(add_cov_data$id,FUN =   function(x){
    event_data_obs$times[which(event_data_obs$side_a == x | event_data_obs$side_b == x)]
  })
  add_cov_data$first_event = lapply(add_cov_data$active,FUN =   function(x){min(x) })
  add_cov_data$last_event = lapply(add_cov_data$active,FUN =   function(x){max(x) })
  # event_data_obs$id = paste(event_data_obs$side_a, event_data_obs$side_b, event_data_obs$times)
  
  return(list(event_data_obs = event_data_obs, exo_cat = exo_cat, exo_dyad = exo_dyad, add_cov_data = add_cov_data))
}

findrow = function(a,b,x){
  if(b>a){
    return((a-1)*(x-1)+ b-1)
  } else{
    return((a-1)*(x-1)+ b)
  }
}

findrow_2 = function(b,a,x){
  if(b>a){
    return((a-1)*(x-1)+ b-1)
  } else{
    return((a-1)*(x-1)+ b)
  }
}

findrows_1 = function(bs,a,x){
  return(as.numeric(sapply(X = bs,a = a,x = x,FUN = findrow_2)))
}

findrows_2 = function(as,b,x){
  return(as.numeric(sapply(X = as,b = b,x = x,FUN = findrow)))
}

findfrom = function(a,x){
  return(((a-1)*(x-1)+1):((a-1)*(x-1)+x-1))
}

findto = function(b,x){
  if(b == 1){
    return(c(((b+1):x-1)*(x-1) +b))
  } else if(b == x){
    return(c(((1:(b-1) -1)*(x-1)+ (b-1))))
  } else {
    return(c(((1:(b-1) -1)*(x-1)+ (b-1)),((b+1):x-1)*(x-1) +b))
  }
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


findrows= function(event_stream,x = 84){
  a = event_stream$from
  b = event_stream$to
  which_1 = which(b>a)
  which_2 = which(b<=a)
  
  a_1 = a[which_1]
  b_1 = b[which_1]
  res_1 = c((a_1-1)*(x-1)+ b_1-1)
  a_2 = a[which_2]
  b_2 = b[which_2]
  res_2 = c((a_2-1)*(x-1)+ b_2)
  res = numeric(length = nrow(event_stream))
  res[which_1] = res_1
  res[which_2] = res_2
  return(res)
}
