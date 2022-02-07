setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../Applications/functions_count.R")
source("simulation_function.R")
Rcpp::sourceCpp('../Applications/helper_fun.cpp')
library(data.table)
library(parallel)
library(fastglm)
library(mgcv)
library(stringr)


# Simulate the covariates Situation 1 ----

set.seed(1343)
# set.seed(123)

n_actors = 40
exo_cov = rnorm(n = n_actors)
exo_cat = sample(size  = n_actors,x = c(1,2,3,4,5,6,7),replace = T)

beta_real_gt = c("Intercept" = -5, "degree_sum" = 0,"degree_abs" = 0.2,
                 "triangle" = 0.1, "repetition" =-0.5,
                 "cov_cont" =2,"cov_cat" = -2)
beta_fake_gt = c("Intercept" = -2.5, "degree_sum" = 0,"degree_abs" = 0,
                 "triangle" = 0, "repetition" = 0,
                 "cov_cont" =0,"cov_cat" = 0)

gc(reset = T,full = T)
no_cores = 10
cl<- makeCluster(no_cores, type="FORK")
set.seed(123)
times = parLapply(cl = cl,X = 1:(1000),fun = simulation_check_number, 
                   n_actors = 40,K = 100,
                   n = 500, 
                   exo_cov = exo_cov, 
                   exo_cat = exo_cat,
                   beta_real_gt = beta_real_gt, 
                   beta_fake_gt =  beta_fake_gt)
save(times, file = "Results/times_simulation_1.RData")
result = parLapply(cl = cl,X = 1:(1000),fun = simulation_complete_final, 
                   n_actors = 40,K = 100,K_2 = 40,
                   n = 500, 
                   exo_cov = exo_cov, 
                   exo_cat = exo_cat,
                   beta_real_gt = beta_real_gt, 
                   beta_fake_gt =  beta_fake_gt)
stopCluster(cl)
save(result, file = "Results/final_result_simulation_1.RData")


# Simulate the covariates Situation 2 ----

beta_real_gt = c("Intercept" = -5, "degree_sum" = 0,"degree_abs" = 0.2,
                 "triangle" = 0.1, "repetition" =-0.5,
                 "cov_cont" =2,"cov_cat" = -3)
beta_fake_gt = c("Intercept" = -15, "degree_sum" = 0,"degree_abs" = 0,
                 "triangle" = 0, "repetition" = 0,
                 "cov_cont" =0,"cov_cat" = 0)
no_cores = 10
cl<- makeCluster(no_cores, type="FORK")
set.seed(123)
result = parLapply(cl = cl,X = 1:(1000)*13,fun = simulation_complete_final, 
                   n_actors = 40,K = 50,K_2 = 5,
                   n = 500, 
                   exo_cov = exo_cov, 
                   exo_cat = exo_cat,
                   beta_real_gt = beta_real_gt, 
                   beta_fake_gt =  beta_fake_gt)
times = parLapply(cl = cl,X = 1:(1000)*13,fun = simulation_check_number, 
                  n_actors = 40,K = 50,
                  n = 500, 
                  exo_cov = exo_cov, 
                  exo_cat = exo_cat,
                  beta_real_gt = beta_real_gt, 
                  beta_fake_gt =  beta_fake_gt)
save(times, file = "Results/times_simulation_2.RData")
stopCluster(cl)
save(result, file = "Results/final_result_simulation_2.RData")

# Get information
load("times_simulation_1.RData")
proportion_real_1 = unlist(lapply(times, function(x){x[[1]]}))
load("times_simulation_2.RData")
proportion_real_2 = unlist(lapply(times, function(x){x[[1]]}))

load("final_result_simulation_1.RData")
result_1 = get_info_simulation(result = result)
load("final_result_simulation_2.RData")
result_2 = get_info_simulation(result = result)


