# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

# Load packages
library(data.table)
library(texreg)
library(parallel)
library(fastglm)
library(Matrix)
library(mgcv)
library(Hmisc)
library(stringr)

# Source scripts
source("../../Simulation/simulation_function.R")
source("../functions.R")
Rcpp::sourceCpp('../helper_fun.cpp')

# Load data 
event_data_3 <- fread("Data/1997-01-01-2021-03-31-Afghanistan-Iraq-Nigeria-Syria.csv") #acled
add_cov_data <- fread("Data/ACLEDgroups_covariates_update.csv") #acled
acronyms <- fread("Data/ACLEDgroups_covariates.csv") #acled
add_cov_data$acronym  = acronyms$acronym
# Preprocess data
data_processed_2018_build = preprocess(event_data = event_data_3, year_events =  2017:2018,
                                       year_actors = 2017:2018, add_cov_data = add_cov_data,
                                       threshold_fatalities  = 0,threshold_participation = 5,
                                       exclude_undefined_militia = T,exclude_rioters = T)

event_data_obs = data_processed_2018_build$event_data_obs
exo_cat = data_processed_2018_build$exo_cat
exo_dyad = data_processed_2018_build$exo_dyad

beta_fake_gt = c("Intercept" = T, "degree_sum" = F,"degree_abs" = F,
                 "triangle" = F, "repetition" = F, "repetition_bin" = F, "triangle_scaled"= F, 
                 "degree_sum_scaled" = F, "degree_abs_scaled"= F)
n_actors = nrow(exo_dyad$common_sponsor_1)
change_all = "all"
beta_real_gt = c("Intercept" = T, "degree_sum" =F,"degree_abs" = T,
                 "triangle" =T, "repetition" = T, "repetition_bin" = T, 
                 "degree_sum_windowed" = F,"degree_abs_windowed" = F,
                 "triangle_windowed" = F, "repetition_windowed" = F, "repetition_bin_windowed" = F, 
                 "triangle_scaled" = F,  "degree_sum_scaled" = F, "degree_abs_scaled" = F)
beg_formula_real = "status ~   s(times, bs = 'ps') + offset(offset) +"


result_conflict_data_remse =  da_estimation_clustered_mgcv(K = 1, K_2 = 3,
                                                           only_beg =F,full_bayes = F,
                                          beta_real_gt = beta_real_gt, 
                                          separable = F,rand_prop_real = 0.5,
                                          beta_fake_gt = beta_fake_gt,
                                          comb_event_stream =event_data_obs,
                                          beg_formula_real = beg_formula_real, 
                                          beg_formula_fake= "status ~  offset(offset) ", 
                                          n_actors = n_actors, 
                                          beg_time = 0,end_time = 730,
                                          seed = 1,eps = 0.01,
                                          exo_cov = list(),
                                          change_all = "all",
                                          exo_dyad = exo_dyad[c(1,5,6)], 
                                          exo_cat  = exo_cat,
                                          exo_cov_fake = list(),
                                          exo_dyad_fake = list(), 
                                          exo_cat_fake  = list(),
                                          build = F,
                                          build_time = 182)
round(mean(result_conflict_data_remse$pis_final), 3)

result_conflict_data_rem =  da_estimation_clustered_mgcv(K = 20, K_2 = 5,only_beg =T,full_bayes = F, 
                                          beta_real_gt = beta_real_gt, 
                                          separable = F,rand_prop_real = 0.5,
                                          beta_fake_gt = beta_fake_gt,
                                          comb_event_stream =event_data_obs,
                                          beg_formula_real = "status ~   s(times, bs = 'ps') + offset(offset) +", 
                                          beg_formula_fake= "status ~ offset(offset) ", 
                                          n_actors = n_actors, 
                                          beg_time = 0,end_time = 730,
                                          seed = 123,eps = 0.01,
                                          exo_cov = list(),
                                          change_all = "all",
                                          exo_dyad = exo_dyad[c(1,5,6)], 
                                          exo_cat  = exo_cat,
                                          exo_cov_fake = list(),
                                          exo_dyad_fake = list(), 
                                          exo_cat_fake  = list(),
                                          build = F,
                                          build_time = 182)

part_ecrem = fun_mi_ecrem_result(result_conflict_data_remse)
fill = part_ecrem
n = nrow(part_ecrem)
fill[,] = ""
part_ecrem = rbind(part_ecrem, fill)[rep(1:n, each= 2) + rep(c(0,n), times = n)]
part_ecrem$coef[seq(from = 2, to = nrow(part_ecrem), by = 2)] = part_ecrem$CI[seq(from = 1, to = nrow(part_ecrem)-1, by = 2)]
part_ecrem$CI = NULL
part_rem =  fun_rem_result(result_conflict_data_rem)[,c(2:5)]
fill = part_rem
n = nrow(part_rem)
fill[,] = ""
part_rem = rbind(part_rem, fill)[rep(1:n, each= 2) + rep(c(0,n), times = n)]
part_rem$coef[seq(from = 2, to = nrow(part_rem), by = 2)] = part_rem$CI[seq(from = 1, to = nrow(part_rem)-1, by = 2)]
part_rem$CI = NULL
part_ecrem$se = part_rem$se = NULL
coefs = cbind(part_ecrem, part_rem)
names(coefs) = c("name",rep(c("Coef./CI","Z Val."),times = 2))
rownames(coefs) = coefs$name
tmp_name = coefs$name
coefs$name = NULL

latex(coefs,label = "tbl:res_sensor",rowname  =tmp_name ,
      cgroup = c("EcREM", "REM"),
      n.cgroup = c(2,2),cgroupTexCmd = " ",rgroupTexCmd = " ",
      caption = "Estimated coefficients with confidence bands noted in bracets. 
      The results of the EcREM are given in the first column, 
      while the coefficients of the REM are depicted in the second column. ")


