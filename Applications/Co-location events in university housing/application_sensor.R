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
library(lubridate)

# Source scripts
source("../functions.R")
source("../../Simulation/simulation_function.R")

Rcpp::sourceCpp('../helper_fun.cpp')

# Load data 
event_data = fread("Data/Proximity.csv")
activities = fread("Data/Activities.csv")
activities = activities[survey.month == "2008.09"]
subjects = fread("Data/Subjects.csv")
political_interest = fread("Data/Politics.csv")
political_interest = political_interest[survey.month == "2008.09"]

friendship_ges = fread("Data/RelationshipsFromSurveys.csv")
friendship = friendship_ges[survey.date == "2008-09-09"]
friendship = friendship_ges[relationship == "CloseFriend"]
friendship = friendship[!id.A == id.B]
friendship = friendship[,c(1,2)]
friendship = friendship[!duplicated(friendship)]

friendship_fb_tag = friendship_ges[relationship == "PoliticalDiscussant"]
friendship_fb_tag = friendship_fb_tag[,c(1,2)]
friendship_fb_tag = friendship_fb_tag[!duplicated(friendship_fb_tag)]

event_data$date = ymd_hms(event_data$time)
# Only look at a subset of the data 
event_data = event_data[(date>ymd("2008-11-01")) & (date<ymd("2008-11-04"))]

# Cut all observations that are repeated between the same people and within 20 minutes 
event_data$id = paste(event_data$user.id, event_data$remote.user.id.if.known)
event_data$id_lag = c("N",event_data$id[-length(event_data$id)])
event_data$change_id = event_data$id != event_data$id_lag
event_data= event_data[order(id,date)]
event_data$date_lag = c(event_data$date[1],event_data$date[-length(event_data$id)])
event_data$time_diff = (event_data$date-event_data$date_lag)/(60*20)
event_data$change_id[event_data$time_diff>0 & event_data$time_diff<1]  = F
event_data = event_data[event_data$change_id,]
event_data = event_data[order(date)]
event_data$Timestamp = event_data$date - min(event_data$date) +1 
event_data$id = paste(event_data$user.id, event_data$remote.user.id.if.known, event_data$Timestamp)


id_to_user = data.table(user = unique(c(event_data$user.id, event_data$remote.user.id.if.known)), 
                        id = 1:(length(unique(c(event_data$user.id, event_data$remote.user.id.if.known)))))
tmp_match = !is.na(match(political_interest$user_id, id_to_user$user))
id_to_user$political_interest = " "
id_to_user$interested_in_politics = " "
id_to_user$liberal_or_conservative = " "


id_to_user$political_interest[match(political_interest$user_id, id_to_user$user)[tmp_match]] = political_interest$voting_for_today[tmp_match]
id_to_user$interested_in_politics[match(political_interest$user_id, id_to_user$user)[tmp_match]] = political_interest$interested_in_politics[tmp_match]
id_to_user$liberal_or_conservative[match(political_interest$user_id, id_to_user$user)[tmp_match]] = political_interest$liberal_or_conservative[tmp_match]
id_to_user$interested_in_politics[is.na(id_to_user$interested_in_politics)] = " "
id_to_user$liberal_or_conservative[is.na(id_to_user$liberal_or_conservative)] = " "
id_to_user$political_interest[grep("Obama",x = id_to_user$political_interest)] = "Obama"
id_to_user$political_interest[grep("John McCain",x = id_to_user$political_interest)] = "John McCain"
id_to_user$liberal_or_conservative_num = factor(id_to_user$liberal_or_conservative, levels = c("Extremely conservative", "Conservative","Slightly conservative",
                                                                                               "Moderate middle of the road","Slightly liberal","Liberal","Extremely liberal")) 
tmp_match = !is.na(match(subjects$user_id, id_to_user$user))
id_to_user$year_school[match(subjects$user_id, id_to_user$user)[tmp_match]] = subjects$year_school[tmp_match]
id_to_user$floor[match(subjects$user_id, id_to_user$user)[tmp_match]] = subjects$floor[tmp_match]


friendship = friendship[(id.A %in% id_to_user$user) & (id.B %in% id_to_user$user)]
friendship$id.A = match(friendship$id.A, id_to_user$user)
friendship$id.B = match(friendship$id.B, id_to_user$user)
friendship_alt = friendship
tmp = friendship_alt$id.A 
friendship_alt$id.A =  friendship_alt$id.B
friendship_alt$id.B = tmp
friendship = rbind(friendship_alt, friendship)
friendship_mat = matrix(F, nrow = nrow(id_to_user), ncol = nrow(id_to_user))
friendship_mat[cbind(friendship$id.A, friendship$id.B)] = T

friendship_fb_tag = friendship_fb_tag[(id.A %in% id_to_user$user) & (id.B %in% id_to_user$user)]
friendship_fb_tag$id.A = match(friendship_fb_tag$id.A, id_to_user$user)
friendship_fb_tag$id.B = match(friendship_fb_tag$id.B, id_to_user$user)
friendship_alt = friendship_fb_tag
tmp = friendship_alt$id.A 
friendship_alt$id.A =  friendship_alt$id.B
friendship_alt$id.B = tmp
friendship_fb_tag = rbind(friendship_alt, friendship_fb_tag)
fb_tag_mat = matrix(F, nrow = nrow(id_to_user), ncol = nrow(id_to_user))
fb_tag_mat[cbind(friendship_fb_tag$id.A, friendship_fb_tag$id.B)] = T

event_data$user.id = match(event_data$user.id, id_to_user$user)
event_data$remote.user.id.if.known = match(event_data$remote.user.id.if.known, id_to_user$user)
event_data$side_a = event_data$user.id
event_data$side_b = event_data$remote.user.id.if.known
event_data$side_a_old = event_data$side_a
event_data$side_b_old = event_data$side_b
event_data$side_a_order = event_data$side_a
event_data$side_b_order = event_data$side_b
# The events are not directed but also not ordered accordingly (side_a<side_b)
tmp = event_data$side_a<event_data$side_b
event_data$side_a_order[tmp] = event_data$side_a[tmp]
event_data$side_a_order[!tmp] = event_data$side_b[!tmp]
event_data$side_b_order[tmp] = event_data$side_b[tmp]
event_data$side_b_order[!tmp] = event_data$side_a[!tmp]

event_data$side_a = event_data$side_a_order
event_data$side_b = event_data$side_b_order
event_data = event_data[side_a!=side_b]

event_data$times = as.numeric(event_data$Timestamp)/(60)
event_data$weight = 1
n_actors = length(unique(c(event_data$side_a, event_data$side_b)))
id_to_user$interested_in_politics_num = as.numeric(factor(id_to_user$interested_in_politics))


# Set up model 
beta_fake_gt = c("Intercept" = T, "degree_sum" = F,"degree_abs" = F,
                 "triangle" = F, "repetition" = F, "repetition_bin" = F, "triangle_scaled"= F, 
                 "degree_sum_scaled" = F, "degree_abs_scaled"= F)

beta_real_gt = c("Intercept" = T, "degree_sum" = F,"degree_abs" = T,
                 "triangle" = T, "repetition" = T, "repetition_bin" = T, "triangle_scaled" =F, 
                 "degree_sum_scaled" = F, "degree_abs_scaled" = F)

result_sensor_data_remse =  da_estimation_clustered_mgcv(K = 10, K_2 = 10,only_beg = F,
                                                         beta_real_gt = beta_real_gt, rand_prop_real = 0.5,
                                                         separable = F,
                                                         beta_fake_gt = beta_fake_gt,
                                                         comb_event_stream = event_data,
                                                         beg_formula_real = "status ~  s(times, bs = 'ps')  + offset(offset) +", 
                                                         beg_formula_fake= "status ~offset(offset)", 
                                                         n_actors = n_actors, 
                                                         seed = 1123,eps = 0.01,
                                                         exo_cov = list(),
                                                         change_all = "all",
                                                         exo_dyad =list("friendship" =friendship_mat, 
                                                                        "interested_in_politics" = outer(id_to_user$interested_in_politics_num, id_to_user$interested_in_politics_num, function(x,y){abs(x-y)}) ), 
                                                         exo_cat  = list("presidential_pref" = id_to_user$political_interest, 
                                                                         "floor" = id_to_user$floor,
                                                                         "year" = id_to_user$year_school),
                                                         exo_cov_fake = list(),
                                                         exo_dyad_fake = list(), 
                                                         exo_cat_fake  = list(),
                                                         build = F,
                                                         build_time = min(event_data_obs$times[year(event_data_obs$date) == 2018]))


result_sensor_data_rem = da_estimation_clustered_mgcv(K = 30, K_2 = 20,only_beg = T,
                                       beta_real_gt = beta_real_gt, rand_prop_real =1,
                                       separable = F,beg_time = 0,end_time = 4318,
                                       beta_fake_gt = beta_fake_gt,
                                       comb_event_stream = event_data,
                                       beg_formula_real = "status ~  s(times, bs = 'ps')  + offset(offset) +", 
                                       beg_formula_fake= "status ~offset(offset)", 
                                       n_actors = n_actors, 
                                       seed = 1123,eps = 0.01,
                                       exo_cov = list(),
                                       change_all = "all",
                                       exo_dyad =list("friendship" =friendship_mat, 
                                                      "interested_in_politics" = outer(id_to_user$interested_in_politics_num, id_to_user$interested_in_politics_num, function(x,y){abs(x-y)}) ), 
                                       exo_cat  = list("presidential_pref" = id_to_user$political_interest, 
                                                       "floor" = id_to_user$floor,
                                                       "year" = id_to_user$year_school),
                                       exo_cov_fake = list(),
                                       exo_dyad_fake = list(), 
                                       exo_cat_fake  = list(),
                                       build = F,
                                       build_time = min(event_data_obs$times[year(event_data_obs$date) == 2018]))


part_ecrem = fun_mi_ecrem_result(result_sensor_data_remse)
fill = part_ecrem
n = nrow(part_ecrem)
fill[,] = ""
part_ecrem = rbind(part_ecrem, fill)[rep(1:n, each= 2) + rep(c(0,n), times = n)]
part_ecrem$coef[seq(from = 2, to = nrow(part_ecrem), by = 2)] = part_ecrem$CI[seq(from = 1, to = nrow(part_ecrem)-1, by = 2)]
part_ecrem$CI = NULL
part_ecrem$se = NULL

part_rem =  fun_rem_result(result_sensor_data_rem)[,c(2:5)]
fill = part_rem
n = nrow(part_rem)
fill[,] = ""
part_rem = rbind(part_rem, fill)[rep(1:n, each= 2) + rep(c(0,n), times = n)]
part_rem$coef[seq(from = 2, to = nrow(part_rem), by = 2)] = part_rem$CI[seq(from = 1, to = nrow(part_rem)-1, by = 2)]
part_rem$CI = NULL
part_rem$se = NULL

coefs = cbind(part_ecrem, part_rem)

names(coefs) = c("name",rep(c("Coef./CI", "Z Val."),times = 2))
tmp_name = coefs$name
tmp_name[tmp_name == "Repetition Bintrue"] = "First Repetition"
coefs$name = NULL
latex(object = coefs,label = "tbl:res_sensor",rowname  =tmp_name ,
      cgroup = c("EcREM", "REM"),first.hline.double= F,rowlabel.just = "l",
      n.cgroup = c(2,2),cgroupTexCmd = " ",rgroupTexCmd = " ",
      caption = "Co-location Events in University Housing: Estimated coefficients with confidence bands noted in bracets in the first column,
      while the Z values are given in the second column. 
      The results of the EcREM are given in the first two columns, 
      while the coefficients of the REM are depicted in the last two columns. ")
