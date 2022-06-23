# creating dataframes of weekly hospitalisations

# we do this in batches of 100 parameter sets, otherwise R runs out of memory
sampled_pars = c("accept"=unique(all_dynamics_accepted$par_id),"reject"=unique(all_dynamics_rejected$par_id))
start_stop_table = cbind(seq(1,length(sampled_pars),by=100),
                         c(seq(1,floor(length(sampled_pars)/100)*100-1,by=100)+99,length(sampled_pars)))

for (k_par in 1:nrow(start_stop_table)){
  
  if (k_par==1) { simul_hosp_rate_weekly_all_broad_agegroups=list()
                  simul_hosp_rate_weekly_under5_over65=list() }
  # selected params
  sel_pars <- sampled_pars[start_stop_table[k_par,1]:start_stop_table[k_par,2]]
  
  if (k_par<21) {
    work_df <- all_dynamics_accepted %>% ungroup %>% filter(par_id %in% sel_pars )
  } else {
    work_df <- all_dynamics_rejected %>% ungroup %>% filter(par_id %in% sel_pars )
  }
  
  print(k_par) # round(100*k_par/length(sampled_pars),1)
  
  # tic()
  processed_df <- work_df %>% filter( (agegroup<=3 | agegroup==5)) %>%
    mutate(hosp=value*hosp_probabilities_broad_age$prob_hosp_per_infection[agegroup]) %>%
    mutate(year_week=paste0(isoyear(date),"-", isoweek(date)),
           year_week=factor(year_week,levels=unique(year_week))) %>%
    group_by(year_week,agegroup,par_id) %>%
    summarise(simul_hosp_sum_full=sum(hosp),date=min(date)) %>%
    relocate(date,.after=year_week) %>% relocate(par_id,.after=date)
  # all age groups separately
  simul_hosp_rate_weekly_all_broad_agegroups[[k_par]] <- processed_df
  
  # into two age groups only
  simul_hosp_rate_weekly_under5_over65[[k_par]] <- processed_df %>% ungroup() %>%
    mutate(broad_age=ifelse(agegroup<=3,"<5y",">65y")) %>%
    group_by(year_week,par_id,broad_age) %>% 
    summarise(simul_hosp_sum_full=sum(simul_hosp_sum_full),date=unique(date) ) %>% # ,accepted=unique(accepted)
    mutate(simul_hosp_rate_100k=ifelse(grepl("65y",broad_age),
                                       1e5*simul_hosp_sum_full/rsv_age_groups$stationary_popul[11],
                                       1e5*simul_hosp_sum_full/sum(rsv_age_groups$stationary_popul[1:7]) ) )
  
  # 
  gc()
  # toc()
  print(paste0(k_par," done"))
}

# list to dataframe
simul_hosp_rate_weekly_all_broad_agegroups <- bind_rows(simul_hosp_rate_weekly_all_broad_agegroups)
simul_hosp_rate_weekly_under5_over65 <- bind_rows(simul_hosp_rate_weekly_under5_over65)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# FOR SIMULATIONs with gradual relaxation of NPI from March 2021

# creating data frames of weekly hospitalisations
# # we do this in batches of 100 parameter sets, otherwise R runs out of memory
# sampled_pars = c("accept"=unique(all_dynamics_accepted$par_id),"reject"=unique(all_dynamics_rejected$par_id))
# start_stop_table = cbind(seq(1,length(sampled_pars),by=100),
#                          c(seq(1,floor(length(sampled_pars)/100)*100-1,by=100)+99,length(sampled_pars)))
# 
# for (k_par in 1:nrow(start_stop_table)){
# 
#   if (k_par==1) { simul_hosp_rate_weekly_all_broad_agegroups_grad_relax=list();
#   simul_hosp_rate_weekly_under5_over65_grad_relax=list() }
#   # selected params
#   sel_pars <- sampled_pars[start_stop_table[k_par,1]:start_stop_table[k_par,2]]
# 
#   print(k_par) # round(100*k_par/length(sampled_pars),1)
# 
#   # tic()
#   processed_df <- all_dynamics_accepted_2e3_gradual_relax %>% ungroup %>% filter(par_id %in% sel_pars) %>%
#     mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
# #     filter(date>=as.Date("2018-09-25") & date<as.Date("2024-10-01")) %>%
#      filter( (agegroup<=3 | agegroup==5)) %>%
#     mutate(hosp=value*hosp_probabilities_broad_age$prob_hosp_per_infection[agegroup]) %>%
#     mutate(year_week=paste0(isoyear(date),"-", isoweek(date)),
#            year_week=factor(year_week,levels=unique(year_week))) %>%
#     group_by(year_week,agegroup,par_id) %>%
#     summarise(simul_hosp_sum_full=sum(hosp),date=min(date)) %>%
#     relocate(date,.after=year_week) %>% relocate(par_id,.after=date)
#   # all age groups separately
#   simul_hosp_rate_weekly_all_broad_agegroups_grad_relax[[k_par]] <- processed_df
# 
#   # into two age groups only
#   simul_hosp_rate_weekly_under5_over65_grad_relax[[k_par]] <- processed_df %>% ungroup() %>%
#     mutate(broad_age=ifelse(agegroup<=3,"<5y",">65y")) %>%
#     group_by(year_week,par_id,broad_age) %>%
#     summarise(simul_hosp_sum_full=sum(simul_hosp_sum_full),date=unique(date) ) %>% # ,accepted=unique(accepted)
#     mutate(simul_hosp_rate_100k=ifelse(grepl("65y",broad_age),
#                                        1e5*simul_hosp_sum_full/rsv_age_groups$stationary_popul[11],
#                                        1e5*simul_hosp_sum_full/sum(rsv_age_groups$stationary_popul[1:7]) ) )
#   #
#   gc()
#   # toc()
#   print(paste0(k_par," done"))
# }
# 
# # list to dataframe
# simul_hosp_rate_weekly_all_broad_agegroups_grad_relax <- bind_rows(simul_hosp_rate_weekly_all_broad_agegroups_grad_relax)
# simul_hosp_rate_weekly_under5_over65_grad_relax <- bind_rows(simul_hosp_rate_weekly_under5_over65_grad_relax)
