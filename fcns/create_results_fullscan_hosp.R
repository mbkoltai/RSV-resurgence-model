mean_ages_broad_age_group <- data.frame(agegroup_broad=1:3,
                                        mean_age=c(0.5,1.5,3.5))

# bin values
results_fullscan_hosp <- bind_rows(results_summ_all_reg_dyn,
                                   results_summ_rejected) %>%
  mutate(hosp_tot=inf_tot*hosp_probabilities$prob_hosp_per_infection[agegroup],
         peak_hosp=peak_inf*hosp_probabilities$prob_hosp_per_infection[agegroup]) %>% 
  mutate(agegroup_broad=c("<1y","1-2y","2-5y","5+y")[findInterval(agegroup,c(2,4,7)+1)+1]) %>%
  group_by(par_id,agegroup_broad,epi_year,exp_dep,age_dep,
           seasforc_width_wks,seasforce_peak,R0,omega) %>% 
  summarise(hosp_tot=sum(hosp_tot),
            peak_hosp=sum(peak_hosp),
            inf_tot=sum(inf_tot)) %>% 
  group_by(agegroup_broad,par_id) %>%
  mutate(hosp_tot_norm=hosp_tot/hosp_tot[epi_year==2019],
         peak_hosp_norm=peak_hosp/peak_hosp[epi_year==2019]) %>%
  mutate(kappa_grouped=findInterval(kappa,seq(-1,1,by=1/5)) ) %>% 
  group_by(kappa_grouped) %>% 
  mutate(kappa_grouped=round(mean(kappa),1)) %>% 
  relocate(c(kappa,kappa_grouped),.after=age_dep)

# mean age should be calculated from narrower age bands
if (!any(grepl("mean_age",colnames(results_fullscan_hosp)))){
  xx = left_join(bind_rows(results_summ_all_reg_dyn,
                           results_summ_rejected) %>%
                   mutate(hosp_tot=inf_tot*hosp_probabilities$prob_hosp_per_infection[agegroup],
                          peak_hosp=peak_inf*hosp_probabilities$prob_hosp_per_infection[agegroup]) %>%
                   filter(agegroup<=7) %>% # 
                   group_by(par_id,epi_year) %>%
                   summarise(mean_age_cumul_hosp_under5y=
                               sum(rsv_age_groups$mean_age_weighted[agegroup]*hosp_tot/sum(hosp_tot))*12,
                             mean_age_peak_hosp_under5y=
                               sum(rsv_age_groups$mean_age_weighted[agegroup]*peak_hosp/sum(peak_hosp))*12),
                 partable %>% select(!const_delta)) %>%
    relocate(c(mean_age_cumul_hosp_under5y,mean_age_peak_hosp_under5y),
             .after=last_col()) %>%
    group_by(par_id) %>%
    mutate(mean_age_cumul_hosp_shift=mean_age_cumul_hosp_under5y - 
             mean_age_cumul_hosp_under5y[epi_year==2017],
           mean_age_cumul_hosp_shift_2018=mean_age_cumul_hosp_under5y - 
             mean_age_cumul_hosp_under5y[epi_year==2018],
           mean_age_cumul_hosp_shift=ifelse(epi_year %% 2 == 0,
                                            mean_age_cumul_hosp_shift_2018,
                                            mean_age_cumul_hosp_shift)) %>%
    select(!mean_age_cumul_hosp_shift_2018)
  # concatenate
  results_fullscan_hosp = left_join(results_fullscan_hosp,xx)
}

# find median parameter set
median_parset <- results_fullscan_hosp %>% ungroup() %>% 
  select(c(kappa_grouped,seasforc_width_wks,seasforce_peak,R0,omega)) %>% 
  group_by(kappa_grouped,seasforc_width_wks,seasforce_peak,R0,omega) %>% 
  summarise_all(unique) %>% ungroup() %>%
  summarise_all(quantile,p=0.5,type=1)