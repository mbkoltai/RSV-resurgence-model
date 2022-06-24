list_prcc_output=list()
for (k_age in 1:length(unique(results_fullscan_hosp$agegroup_broad))) {
  
  n_agegr=unique(results_fullscan_hosp$agegroup_broad)[k_age]
  # 2020 and 2021 epi-years are treated together. bc if there was an outbreak in eg. June 2021, 
  # this is epi-year 2020, but should be counted as an off-season outbreak
  subset_prcc <- results_fullscan_hosp %>% 
    filter(agegroup_broad %in% n_agegr & (epi_year %in% c(2020,2021)) ) %>% 
    select(par_id,agegroup_broad,epi_year,exp_dep,age_dep,seasforc_width_wks,seasforce_peak,R0,omega,
           hosp_tot_norm,peak_hosp_norm,mean_age_cumul_hosp_shift) %>%
    group_by(par_id,exp_dep,age_dep,seasforc_width_wks,seasforce_peak,R0,omega) %>%
    summarise(hosp_tot_norm=sum(hosp_tot_norm),peak_hosp_norm=max(peak_hosp_norm),
              mean_age_shift_2020=mean_age_cumul_hosp_shift[epi_year==2020],
              mean_age_shift_2021=mean_age_cumul_hosp_shift[epi_year==2021]) %>%
    ungroup() %>% select(!c(par_id))
  
  n_var=4
  list_prcc_output[[k_age]] = bind_rows(
    epi.prcc(subset_prcc %>% 
               select(!c(peak_hosp_norm,mean_age_shift_2020,mean_age_shift_2021)),sided.test=2,conf.level=0.95) %>% 
      mutate(output="cumulative"), #[,-ncol(subset_prcc)]
    epi.prcc(subset_prcc %>% 
               select(!c(hosp_tot_norm,mean_age_shift_2020,mean_age_shift_2021)),sided.test=2,conf.level=0.95) %>% 
      mutate(output="peak"), #[,-ncol(subset_prcc)]
    epi.prcc(subset_prcc %>% 
               select(!c(hosp_tot_norm,peak_hosp_norm,mean_age_shift_2021)),sided.test=2,conf.level=0.95) %>% 
      mutate(output="mean_age_shift_2020"),
    epi.prcc(subset_prcc %>% 
               select(!c(hosp_tot_norm,peak_hosp_norm,mean_age_shift_2020)),sided.test=2,conf.level=0.95) %>% 
      mutate(output="mean_age_shift_2021")) %>% 
    mutate(name=rep(colnames(subset_prcc)[-((ncol(subset_prcc)-(n_var-1)):ncol(subset_prcc))],n_var) ) %>% 
    mutate(agegroup_broad=unique(results_fullscan_hosp$agegroup_broad)[k_age]) %>% 
    relocate(name,.before=est) %>% relocate(agegroup_broad,.after=name) %>%
    mutate(name=case_when(grepl("exp_dep",name) ~ "exposure-dependence", 
                          grepl("age_dep",name) ~ "age-dependence",
                          grepl("seasforce_peak",name) ~ "seasonal forcing (peak relative to baseline)",
                          grepl("seasforc_width_wks",name) ~ "seasonal forcing (duration)",
                          grepl("R0",name) ~ "R0 (baseline)",
                          grepl("omega",name) ~ "waning rate")) %>%
    mutate(filter_crit=
             ifelse((grepl("mean_age",output) & agegroup_broad %in% "<1y") | !grepl("mean_age",output),T,F)) %>%
    filter(filter_crit) %>% select(!filter_crit)
}
