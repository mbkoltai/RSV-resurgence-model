startweek_startyr=isoweek(min(simul_hosp_rate_weekly_all_broad_agegroups$date))-23
start_yr=isoyear(min(simul_hosp_rate_weekly_all_broad_agegroups$date))

dyn_hosp_weekly_norm_to_peak <- simul_hosp_rate_weekly_all_broad_agegroups %>% ungroup() %>% 
  filter(par_id %in% partable_regular_dyn$par_id & agegroup %in% 1:3) %>% 
  rename(value=simul_hosp_sum_full) %>%
  mutate(epi_year=ifelse(week(date)>=epi_year_week_start,
                         year(date),year(date)-1)) %>% 
  group_by(epi_year,par_id,agegroup) %>% 
  mutate(epi_week=ifelse(epi_year==start_yr,
                         row_number()+startweek_startyr,row_number())) %>%    
  ungroup() %>% group_by(agegroup,par_id) %>%
  mutate(peak_value_2017=max(value[epi_year %in% 2017]),
         peak_value_2018=max(value[epi_year %in% 2018])) %>% # peak level of selected year
  group_by(agegroup,par_id) %>% 
  # peak week of year of reference  
  mutate(peak_week_2017=min(epi_week[value==peak_value_2017 & (epi_year %in% 2017)]),
         peak_week_2018=min(epi_week[value==peak_value_2018 & (epi_year %in% 2018)])) %>% 
  # distance from peak week of sel year    
  ungroup() %>% 
  mutate(peak_week_distance=ifelse(epi_year %% 2==0,
                                   epi_week-peak_week_2018,epi_week-peak_week_2017)) %>% 
  mutate(value_norm=ifelse(epi_year %% 2==0,
                           value/peak_value_2018,value/peak_value_2017))

# bin parameter values 
binned_partable_regular_dyn <- partable_regular_dyn %>% 
  select(!c(const_delta,seasforc_width_wks,peak_week)) %>%
  filter(par_id %in% unique(dyn_hosp_weekly_norm_to_peak$par_id)) %>%
  mutate(age_exp_ratio=age_dep/exp_dep) %>%
  pivot_longer(!par_id) %>% group_by(name) %>% 
  mutate(par_bin=ntile(value,8)) %>% 
  group_by(name,par_bin) %>% mutate(median_parbin=median(value))