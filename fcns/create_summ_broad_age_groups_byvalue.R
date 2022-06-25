
summ_broad_age_groups_byvalue <- results_fullscan_hosp %>% 
  filter(!agegroup_broad %in% "5+y" & 
           par_id %in% partable_regular_dyn$par_id & 
           epi_year>2019) %>%
  select(agegroup_broad,par_id,epi_year,exp_dep,age_dep,
         seasforc_width_wks,seasforce_peak,R0,omega,
         hosp_tot_norm,peak_hosp_norm,mean_age_cumul_hosp_shift) %>%
  mutate(age_exp_ratio=age_dep/exp_dep) %>%
  pivot_longer(!c(agegroup_broad,par_id,epi_year,hosp_tot_norm,
                  peak_hosp_norm,mean_age_cumul_hosp_shift)) %>%
  rename(parname=name,parvalue=value) %>% 
  relocate(c(parname,parvalue),.before=hosp_tot_norm) %>%
  mutate(epi_year_20_21_merged=
           ifelse(epi_year %in% c(2020,2021),"2020-21",epi_year)) %>% 
  # off-season outbreaks fall into epi-year 2020 -> 
  # take sum of 2020+2021 for cumulative and 20-21 epiyear maximum for peak
  group_by(agegroup_broad,par_id,epi_year_20_21_merged,parname) %>%
  summarise(parvalue=unique(parvalue),hosp_tot_norm=sum(hosp_tot_norm),
            peak_hosp_norm=max(peak_hosp_norm),
            mean_age_cumul_hosp_shift=max(mean_age_cumul_hosp_shift)) %>%
  group_by(parname) %>%
  mutate(par_bin=ntile(parvalue,10)) %>% relocate(par_bin,.after=parvalue) %>%
  pivot_longer(!c(agegroup_broad,par_id,epi_year_20_21_merged,
                  parname,parvalue,par_bin)) %>% rename(varname=name) %>%
  # calculate summary stats (median, ci50, ci95)
  group_by(agegroup_broad,epi_year_20_21_merged,parname,par_bin,varname) %>%
  summarise(parvalue=median(parvalue),mean_parvalue=mean(parvalue),
            median=median(value,na.rm=T),mean=mean(value,na.rm=T),
            ci50_l=quantile(value,probs=0.25,na.rm=T),
            ci50_u=quantile(value,probs=0.75,na.rm=T),
            ci95_l=quantile(value,probs=0.025,na.rm=T),
            ci95_u=quantile(value,probs=0.975,na.rm=T)) %>%
  group_by(parname,par_bin) %>% 
  mutate(parvalue=median(parvalue),mean_parvalue=mean(parvalue)) %>%
  filter(!(varname %in% "mean_age_cumul_hosp_shift" & !agegroup_broad %in% "<1y")) %>%
  mutate(agegroup_broad=ifelse(varname %in% "mean_age_cumul_hosp_shift",
                               "<5y",agegroup_broad)) %>%
  rename(epi_year=epi_year_20_21_merged)