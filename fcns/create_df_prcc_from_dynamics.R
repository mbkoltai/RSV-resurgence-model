# summ_dyn_peak_cumul_meanage

# this script calculates PRCC values between parameters and cumul/peak hosp post-NPI calculated from dynamics

df_dyn_prcc = left_join(partable_regular_dyn %>% select(!const_delta) %>% 
                          filter(par_id %in% unique(summ_dyn_peak_cumul_meanage$par_id)) %>% 
                          mutate(age_exp_dep_ratio=age_dep/exp_dep),
                        summ_dyn_peak_cumul_meanage) %>% select(!value) %>%
  mutate(agegroup=ifelse(is.na(agegroup),1,agegroup)) %>% filter(agegroup %in% 1:3 & epi_year>2019)

list_prcc_dyn=list()
# "peak_hosp"        "cumul_hosp"       "peak_yday"        "mean_age_under5y"
n_prcc_cnt=0
for (k_age in unique(df_dyn_prcc$agegroup)) {
  for (k_name in unique(df_dyn_prcc$name)) {
    for (k_epi_year in unique(df_dyn_prcc$epi_year)) {
      for (k_epiyear_start_wk in unique(df_dyn_prcc$epi_year_wk_start)){
        
        subset_prcc_dyn <- df_dyn_prcc %>% filter(agegroup %in% k_age & 
                                                    name %in% k_name &
                                                    epi_year %in% k_epi_year &
                                                    epi_year_wk_start %in% k_epiyear_start_wk) %>% 
          select(!c(epi_year,epi_year_wk_start,agegroup,par_id,name,peak_week))
        
        if (nrow(subset_prcc_dyn>0)) {
          # counter
          n_prcc_cnt=n_prcc_cnt+1
          list_prcc_dyn[[n_prcc_cnt]] = epi.prcc(subset_prcc_dyn,sided.test=2,conf.level=0.95) %>% 
            mutate(parname=colnames(subset_prcc_dyn)[-ncol(subset_prcc_dyn)] ) %>% 
            mutate(agegroup=k_age, output=k_name, epi_year=k_epi_year, epi_year_wk_start=k_epiyear_start_wk) %>% 
            relocate(c(parname,agegroup,epi_year,epi_year_wk_start,output),.before=est) %>% 
            mutate(parname=case_when(parname %in% "exp_dep" ~ "exposure-dependence", 
                                     grepl("age_dep",parname) ~ "age-dependence",
                                     grepl("seasforce_peak",parname) ~ "seasonal forcing (peak relative to baseline)",
                                     grepl("seasforc_width_wks",parname) ~ "seasonal forcing (duration)",
                                     grepl("R0",parname) ~ "R0 (baseline)",
                                     grepl("omega",parname) ~ "waning rate",
                                     parname %in% "age_exp_dep_ratio" ~ "age- to exposure-dependence ratio"),
                   `p_below_0.05`=p.value<0.05) 
          # progress
          message(paste0("calculating PRCC #",n_prcc_cnt,"/100")) }
        
      }
    }
  }
}

rm(df_dyn_prcc)