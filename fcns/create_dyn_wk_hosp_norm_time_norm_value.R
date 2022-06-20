# calculate statistics at time points relative to peak week of selected reference year,
# for each bin of the selected parameters

n_cnt=0; summ_dyn_all_parsets_broad_age_relat_time=list()
for (sel_parname in unique(binned_partable_regular_dyn$name)) {
  for (par_bin_val in sort(unique(binned_partable_regular_dyn$par_bin))) {
    sel_pars = binned_partable_regular_dyn$par_id[binned_partable_regular_dyn$name %in% sel_parname & 
                                                    binned_partable_regular_dyn$par_bin %in% par_bin_val]
    n_cnt=n_cnt+1; 
    med_parval=(binned_partable_regular_dyn %>% select(name,par_bin,median_parbin) %>% distinct() %>% 
                    filter(name %in% sel_parname & par_bin %in% par_bin_val))$median_parbin
    # calculate median and CIs
    summ_dyn_all_parsets_broad_age_relat_time[[n_cnt]] = dyn_hosp_weekly_norm_to_peak %>% 
      ungroup() %>% filter(par_id %in% sel_pars & epi_year<2024 & epi_year>=2017) %>%
      group_by(epi_year,peak_week_distance,agegroup) %>%
      summarise(median_date=median(date), # unique(paste0(epi_year,"-",isoweek(median(date)))),
                median=median(value_norm),mean=mean(value_norm),
                ci50_low=quantile(value_norm,c(0.25,0.75))[1],
                ci50_up=quantile(value_norm,c(0.25,0.75))[2],
                ci95_low=quantile(value_norm,c(0.025,0.975))[1],
                ci95_up=quantile(value_norm,c(0.025,0.975))[2]) %>%
      mutate(par_bin=par_bin_val,parname=sel_parname,par_median=med_parval,
             agegroup=c("<1y","1-2y","2-5y",NA,"65+y")[agegroup]) %>% 
      relocate(c(parname,par_bin,par_median),.before=epi_year)
    # print progress
    message(sel_parname," ",par_bin_val)
  }
}

summ_dyn_all_parsets_broad_age_relat_time = bind_rows(summ_dyn_all_parsets_broad_age_relat_time)
