for (k in 1:ncol(median_parset)){ 
  if (k==1) { output_ranges_full_scan <- data.frame(); mean_age_shift_ranges<-data.frame() }
  
  # find parameter set closest to the median (closest by mean abs percentage error)
  error_threshold=5/100; n_par_scan=length(colnames(median_parset)[-k])
  parset_close_to_median <- 
    map2_df(abs(map2_df(results_fullscan_hosp[,colnames(median_parset)[-k]] %>% distinct(),
                        median_parset[,-k],`-`)),median_parset[,-k],`/`) %>% # this is calculating relative error (%)
    mutate(mean_abs_per_error=rowSums(across(where(is.numeric)))/n_par_scan,
           max_abs_per_error=pmax(!!!syms(colnames(median_parset)[-k])),
           par_id=unique(results_fullscan_hosp$par_id)[row_number()]) %>% 
    pivot_longer(!c(mean_abs_per_error,max_abs_per_error,par_id)) %>% group_by(name) %>% 
    mutate(error_ci95_low=quantile(value,probs=0.25)) %>% 
    group_by(par_id) %>% summarise(sum_threshold=sum(value<error_ci95_low),
                                   map_error=unique(mean_abs_per_error),
                                   max_abs_per_error=unique(max_abs_per_error)) %>%
    filter(max_abs_per_error<0.1) # sum_threshold==n_par_scan
  # outputs of median param set
  sel_all <- results_fullscan_hosp %>% filter(par_id %in% parset_close_to_median$par_id) %>%
    group_by(agegroup_broad,epi_year) %>% filter(epi_year==2021)
  # range of outputs for all param sets
  x_full_scan <- bind_rows(sel_all %>% summarise(norm_median=median(hosp_tot_norm,na.rm=T),
                                                 norm_min=min(hosp_tot_norm,na.rm=T),norm_max=max(hosp_tot_norm,na.rm=T),
                                                 norm_ci95_low=quantile(hosp_tot_norm,probs=0.05,na.rm=T),
                                                 norm_ci95_up=quantile(hosp_tot_norm,probs=0.95,na.rm=T)) %>% 
                             mutate(scan_param=colnames(median_parset)[k],range="full",vartype="cumulative"),
                           sel_all %>% summarise(norm_median=median(peak_hosp_norm,na.rm=T),
                                                 norm_min=min(peak_hosp_norm,na.rm=T),norm_max=max(peak_hosp_norm,na.rm=T),
                                                 norm_ci95_low=quantile(peak_hosp_norm,probs=0.05,na.rm=T),
                                                 norm_ci95_up=quantile(peak_hosp_norm,probs=0.95,na.rm=T)) %>%
                             mutate(scan_param=colnames(median_parset)[k],range="full",vartype="peak"))
  
  # accepted parameterisations
  sel_sel_parsets <- results_fullscan_hosp %>% 
    filter(par_id %in% parset_close_to_median$par_id) %>% 
    filter(par_id %in% partable_regular_dyn$par_id) %>%
    # right_join(median_parset[,-k],by=colnames(median_parset)[-k]) %>% 
    group_by(agegroup_broad,epi_year) %>% filter(epi_year==2021)
  # evaluate min, max, median
  x_parset_sel <- bind_rows(sel_sel_parsets %>% 
                              summarise(norm_median=median(hosp_tot_norm,na.rm=T),
                                        norm_min=min(hosp_tot_norm,na.rm=T),norm_max=max(hosp_tot_norm,na.rm=T),
                                        norm_ci95_low=quantile(hosp_tot_norm,probs=0.05,na.rm=T),
                                        norm_ci95_up=quantile(hosp_tot_norm,probs=0.95,na.rm=T)) %>%
                              mutate(scan_param=colnames(median_parset)[k],range="sel",vartype="cumulative"),
                            sel_sel_parsets %>% 
                              summarise(norm_median=median(peak_hosp_norm,na.rm=T),
                                        norm_min=min(peak_hosp_norm,na.rm=T),norm_max=max(peak_hosp_norm,na.rm=T),
                                        norm_ci95_low=quantile(peak_hosp_norm,probs=0.05,na.rm=T),
                                        norm_ci95_up=quantile(peak_hosp_norm,probs=0.95,na.rm=T)) %>%
                              mutate(scan_param=colnames(median_parset)[k],range="sel",vartype="peak"))
  # collect all results (cumul, peak)
  output_ranges_full_scan <- bind_rows(output_ranges_full_scan,x_full_scan,x_parset_sel)
  # mean age shift 
  x_mean_age <- bind_rows(sel_all %>% 
                            summarise(norm_median=median(mean_age_cumul_hosp_shift,na.rm=T),
                                      norm_min=min(mean_age_cumul_hosp_shift,na.rm=T),norm_max=max(mean_age_cumul_hosp_shift,na.rm=T),
                                      norm_ci95_low=quantile(mean_age_cumul_hosp_shift,probs=0.05,na.rm=T),
                                      norm_ci95_up=quantile(mean_age_cumul_hosp_shift,probs=0.95,na.rm=T)) %>%
                            mutate(scan_param=colnames(median_parset)[k],range="full"),
                          sel_sel_parsets %>% summarise(norm_median=median(mean_age_cumul_hosp_shift,na.rm=T),
                                                        norm_min=min(mean_age_cumul_hosp_shift,na.rm=T),norm_max=max(mean_age_cumul_hosp_shift,na.rm=T),
                                                        norm_ci95_low=quantile(mean_age_cumul_hosp_shift,probs=0.05,na.rm=T),
                                                        norm_ci95_up=quantile(mean_age_cumul_hosp_shift,probs=0.95,na.rm=T)) %>%
                            mutate(scan_param=colnames(median_parset)[k],range="sel")) %>% 
    group_by(scan_param,range) %>% 
    select(!agegroup_broad) %>% 
    summarise_all(unique)
  # 
  mean_age_shift_ranges <- bind_rows(x_mean_age,mean_age_shift_ranges)
  
  if (k==ncol(median_parset)) {
    output_ranges_full_scan <- output_ranges_full_scan %>% 
      relocate(c(scan_param,range,vartype),.after=epi_year) }
}
### END OF LOOP (building dataframe for Fig 2)