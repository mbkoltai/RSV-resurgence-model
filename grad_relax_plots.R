# summary by year
epi_year_week_start <- isoweek(as.Date("2020-06-01"))
# normalise by pre-pandemic epi-years 2017 and 2018
summ_dyn_peak_cumul_meanage <- simul_hosp_rate_weekly_all_broad_agegroups_grad_relax %>%
  filter(par_id %in% partable_regular_dyn$par_id) %>%
  mutate(epi_year=ifelse(week(date)>=40,year(date),year(date)-1),
         epi_year_june=ifelse(week(date)>=epi_year_week_start,year(date),year(date)-1)) %>% ungroup() %>%
  filter(epi_year>2016) %>%
  pivot_longer(cols=c(epi_year,epi_year_june)) %>% rename(epi_year_wk_start=name,epi_year=value) %>%
  mutate(epi_year_wk_start=ifelse(grepl("june",epi_year_wk_start),epi_year_week_start,40)) %>%
  group_by(epi_year,epi_year_wk_start,par_id,agegroup) %>%
  summarise(peak_hosp=max(simul_hosp_sum_full), cumul_hosp=sum(simul_hosp_sum_full),
            peak_yday=mean(yday(date[simul_hosp_sum_full %in% max(simul_hosp_sum_full)]))) %>%
  pivot_longer(!c(epi_year,epi_year_wk_start,par_id,agegroup)) %>%
  group_by(par_id,agegroup,name,epi_year_wk_start) %>% # calculate pre-pandemic mean
  # because of some biennial patterns need to normalise by the right pre-pandemic year (2022/2018, 2021/2017)
  mutate(norm_value=case_when(grepl("hosp",name) ~ value/value[epi_year==ifelse(epi_year %% 2==0,2018,2017)],
                              !grepl("hosp",name) ~ value-value[epi_year==ifelse(epi_year %% 2==0,2018,2017)])) %>%
  relocate(par_id,.before=epi_year)
# add mean age
if (!any(grepl("age",unique(summ_dyn_peak_cumul_meanage$name)))) {
  summ_dyn_peak_cumul_meanage = 
    bind_rows(summ_dyn_peak_cumul_meanage,
              summ_dyn_peak_cumul_meanage %>% group_by(par_id,epi_year_wk_start,name,epi_year) %>%
                mutate(mean_age_under5y=ifelse(agegroup==1 & name %in% "cumul_hosp",sum(value[agegroup==1]*1/2+
                                              value[agegroup==1]*1.5+value[agegroup==1]*3.5)/sum(value),NA)) %>% 
                ungroup() %>% select(epi_year,epi_year_wk_start,par_id,mean_age_under5y) %>% 
                filter(!is.na(mean_age_under5y)) %>% distinct() %>% 
                group_by(par_id,epi_year_wk_start) %>%
                mutate(mean_age_shift=mean_age_under5y-mean_age_under5y[epi_year==ifelse(epi_year %% 2==0,2018,2017)]) %>%
                rename(value=mean_age_under5y,norm_value=mean_age_shift) %>% mutate(name="mean_age_under5y") )
}

# left-join with partable to have input parameters, calculate summary stats by quantiles of parameters
# summ_max_incid_seas_length_byvalue
summ_dyn_peak_cumul_meanage_byparvalue <- left_join(
  summ_dyn_peak_cumul_meanage %>% filter(epi_year>2019),
  left_join(partable_regular_dyn,pred_pca %>% select(par_id,PC1)) %>% 
    select(!const_delta) %>% rename(waning=omega,age_exp_PC1=PC1) %>% mutate(age_exp_ratio=age_dep/exp_dep) %>% 
    pivot_longer(!par_id) %>% rename(parname=name,parvalue=value), by="par_id") %>% 
  relocate(c(name,value,norm_value),.after=last_col()) %>% relocate(par_id,.before=epi_year) %>% rename(varname=name) %>%
  group_by(parname) %>% mutate(par_bin=ntile(parvalue,10)) %>% relocate(par_bin,.after=parvalue) %>% 
  # calculate summary stats (median, ci50, ci95)
  group_by(agegroup,epi_year,epi_year_wk_start,parname,par_bin,varname) %>%
  summarise(parvalue=median(parvalue),mean_parvalue=mean(parvalue),
            median=median(norm_value),mean=mean(norm_value),
            ci50_l=quantile(norm_value,probs=0.25),ci50_u=quantile(norm_value,probs=0.75),
            ci95_l=quantile(norm_value,probs=0.025),ci95_u=quantile(norm_value,probs=0.975))

##############################################################
# plot PEAK VALUE of cases/hospitalisations (+ season length), disaggregated by parameter values

# Figure 3B, SI 7-10
sel_pars <- unique(summ_dyn_peak_cumul_meanage_byparvalue$parname)
sel_varnames <- unique(summ_dyn_peak_cumul_meanage_byparvalue$varname)
start_year=2020
# this file also saves the plots!
source("fcns/create_summ_plots_from_dyn.R")

##############################################################
# calculate PRCC values wrt parameters

source("fcns/create_df_prcc_from_dynamics.R")
# creates the list: list_prcc_dyn

# plot PRCC on cumulative/peak hospitalisations, timing of peak (yday relative to pre-NPI), mean age of hosp
source("fcns/create_prcc_plot_from_dyn.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Evaluating DYNAMICS as a function of parameters (relative timing of peak and duration/shape of season)

# all_dynamics_accepted <- read_csv(here(foldername,"all_dynamics_accepted_2e3.csv")) %>%
#  mutate(date=as.Date(start_date_dyn_save)+t-min(t)) %>%
#   filter(date>=as.Date("2018-06-01") & date<as.Date("2024-10-01"))

epi_year_week_start=23
startweek_startyr=isoweek(min(simul_hosp_rate_weekly_all_broad_agegroups_grad_relax$date))-23
start_yr=isoyear(min(simul_hosp_rate_weekly_all_broad_agegroups_grad_relax$date))

dyn_hosp_weekly_norm_to_peak <- simul_hosp_rate_weekly_all_broad_agegroups_grad_relax %>% ungroup() %>% 
  filter(par_id %in% partable_regular_dyn$par_id & agegroup %in% 1:3) %>% rename(value=simul_hosp_sum_full) %>%
  # filter those params where the off-season is earlier
  filter(par_id %in% 
    (dist_hosp_2021_22 %>% filter(mean_sqrd_dist<=quantile(dist_hosp_2021_22$mean_sqrd_dist,probs=0.2)))$par_id) %>%
  mutate(epi_year=ifelse(week(date)>=epi_year_week_start,isoyear(date),isoyear(date)-1)) %>% 
  group_by(epi_year,par_id,agegroup) %>% 
  mutate(epi_week=ifelse(epi_year==start_yr, row_number()+startweek_startyr,row_number())) %>%    
  ungroup() %>% group_by(agegroup,par_id) %>%
  mutate(peak_value_2017=max(value[epi_year %in% 2017]),
         peak_value_2018=max(value[epi_year %in% 2018])) %>% # peak level of selected year
  group_by(agegroup,par_id) %>% 
  # peak week of year of reference  
  mutate(peak_week_2017=min(epi_week[value==peak_value_2017 & (epi_year %in% 2017)]),
         peak_week_2018=min(epi_week[value==peak_value_2018 & (epi_year %in% 2018)])) %>% 
  # distance from peak week of sel year    
  ungroup() %>% mutate(peak_week_distance=ifelse(epi_year %% 2==0,epi_week-peak_week_2018,epi_week-peak_week_2017)) %>% 
  mutate(value_norm=ifelse(epi_year %% 2==0,value/peak_value_2018,value/peak_value_2017))

# bin parameter values
binned_partable_regular_dyn <- partable_regular_dyn %>% select(!c(const_delta,seasforc_width_wks,peak_week)) %>%
  filter(par_id %in% unique(dyn_hosp_weekly_norm_to_peak$par_id)) %>%
  mutate(age_exp_ratio=age_dep/exp_dep) %>%
  pivot_longer(!par_id) %>% group_by(name) %>% mutate(par_bin=ntile(value,8)) %>% 
  group_by(name,par_bin) %>% mutate(median_parbin=median(value))

#######
# calculate statistics at time points relative to peak week of selected reference year,
# for each bin of the selected parameters
# this takes about 1 min (progress is printed)
source("fcns/create_dyn_wk_hosp_norm_time_norm_value.R")

# this creates the dataframe `summ_dyn_all_parsets_broad_age_relat_time`

##############################################################
# Plot dynamics faceted by the age vs immunity dependence of susceptibility (Figure 4)

# FIGURE 4 and SI FIG : plot faceted by years, color-coded by parameter values
source("fcns/create_norm_dyn_plots.R")

##############################################################
# plot parameter distribs for early off season outbreaks
cutoff_llh = quantile((all_likelihoods %>% filter(name %in% "complete likelihood" & accepted))$value,na.rm=T,probs=0.5)
subsample_par = (all_likelihoods %>% filter(name %in% "complete likelihood" & accepted & value<cutoff_llh))$par_id

# we calculate euclidean distance from observed hosp rate in 2021-22
dist_hosp_2021_22 = left_join(
  simul_hosp_rate_weekly_under5_over65_grad_relax %>%
    select(broad_age,par_id,date,simul_hosp_rate_100k) %>% ungroup() %>%
    filter(broad_age %in% "<5y") %>% 
    mutate(sel_par=ifelse(par_id %in% partable_regular_dyn$par_id,"accepted","rejected")) %>%
    # select best LLH params?
    filter((par_id %in% subsample_par) & date>=as.Date("2021-01-01") & date<=max(SARI_watch_all_hosp$date)),
  SARI_watch_all_hosp,by=c("date","year_week","broad_age")) %>%
  mutate(abs_dist=abs(rate_under5yrs-simul_hosp_rate_100k),sqrd_dist=abs_dist^2) %>%
  group_by(par_id,broad_age) %>% summarise(mean_abs_dist=mean(abs_dist),mean_sqrd_dist=mean(sqrd_dist))

early_off_season = (dist_hosp_2021_22 %>% filter(mean_sqrd_dist<=quantile(dist_hosp_2021_22$mean_sqrd_dist,probs=0.1)))$par_id



plot_partable_histogram_offseas = partable_regular_dyn %>% filter(par_id %in% subsample_par) %>%
  mutate(`early off season`=par_id %in% early_off_season) %>% select(!const_delta) %>%
  mutate(`waning (days)`=1/omega,`maximal forcing (% above baseline)`=1e2*seasforce_peak) %>% rename(`R0 baseline`=R0) %>%
  select(!c(omega,seasforce_peak)) %>% # `R0 peak`=R0*(1+seasforce_peak)
  rename(`age-dependence`=age_dep,`exposure-dependence`=exp_dep,
         `peak forcing (week)`=peak_week,`season width (weeks)`=seasforc_width_wks) %>%
  group_by(`early off season`) %>% pivot_longer(!c(par_id,`early off season`))

# plot
ggplot(plot_partable_histogram_offseas %>% filter(!name %in% "peak forcing (week)"),
       aes(x=`early off season`,y=value,color=`early off season`)) + 
  geom_jitter(width=0.4,alpha=1/4) + # geom_violin(fill=NA,show.legend=F) + 
  geom_boxplot(fill=NA,width=0.88,size=3/4,outlier.colour=NA,color="black") + # 
  facet_wrap(~name,scales="free_x",nrow = 4) + # scale_y_log10() + # scale_x_discrete(expand=expansion(0.03,0)) +
  scale_color_manual(values=c("grey","blue"),guide=guide_legend(override.aes=list(size=3))) + 
  xlab("") + ylab("parameter values") + theme_bw() + standard_theme + 
  theme(strip.text=element_text(size=13),axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),
        legend.position="top",legend.text=element_text(size=13),legend.title=element_text(size=13)) + coord_flip()
# save
ggsave("simul_output/2e3_accepted_linear_relaxing/param_distrib_offseas_timing_jitter.png",width=28,height=18,units="cm")
