create_summ_dyn_peak=F

if (create_summ_dyn_peak) {

unzip("repo_data/simul_hosp_rate_weekly_all_broad_agegroups.zip",exdir="repo_data/")
simul_hosp_rate_weekly_all_broad_agegroups <- read_csv("repo_data/simul_hosp_rate_weekly_all_broad_agegroups.csv")

# summary by year
epi_year_week_start <- isoweek(as.Date("2020-06-01"))
# normalise by pre-pandemic epi-years 2017 and 2018
summ_dyn_peak_cumul_meanage <- simul_hosp_rate_weekly_all_broad_agegroups %>%
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
  # because of some biennial patterns need to normalise 
  # by the right pre-pandemic year (2022/2018, 2021/2017)
  mutate(norm_value=case_when(
    grepl("hosp",name) ~ value/value[epi_year==ifelse(epi_year %% 2==0,2018,2017)],
    !grepl("hosp",name) ~ value-value[epi_year==ifelse(epi_year %% 2==0,2018,2017)])) %>%
  relocate(par_id,.before=epi_year)
# add mean age
if (!any(grepl("age",unique(summ_dyn_peak_cumul_meanage$name)))) {
  summ_dyn_peak_cumul_meanage = 
    bind_rows(summ_dyn_peak_cumul_meanage,
              summ_dyn_peak_cumul_meanage %>% 
                group_by(par_id,epi_year_wk_start,name,epi_year) %>%
                mutate(mean_age_under5y=ifelse(agegroup==1 & 
                          name %in% "cumul_hosp",sum(value[agegroup==1]*1/2+
                          value[agegroup==1]*1.5+value[agegroup==1]*3.5)/sum(value),NA)) %>% 
                ungroup() %>% select(epi_year,epi_year_wk_start,par_id,mean_age_under5y) %>% 
                filter(!is.na(mean_age_under5y)) %>% distinct() %>% 
                group_by(par_id,epi_year_wk_start) %>%
                mutate(mean_age_shift=
                         mean_age_under5y-
                         mean_age_under5y[epi_year==ifelse(epi_year %% 2==0,2018,2017)]) %>%
                rename(value=mean_age_under5y,
                       norm_value=mean_age_shift) %>% mutate(name="mean_age_under5y") )
}

# left-join with partable to have input parameters, calculate summary stats by quantiles of parameters
# summ_max_incid_seas_length_byvalue
summ_dyn_peak_cumul_meanage_byparvalue <- left_join(
  summ_dyn_peak_cumul_meanage %>% filter(epi_year>2019),
  partable_regular_dyn %>% 
    select(!const_delta) %>% 
    rename(waning=omega) %>% 
    mutate(age_exp_ratio=age_dep/exp_dep) %>% 
    pivot_longer(!par_id) %>% rename(parname=name,parvalue=value), by="par_id") %>% 
  relocate(c(name,value,norm_value),.after=last_col()) %>% 
  relocate(par_id,.before=epi_year) %>% rename(varname=name) %>%
  group_by(parname) %>% mutate(par_bin=ntile(parvalue,10)) %>% 
  relocate(par_bin,.after=parvalue) %>% 
  # calculate summary stats (median, ci50, ci95)
  group_by(agegroup,epi_year,epi_year_wk_start,parname,par_bin,varname) %>%
  summarise(parvalue=median(parvalue),mean_parvalue=mean(parvalue),
            median=median(norm_value),mean=mean(norm_value),
            ci50_l=quantile(norm_value,probs=0.25),
            ci50_u=quantile(norm_value,probs=0.75),
            ci95_l=quantile(norm_value,probs=0.025),
            ci95_u=quantile(norm_value,probs=0.975))
} else {
  summ_dyn_peak_cumul_meanage <- read_csv("repo_data/summ_dyn_peak_cumul_meanage.csv")
  summ_dyn_peak_cumul_meanage_byparvalue <- read_csv("repo_data/summ_dyn_peak_cumul_meanage_byparvalue.csv")
}

# summ_dyn_peak_cumul_meanage
# summ_dyn_peak_cumul_meanage_byparvalue