# Run individual simulations
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
rm(list=ls())
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
# here::i_am()
# currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# load libraries
x1 <- c("tidyverse","deSolve","lubridate","gtools","ISOweek","lhs")
# "rstudioapi","devtools","wpp2019","Rcpp","tsibble","pracma","qs","ungeviz","zoo","RcppRoll"
# if (!any(grepl("ungeviz",row.names(installed.packages())))) { devtools::install_github("wilkelab/ungeviz") }
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
# Load all packages
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
# select <- dplyr::select; # row_number <- dplyr::row_number; summarise <- dplyr::summarise
source('fcns/essential_fcns.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
source("load_params.R")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SARI-Watch data with NUMBER of hospitalisations
# library(ISOweek)
SARIwatch_RSVhosp_under5_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_under5_2018_2020_weekly_counts.csv",col_types="ffddd") %>% 
  mutate(wk_n=gsub("\\.","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),
                                             date=ISOweek2date(paste0(gsub("\\.","-W",wk_n),"-1")))
# compare SARI-Watch with literature estimates 
# annual RATE of hospitalisations: 500/1e5 (18-19), 494/1e5
SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
  summarise(annual_cumul_rate=sum(rate_under5yrs))

# annual RATE from Reeves 2017: <1y: 35.1/1000, 1-4y: 5.31/1000
# = 1083.8/1e5
(35.1*sum(rsv_age_groups$value[1:2]) + 5.31*sum(rsv_age_groups$value[3:7]))*100/sum(rsv_age_groups$value[1:7])

# annual RATE from Taylor 2016: <0.5y: 4184/1e5, 6-23mts: 1272/1e5, 2-4y: 114/1e5
# = 822.7/1e5
(4184*sum(rsv_age_groups$value[1])+1272*sum(rsv_age_groups$value[2:4]) + 
    114*sum(rsv_age_groups$value[5:7]))/sum(rsv_age_groups$value[1:7])

# underreporting should be ratio of (SARI-Watch rate) to (lit estim rate)
under_report_factor_under5 = mean((SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% group_by(year) %>% 
                            summarise(annual_cumul_rate=sum(rate_under5yrs)))$annual_cumul_rate) / 
  mean(c( (35.1*sum(rsv_age_groups$value[1:2]) + 
            5.31*sum(rsv_age_groups$value[3:7]))*100/sum(rsv_age_groups$value[1:7]),
         (4184*sum(rsv_age_groups$value[1])+1272*sum(rsv_age_groups$value[2:4])+
            114*sum(rsv_age_groups$value[5:7]))/sum(rsv_age_groups$value[1:7])))

### ### ### ### ### ### ### ### ### ### ### ### 
# for 65+y
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts <- 
  read_csv("data/SARIwatch_RSVhosp_over65y_2018_2020_weekly_counts.csv") %>%
    mutate(wk_n=gsub("-","-W",wk_n),wk_n=factor(wk_n,levels=unique(wk_n)),date=ISOweek2date(paste0(wk_n,"-1")))

# plot
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% pivot_longer(!c(year,wk_n,pop_AGE65PLUS,date)) %>% 
ggplot(aes(x=wk_n,y=value,group=1)) + geom_line() + geom_point() + facet_grid(name~year,scales="free") + 
  xlab("") + ylab("count or rate per 100K") + theme_bw() + standard_theme

# Estimates from literature (https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-015-1218-z/tables/3)
# 65-74: 86/100e3 (62-101); 75+: 234/100e3 (180-291)
# from SARI-Watch
SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% group_by(year) %>% summarise(annual_rate=sum(rate_65yplus))
# from lit estimate
over65_hosp_rate_100k_lit_estim = (sum(ons_2020_midyear_estimates_uk$value[ons_2020_midyear_estimates_uk$age %in% 65:74])*86 + 
  sum(ons_2020_midyear_estimates_uk$value[76:nrow(ons_2020_midyear_estimates_uk)])*234)/
  sum(ons_2020_midyear_estimates_uk$value[66:nrow(ons_2020_midyear_estimates_uk)])
# under-reporting rate
under_report_factor_over65y <- median((SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% 
      group_by(year) %>% summarise(annual_rate=sum(rate_65yplus)))$annual_rate)/over65_hosp_rate_100k_lit_estim

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set variable/control parameters
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
exp_dep <- 2/3; age_dep <- 1/10  # partable$exp_dep[k_par] # exp=1.75,agedep=0.3125 --> const delta: 2.55
const_delta <- 0.3; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # partable$const_delta[k_par]
# matrix of susceptibility parameters
delta_susc <- round(sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))}),5)

# calculate R0 *at the baseline* (no forcing)
R0_calc_SIRS(C_m,delta_susc,rho,n_inf)

### ### ### ### ### ### ### ### ### ### ### ### 
# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) 
# npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-05-17"))
npi_strengh=0 # 1 = 100% reduction
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,
                                   years_pre_post_npi=c(15,0),season_width_wks=8,init_mt_day="01-01", # "06-01"
                                   peak_week=44,forcing_above_baseline=1/2,npireduc_strength=npi_strengh)
# plot seasonal forcing
# fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))

# check forcing
# ggplot(data.frame(date=unique(df_cases_infs$date),forcing=forcing_vector_npi) %>%
#          mutate(nday=row_number(),year=year(date))  ) + # %>% filter(nday>200&nday<2.25e3)
#   geom_line(aes(x=yday(date),y=forcing,color=factor(year))) + scale_x_continuous(breaks=(0:52)*7,expand=expansion(0.01,0)) + 
#   geom_vline(xintercept=yday("2020-11-01"),color="red") + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### 
# all parameters collected as a list (sent to fcn running ODE solver)
deaths_vector <- matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))
params <- list(list(birth_rates,deaths_vector),K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc,mat_imm_inds)

# INITIAL CONDITION
# init_set can be "from scratch" (all susceptible) OR last state of a previous simulation
# initvals_sirs_model = round(matrix(ode_sol[1700,(1:(ncol(ode_sol)-1))+1])); initvals_sirs_model[100:132]=0
# initvals_sirs_model <- fcn_set_initconds(rsv_agegroup_sizes=rsv_age_groups$stationary_popul,
#                           init_set=c("previous","fromscratch")[ifelse(exists("ode_solution"),1,2)],init_cond_src=c("output","file")[1],
#                           input_from_prev_simul=ode_solution,init_seed=10,seed_vars="all",filename="")
# take last time point of existing simulation
# if (exists("ode_sol")){
#   initvals_sirs_model=round(matrix(ode_sol[nrow(ode_sol)-1,(1:(ncol(ode_sol)-1))+1])); initvals_sirs_model[100:132]=0
# } else { # or read in a stationary solution from file
#   initvals_sirs_model<-readRDS("extras/stationary_sol.RDS")
# }
# read in a demographically stable solution
initvals_sirs_model <- readRDS("extras/stationary_sol.RDS")

# SOLVE ODE SYSTEM
# can use 'lsoda' (faster) or 'lsodes'
tm <- proc.time(); ode_sol <- lsoda(initvals_sirs_model,timesteps,
                                    func=sirs_seasonal_forc_mat_immun,parms=params); round(proc.time()-tm,2)
# Extract infection *incidence*
df_cases_infs <- fcn_process_odesol_incid(ode_sol,n_age,n_inf,n_compartment,simul_start_end) 
# %>% mutate(value=round(value,2))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# this step only needed if want to investigate other variables or want to use it as input for next simulation (initial conds)
# reshape data: get all variables in `ode_solution_tidy` 
# g(final_pop,ode_solution,ode_solution_tidy) %=% fun_process_simul_output(ode_sol,varname_list,
#                incidvar="newinf",incid_only=F,init_date=simul_start_end[1],n_age,n_inf,rsv_age_groups,neg_thresh=-0.01)
# # Has population size converged? (PERCENTAGE difference!)
# c(round(100*abs(final_pop[,3]-final_pop[,2])/final_pop[,2],2))
round(100*abs(sapply(1:n_age, function(x) ( sum(ode_sol[1,(x-1)*9+(1:9)+1])-
                                             sum(ode_sol[nrow(ode_sol)-1,(x-1)*9+(1:9)+1]) ) ))/
                                          sapply(1:n_age, function(x) sum(ode_sol[1,(x-1)*9+(1:9)+1]) ),3)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# 1/2/3rd infections together
df_cases_infs %>% filter(agegroup<=7) %>% group_by(date) %>% 
  summarise(value=sum(value)) %>% filter(date>as.Date("2014-09-01")) %>%
# mutate(epi_date=ifelse(yday(date)>220,yday(date)-220,yday(date)+145)) %>% filter(epi_date<220 & epi_date>40 & year(date)<=2017) %>%
ggplot(aes(x=date,y=value)) + geom_line() +
  # ggplot(aes(x=epi_date,y=value,color=factor(year(date)))) + geom_line() + # facet_wrap() +
  scale_x_date(date_breaks="3 month",expand=expansion(0.01/3,0)) + # scale_x_continuous(breaks=(0:20)*20,expand=expansion(0.01,0)) + 
  scale_y_continuous(breaks=(0:30)*2.5e3,expand=expansion(0,0)) + 
  geom_vline(data=df_cases_infs %>% filter(yday(date)==yday("2017-12-05")) %>% 
               ungroup() %>% select(date) %>% unique(),
              aes(xintercept=date),size=1/3,linetype="dashed") + 
  xlab("") + ylab("new infections/day") +
  ggtitle("daily incident infections in <5y") + theme_bw() + standard_theme + theme(legend.position="bottom") + labs(color="year")

# show years overlaid with years as diff colours
df_cases_infs %>% filter(agegroup<=7) %>% group_by(date) %>% summarise(value=sum(value)) %>% 
  filter(date>as.Date("2014-09-01")) %>% mutate(year=year(date),yday=yday(date)) %>%
ggplot(aes(x=yday,y=value,color=factor(year),group=year)) + geom_line() + 
  scale_x_continuous(breaks=(0:52)*7,expand=expansion(0.01,0)) +
  theme_bw() + standard_theme # scale_color_gradient(low="blue",high="red") + 

# what is the attack rate?
attack_rates <- left_join(df_cases_infs %>% filter(date>as.Date("2019-10-01") & date<as.Date("2020-04-01")) %>%
            group_by(t,date,agegroup) %>% summarise(value=sum(value)) %>% 
              group_by(agegroup) %>% summarise(cumul_inf=sum(value)),
          rsv_age_groups %>% 
            select(agegroup_name,stationary_popul,value) %>% rename(real_pop=value) %>%
            mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))),by="agegroup" ) %>%
            mutate(attack_rate_stat_pol=cumul_inf/stationary_popul,attack_rate_real_pol=cumul_inf/real_pop)
# plot
ggplot(attack_rates) + 
  geom_segment(aes(x=agegroup-1/2,xend=agegroup+1/2,y=attack_rate_stat_pol*1e2,yend=attack_rate_stat_pol*1e2),size=1.4) +
  scale_x_continuous(breaks=1:11,expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:15)*10) + 
  geom_text(aes(x=agegroup,y=5,label=paste0(agegroup_name,"y"))) +
  geom_vline(xintercept=(0:11)+1/2,size=1/3,linetype="dashed") + xlab("age groups") + ylab("% attack rate") +
  theme_bw() + standard_theme

# for entire population (in millions)
attack_rates %>% summarise(cumul_inf=round(sum(cumul_inf)/1e6,3),stationary_popul=round(sum(stationary_popul)/1e6,3),
                           real_pop=round(sum(real_pop)/1e6,3),AR_percent=100*cumul_inf/stationary_popul)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# load hospitalisation probability data
source(here("fcns", "calc_hosp_rates.R"))

# broad agegroups from original agegroups:
# findInterval(1:11,c(2,4,7,10)+1)+1

# calculate weekly under-5y hospitalisations
simul_hosp_rate <- left_join(
  all_dynamics_accepted %>% 
    filter(agegroup<=3 & date>=min(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts$date)) %>% 
    mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],unique(rsv_age_groups$agegroup_name))) %>% 
    group_by(t,date,agegroup_name) %>% summarise(value=sum(value)),
    # data on hosp
    hosp_probabilities %>% 
    select(agegroup_name,prob_hosp_per_infection_adj,prob_hosp_per_infection),by=c("agegroup_name")) %>%
    mutate(hosp=value*prob_hosp_per_infection*under_report_factor_under5,
         year_week=paste0(isoyear(date),"-",ifelse(nchar(isoweek(date))>1,isoweek(date),paste0("0",isoweek(date))) ),
         year_week=factor(year_week,levels=unique(year_week))) %>% 
  group_by(year_week,agegroup_name) %>% summarise(cases=sum(value),hosp=sum(hosp),date=min(date)) %>% 
  mutate(hosp_rate_per_100k=hosp/(sum(rsv_age_groups$stationary_popul[1:7])/1e5)) %>% 
  relocate(date,.before=year_week)
# rate per 100K is for ALL <5yr children, not by subgroups within <5yr

# PLOT simulated hospitalisations
hosp_plot_flag=F
if (hosp_plot_flag){
ggplot(simul_hosp_rate %>% filter(date<as.Date("2020-04-01") & date>as.Date("2018-10-01") ) ) +
  geom_bar(aes(x=year_week,y=hosp_rate_per_100k,fill=fct_rev(agegroup_name)),stat="identity",color="black",size=1/6) +
  scale_x_discrete(breaks=every_nth(n=2)) + scale_y_continuous(expand=expansion(0.01,0)) + 
  xlab("") + ylab("weekly hospit per 100k") + theme_bw() + standard_theme + theme(legend.title=element_blank())
}

# Merge simulated hospitalisations with SARI-Watch counts data, calculate likelihood
# calculate negative log-likelihood: -sum(dpois(x=x_data,lambda=x_model,log=T))
model_data_merged <- left_join(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts,
                        simul_hosp_rate %>% group_by(date) %>% 
                          summarise(simul_hosp_sum_full=sum(hosp),
                              simul_hosp_rate_per_100k=1e5*simul_hosp_sum_full/sum(rsv_age_groups$stationary_popul[1:7])) ) %>%
  mutate(simul_hosp_scaled=simul_hosp_sum_full*pop_AGEUNDER5/sum(rsv_age_groups$stationary_popul[1:7]),
         log_lklh_poiss=dpois(x=casesunder5total,lambda=simul_hosp_scaled,log=T)) %>% relocate(date,.after=wk_n)

#
# PLOT COUNTS (data) VS simulated hospits (scaled to catchment areas): under-5Y
model_data_merged %>% mutate(week=week(date)) %>% # filter(week>=40 | week<=12) %>%
ggplot(aes(x=wk_n)) + 
  geom_point(aes(y=casesunder5total,color="data"),size=2) + geom_line(aes(y=casesunder5total),group=1,color="red") +
  geom_point(aes(y=simul_hosp_scaled,color="simulations")) + geom_line(aes(y=simul_hosp_scaled),group=1,color="black") +
  geom_segment(aes(y=casesunder5total,yend=simul_hosp_scaled,xend=wk_n),size=1/3,linetype="dashed") +
  facet_wrap(~year,scales="free_x") + scale_color_manual(values=c("data"="red","simulations"="black")) + labs(color="") +
  geom_text(label=paste0("Neg. log. llh.=",round(-sum(model_data_merged$log_lklh_poiss,na.rm=T))),x=22,y=150,check_overlap=TRUE) +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*20) + # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  xlab("") + ylab("reported hospitalisations <5y") + theme_bw() + standard_theme + 
  theme(legend.text=element_text(size=14),legend.position="top")

# PLOT RATES/100K
model_data_merged %>% mutate(week=week(date)) %>% filter(week>=40 | week<=12) %>%
ggplot(aes(x=wk_n)) + 
  geom_point(aes(y=simul_hosp_rate_per_100k,color="simulations")) + geom_line(aes(y=simul_hosp_rate_per_100k),group=1,color="black") +
  # data
  geom_point(aes(y=rate_under5yrs,color="data"),size=2) + geom_line(aes(y=rate_under5yrs),group=1,color="red") +
  geom_segment(aes(y=rate_under5yrs,yend=simul_hosp_rate_per_100k,xend=wk_n),size=1/3,linetype="dashed") +
  facet_wrap(~year,scales="free_x") + scale_color_manual(values=c("data"="red","simulations"="black")) + labs(color="") +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*5) + # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  xlab("") + ylab("rate hospitalisations/100K <5y") + theme_bw() + standard_theme + theme(legend.text=element_text(size=14))

### ### ### ### ### ### ### 
# hosp counts for 65+y
simul_hosp_over65 = left_join(
  SARIwatch_RSVhosp_over65_2018_2020_weekly_counts,
  left_join(df_cases_infs %>% filter(agegroup==11 & date>as.Date("2018-01-01")),
          hosp_probabilities %>% mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))) %>%
            select(c(agegroup,agegroup_name,prob_hosp_per_infection)) ) %>% 
      mutate(hosp=value*prob_hosp_per_infection*under_report_factor_over65y,
         wk_n=paste0(isoyear(date),"-W",ifelse(nchar(isoweek(date))>1,isoweek(date),paste0("0",isoweek(date))) ),
         wk_n=factor(wk_n,levels=unique(wk_n))) %>% group_by(wk_n,agegroup_name) %>% 
      summarise(cases=sum(value),simul_hosp_sum_full=sum(hosp),date=min(date))  ) %>%
  mutate(simul_hosp_scaled=simul_hosp_sum_full*pop_AGE65PLUS/rsv_age_groups$stationary_popul[11],
         simul_hosp_rate_per100k=1e5*simul_hosp_sum_full/rsv_age_groups$stationary_popul[11],
         log_lklh_poiss=dpois(x=cases65plustotal,lambda=simul_hosp_scaled,log=T)) %>% relocate(date,.after=wk_n)

# plot 65+y hospitalisation COUNTS 
simul_hosp_over65 %>% 
  ggplot(aes(x=wk_n)) + 
  geom_point(aes(y=cases65plustotal,color="data"),size=2) + geom_line(aes(y=cases65plustotal),group=1,color="red") +
  geom_point(aes(y=simul_hosp_scaled,color="simulations")) + geom_line(aes(y=simul_hosp_scaled),group=1,color="black") +
  geom_segment(aes(y=cases65plustotal,yend=simul_hosp_scaled,xend=wk_n),size=1/3,linetype="dashed") +
  facet_wrap(~year,scales="free_x") + scale_color_manual(values=c("data"="red","simulations"="black")) + labs(color="") +
  geom_text(label=paste0("Neg. log. llh.=",round(-sum(simul_hosp_over65$log_lklh_poiss,na.rm=T))),x=22,y=50,check_overlap=TRUE) +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*20) + # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  xlab("") + ylab("reported hospitalisations 65+y") + theme_bw() + standard_theme + 
  theme(legend.text=element_text(size=14),legend.position="top")

# plot 65+y hospitalisation RATES (per 100K)
simul_hosp_over65 %>% ggplot(aes(x=wk_n)) +
  geom_point(aes(y=rate_65yplus,color="data"),size=2) + geom_line(aes(y=rate_65yplus),group=1,color="red") +
  geom_point(aes(y=simul_hosp_rate_per100k,color="simulations")) + geom_line(aes(y=simul_hosp_rate_per100k),group=1,color="black") +
  geom_segment(aes(y=rate_65yplus,yend=simul_hosp_rate_per100k,xend=wk_n),size=1/3,linetype="dashed") +
  facet_wrap(~year,scales="free_x") + scale_color_manual(values=c("data"="red","simulations"="black")) + labs(color="") +
  geom_text(label=paste0("Neg. log. llh.=",round(-sum(simul_hosp_over65$log_lklh_poiss,na.rm=T))),x=22,y=50,check_overlap=TRUE) +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*20) + # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  xlab("") + ylab("65+y hospitalisations/100e3 population") + theme_bw() + standard_theme + 
  theme(legend.text=element_text(size=14),legend.position="top")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# compare to DATAMART RSV detections
datamart_0_5y_2012_2022 <- read_csv("data/datamart_0_5y_2012_2022.csv")
# PLOT
simul_hosp_rate %>% filter(date<as.Date("2020-04-01") & date>as.Date("2012-07-01")) %>% group_by(date) %>% 
  summarise(simul_hosp_sum_full=sum(hosp)) %>%
  ggplot(aes(x=date,y=simul_hosp_sum_full*1/5,color="simulations")) + 
  # simulation
  geom_point(size=1) + geom_line(group=1) + 
  # data
  geom_point(data=datamart_0_5y_2012_2022,aes(x=date,y=num_positives_0_5y,color="data"),shape=21,size=1) +
  geom_line(data=datamart_0_5y_2012_2022,aes(x=date,y=num_positives_0_5y,color="data")) +
  scale_x_date(date_breaks="3 month",expand=expansion(0.01,0)) + xlab("") + ylab("DataMart vs 0.2*(simulated hospitalisations) <5y") +
  scale_color_manual(values=c("data"="red","simulations"="black")) + labs(color="") + theme_bw() + standard_theme

# calculate likelihood
# l_neg_llh = sapply((1:20)/40, function(x) 
#   left_join(simul_hosp_rate %>% filter(date<as.Date("2020-04-01") & date>as.Date("2012-07-01")) %>% group_by(date) %>% 
#             summarise(simul_hosp_sum_full=sum(hosp)), datamart_0_5y_2012_2022, by="date") %>%
#   mutate(log_lklh_poiss=dpois(x=num_positives_0_5y,lambda=x*simul_hosp_sum_full,log=T)) %>%
#   summarise(neg_LLH=-sum(log_lklh_poiss,na.rm=T)))
# 
# plot((1:20)/40,unlist(l_neg_llh),type="b",xlab="reporting rate",ylab = "negative LogLklhd")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# attack rates from 2015 Kenya study (https://academic.oup.com/jid/article/212/11/1711/2911898)
kenya_attack_rates = read_csv(file="data/kenya_attack_rates.csv") %>% 
  mutate(attack_rate=RSV_posit/n_test,sympt_attack_rate=RSV_sympt_posit/n_test)
# data.frame(age_group=c("0-1","1-5","5-15","15-40","40-100"),
  #            n_test=c(55,82,165,147,44),RSV_posit=c(33,52,73,38,9),RSV_sympt_posit=c(30,43,35,9,2)) %>%
  #            mutate(attack_rate=RSV_posit/n_test,sympt_attack_rate=RSV_sympt_posit/n_test)
# write_csv(kenya_attack_rates,file="data/kenya_attack_rates.csv")

# attack rates from Kenya aligned with model age groups
estim_attack_rates_kenya <- kenya_attack_rates[c(rep(1,2),rep(2,5),3,4,rep(5,2)),] %>% 
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name,levels=rsv_age_groups$agegroup_name)) %>% 
  rename(original_agegroups=age_group) %>% relocate(c(agegroup_name),.before=original_agegroups)

# load attack rate estimates from Hodgson
# estim_attack_rates <- read_csv("repo_data/estim_attack_rates.csv") %>% 
#   mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name)))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Latin Hypercube Sampling
# library(lhs)
# LHS -> simulate ODEs -> get incidence -> get hospitalisations -> calculate likelihood
A_opt <- randomLHS(n=10e3,k=6) # optimumLHS(n=100,k=5,maxSweeps=4,eps=0.01) # randomLHS()
# exp_dep=seq(1/4,2,1/8), age_dep=seq(1/8,1,1/16)
# convert to relevant ranges, columns: 
# 1) expdep 2) agedep 3) peak_width 4) peak_height 5) waning 6) peak week
A_opt[,1] <- qunif(A_opt[,1],min=1/3,max=1.5); A_opt[,2]<-qunif(A_opt[,2],min=0.1,max=1/2)
A_opt[,3] <- qgamma(A_opt[,3],shape=10,rate=2) # qnorm(A_opt[,3],mean=5,sd=3); 
A_opt[,4] <- qunif(A_opt[,4],min=0.3,max=1.4)
A_opt[,5] <- qnorm(A_opt[,5],mean=350,sd=50); A_opt[,5] <- 1/A_opt[,5]; 
A_opt[,6] <- round(qunif(A_opt[,6],min=43,max=48))
if (any(A_opt<0)) {A_opt=abs(A_opt); message("negative values")}
# plot distribs
# data.frame(A_opt) %>% mutate(n=row_number()) %>% pivot_longer(!n) %>%
#   ggplot() + geom_histogram(aes(x=value)) + facet_wrap(~name,scales = "free") + theme_bw()

# LOOP
# cl <- parallel::makeForkCluster(nnodes = 6); doParallel::registerDoParallel(cl)
list_AR=list_hosp_with_data=list(); k_valid=0
for (k_n in 1:nrow(A_opt)) {
# foreach(k_n = 1:nrow(A_opt), .combine='c') %dopar%  {
  exp_dep <- A_opt[k_n,1]; age_dep <- A_opt[k_n,2]
  const_delta <- 1/2; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # partable$const_delta[k_n]
  # matrix of susceptibility parameters
  delta_susc <- round(sapply(1:n_age, function(x) { delta_primary/(exp(age_dep*x)) }),6)
  # R0_calc_SIRS(C_m,delta_susc,rho,n_inf)
  if (R0_calc_SIRS(C_m,delta_susc,rho,n_inf)*(1+A_opt[k_n,4])>1) {
    tm <- proc.time(); 
    
    k_valid=k_valid+1
  message(paste0("#",k_n,", R0=",round(R0_calc_SIRS(C_m,delta_susc,rho,n_inf)*(1+A_opt[k_n,4]),2)))
  message(paste0(k_valid,"th accepted par. sets"))
  # forcing: PEAK HEIGHT + WIDTH
  # peak_week_val=round(runif(1,min=44,max=48))
  g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,
                                years_pre_post_npi=c(20,0),season_width_wks=A_opt[k_n,3],init_mt_day="07-01",
                                peak_week=A_opt[k_n,6],forcing_above_baseline=A_opt[k_n,4],npireduc_strength=0)
  approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
  approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))
  
  # WANING
  omega=A_opt[k_n,5] # 1/rho=rweibull(1, shape=4.1,scale=8.3) # rho=1/7
  # KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
  K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,
                            agegroup_durations=rsv_age_groups$duration)
  
  # all input parameters
  params <- list(list(birth_rates,deaths_vector,rsv_age_groups$stationary_popul),K_m,
                 contmatr_rowvector,inf_vars_inds, susc_vars_inds,delta_susc,mat_imm_inds)
  # SOLVE ODE SYSTEM
  tm <- proc.time(); ode_sol <- lsoda(initvals_sirs_model,timesteps,
                                      func=sirs_seasonal_forc_mat_immun,parms=params); round(proc.time()-tm,2)
  
  # Extract infection *incidence*
  df_cases_infs <- fcn_process_odesol_incid(ode_sol[(nrow(ode_sol)-1000):nrow(ode_sol),],
                                            n_age,n_inf,n_compartment,simul_start_end)

  # 1) expdep 2) agedep 3) peak_width 4) peak_height 5) waning
  # attack rate in a 26-week period (compare to Kenya): all infs from week 40 to 13
  attack_rates <- left_join(df_cases_infs %>% filter((date>=as.Date("2019-09-30") & date<as.Date("2020-03-29")) | 
                                                       (date>=as.Date("2018-10-01") & date<as.Date("2019-03-31")) ) %>%
            # for the 1st age group we count only 1st infections for the attack rate 
            # (bc in empirical studies they also count <at least one> infection, baby reinfected in a season would not be counted twice)
                  filter(!(agegroup==1 & infection>1)) %>%
                  group_by(t,date,agegroup) %>% summarise(value=sum(value)) %>% 
                  group_by(agegroup) %>% summarise(cumul_inf=sum(value)/2), # average of 2 years!
                            rsv_age_groups %>% select(agegroup_name,stationary_popul,value) %>% 
                                rename(real_pop=value) %>%
                                mutate(agegroup=as.numeric(factor(agegroup_name,levels=unique(agegroup_name)))),
                  by="agegroup") %>%
                  mutate(attack_rate_stat_pol=cumul_inf/stationary_popul, 
                         attack_rate_real_pol=cumul_inf/real_pop,k_par=k_n,
                    expdep=A_opt[k_n,1],agedep=A_opt[k_n,2],R0=R0_calc_SIRS(C_m,delta_susc,rho,n_inf),
                    peakweek=A_opt[k_n,6],peakwidth=A_opt[k_n,3],peakheight=A_opt[k_n,4],waning=A_opt[k_n,5])
  # likelihood wrt attack rates
  attack_rates_aligned_kenya = c( 
    sum(attack_rates$attack_rate_stat_pol[1:2]*attack_rates$stationary_popul[1:2]/
          sum(attack_rates$stationary_popul[1:2])), # agegroup 1-2
    sum(attack_rates$attack_rate_stat_pol[3:7]*attack_rates$stationary_popul[3:7]/
          sum(attack_rates$stationary_popul[3:7])),
    attack_rates$attack_rate_stat_pol[8],
    attack_rates$attack_rate_stat_pol[9],
    sum(attack_rates$attack_rate_stat_pol[10:11]*attack_rates$stationary_popul[10:11]/
          sum(attack_rates$stationary_popul[10:11]))
    )
  AR_log_LLH = dbinom(x=kenya_attack_rates$RSV_sympt_posit,size=kenya_attack_rates$n_test,
                      prob=attack_rates_aligned_kenya,log=T)
  # append as column: we can count LogLLH from attack rates 2x because it's for 2 seasons
  attack_rates <- attack_rates %>% mutate(AR_dbinom_neg_LLH=-sum(AR_log_LLH))
  # HOSPITALISATIONS
  # list_hosp_with_data[[k_valid]]
  hosp_daily=left_join(
    df_cases_infs %>% 
        filter( (agegroup<=7 | agegroup==11) & date>as.Date("2018-09-15")) %>%
        # mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
        #   unique(rsv_age_groups$agegroup_name))) %>%
        group_by(t,date,agegroup) %>% summarise(value=sum(value)),
    # data on hospitalisation probs
    hosp_probabilities %>% select(agegroup,prob_hosp_per_infection),by=c("agegroup")) %>%
      mutate(hosp=value*prob_hosp_per_infection) %>% select(!c(prob_hosp_per_infection)) %>%
      mutate(broad_age=ifelse(agegroup<=7,"<5y",">65y")) %>%
      group_by(date,broad_age) %>% summarise(value=sum(value),hosp=sum(hosp)) %>%
      mutate(k_par=k_n,R0=R0_calc_SIRS(C_m,delta_susc,rho,n_inf))
  # aggreg to weeks
  list_hosp_with_data_weekly[[k_valid]] = hosp_daily %>% 
    mutate(year_week=paste0(isoyear(date),"-", isoweek(date)),
           # ifelse(nchar(isoweek(date))>1,isoweek(date),paste0("0",isoweek(date))) # 
           year_week=factor(year_week,levels=unique(year_week))) %>% group_by(year_week,k_par,broad_age) %>% 
    summarise(simul_cases=sum(value),simul_hosp_sum_full=sum(hosp),date=min(date) ) %>%
    mutate(simul_hosp_rate_100k=ifelse(grepl("65y",broad_age),
                                       1e5*simul_hosp_sum_full/rsv_age_groups$stationary_popul[11],
                                       1e5*simul_hosp_sum_full/sum(rsv_age_groups$stationary_popul[1:7]) ) )  
  
  # collect outputs
  list_AR[[k_valid]] = attack_rates
  # round(proc.time()-tm,2)
  }
}
# parallel::stopCluster(cl)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# aggregate everything into dataframes

# aggregate results
all_hosp <- left_join(
    left_join(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts,
              SARIwatch_RSVhosp_over65_2018_2020_weekly_counts,by=c("wk_n","date","year")), 
    bind_rows(list_hosp_with_data_weekly), by="date") %>%
    mutate(simul_hosp_scaled=ifelse(grepl("65y",broad_age),
                            simul_hosp_sum_full*under_report_factor_over65y*(
                              pop_AGE65PLUS/rsv_age_groups$stationary_popul[11]),
                            simul_hosp_sum_full*under_report_factor_under5*(
                              pop_AGEUNDER5/sum(rsv_age_groups$stationary_popul[1:7])) ),
           log_lklh_poiss=dpois(x=ifelse(grepl("65y",broad_age),cases65plustotal,casesunder5total),
                                lambda=simul_hosp_scaled,log=T)) %>% relocate(date,.after=wk_n) %>%
    group_by(k_par,broad_age) %>% mutate(sum_neg_llh=-sum(log_lklh_poiss,na.rm=T)) # %>% distinct()
  
# all_attack_rates=bind_rows(list_AR)
all_attack_rates <- left_join(bind_rows(list_AR), 
                                all_hosp %>% select(c(k_par,sum_neg_llh,broad_age)) %>% distinct() %>% 
                                  mutate(broad_age=gsub("<|>","",broad_age)) %>% 
                          pivot_wider(id_cols=k_par,names_from=broad_age,values_from=sum_neg_llh,names_prefix="sum_neg_llh_") ) %>% 
    mutate(agegroup_name=factor(agegroup_name,levels=unique(agegroup_name)), 
           ALL_neg_log_LLH=AR_dbinom_neg_LLH+sum_neg_llh_5y+sum_neg_llh_65y)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# show all hosps
all_hosp %>% group_by(broad_age) %>% mutate(median_negLLH=median(sum_neg_llh)) %>% filter(sum_neg_llh<median_negLLH/3) %>%
  ggplot() + facet_grid(broad_age~year,scales="free") + 
  # simul
  geom_line(aes(x=date,y=simul_hosp_scaled,group=k_par,color=sum_neg_llh),size=1/2,alpha=1/3) + 
  # data
  geom_line(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                           SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)),
            aes(x=date,y=casenumber),group=1,color="red",size=1/2,linetype="dashed")+
  geom_point(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                            SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)),
             aes(x=date,y=casenumber),color="red",size=2) + 
  scale_color_gradient(low="black",high="grey") + scale_x_date(date_breaks="1 month",expand=expansion(0.01,0)) + 
  ylab("hospitalisation count") + theme_bw() + standard_theme

# expdep-agedep scatter plot
all_attack_rates %>% select(k_par,ALL_neg_log_LLH,agedep,expdep) %>% distinct() %>%
ggplot() + 
  geom_point(aes(x=agedep,y=expdep,fill=ALL_neg_log_LLH,size=log10(ALL_neg_log_LLH)),color="black",stroke=1/3,shape=21,alpha=2/3) +
  scale_fill_gradient2(low="blue",high="red",midpoint=4e3) + scale_size(trans='reverse',range=c(0.2,5)) + # breaks=(1/2+(15:0))*1e3,
  xlab("age dependence") + ylab("exposure dependence") + theme_bw() + standard_theme
# save
ggsave("extras/agedep_expdep_logLLH_scatterplot.png",width=28,height=20,units="cm")

# plot attack rates & negative Log-Lklhd for ALL components of LogLLH
ggplot(all_attack_rates %>% filter(attack_rate_stat_pol>=1/100)) + 
  facet_wrap(~agegroup_name,scales="free_y") +
  geom_point(aes(x=ALL_neg_log_LLH,y=attack_rate_stat_pol*100,color=agedep),alpha=0.6,size=1) + # shape=21, # scale(expdep)-scale(agedep)
  # geom_point(aes(x=ALL_neg_log_LLH,y=attack_rate_stat_pol*100,fill=expdep),shape=21,alpha=1/2,size=1) + #  # scale(expdep)-scale(agedep)
  geom_hline(data=estim_attack_rates_kenya,aes(yintercept=100*sympt_attack_rate),size=1/2) +
  geom_hline(data=estim_attack_rates_kenya,aes(yintercept=100*sympt_attack_rate*1/2),size=1/2,linetype="dashed") +
  geom_hline(data=estim_attack_rates_kenya,aes(yintercept=100*sympt_attack_rate*2),size=1/2,linetype="dashed") + # ,color="darkgrey"
  # geom_rect(data=estim_attack_rates_kenya,xmin=-Inf,xmax=Inf,
  #           aes(ymin=100*sympt_attack_rate/2,ymax=100*sympt_attack_rate*2),alpha=1/5,fill="grey") +
  scale_x_log10(expand=expansion(0.02,0)) + scale_y_log10() + # scale_y_continuous(breaks=(0:9)*10) + xlab("age group") + 
  scale_color_gradient2(low="blue",high="red",mid="grey",midpoint=median(all_attack_rates$agedep)) + # scale(expdep)-scale(agedep)
  # scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(all_attack_rates$expdep)) + # scale(expdep)-scale(agedep)
  ylab("% attack rate") + # scale_color_gradient(low="blue",high="red") + 
  theme_bw() + standard_theme + # theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) + 
  labs(color="agedep") + xlab("Negative log-likelihood")
# save
# ggsave(paste0("extras/LHS_attack_rate_loglkls_5e3samples.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/LHS_attack_rate_loglkls_5e3samples_expdep.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/LHS_attack_rate_loglkls_5e3samples_agedep.png"),width=28,height=16,units="cm")

# attack rates + LogLLH = <5y + attack rate
# ggplot(all_attack_rates %>% filter(attack_rate_stat_pol>=0.3/100)) + 
#   facet_wrap(~agegroup_name,scales="free_y") +
#   geom_jitter(aes(x=sum_neg_llh_5y+sum_neg_llh_65y,y=attack_rate_stat_pol*100),alpha=1/3) +
#   geom_hline(data=estim_attack_rates_kenya,aes(yintercept=100*sympt_attack_rate),size=1/2) + 
#   geom_rect(data=estim_attack_rates_kenya,xmin=-Inf,xmax=Inf,
#             aes(ymin=100*sympt_attack_rate/2,ymax=100*sympt_attack_rate*2),alpha=1/5,fill="red") +
#   scale_x_log10(expand=expansion(0.02,0)) + scale_y_log10() + # scale_y_continuous(breaks=(0:9)*10) + xlab("age group") + 
#   scale_color_gradient(low="blue",high="red") + xlab("Negative log-likelihood") + ylab("% attack rate") +
#   theme_bw() + standard_theme + labs(color="Negative Log-Likelihood")

# plot hospitalisation COUNTS
left_join(all_hosp %>% filter(sum_neg_llh<median(unique(all_hosp$sum_neg_llh))/2), 
          all_attack_rates %>% select(c(k_par,peakweek,R0)),by="k_par") %>% 
ggplot(aes(x=wk_n,group=k_par)) + 
  geom_line(aes(y=simul_hosp_scaled,color=sum_neg_llh),alpha=1/2) + # simulations
  facet_grid(broad_age~year,scales="free") + # ,scales="free" # scale_color_manual(values=c("data"="red","simulations"="black")) + 
  # data
  geom_line(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                           SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)),
            aes(y=casenumber),group=1,color="red",size=1/2,linetype="dashed") +
  geom_point(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                            SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)),
             aes(y=casenumber),group=1,color="red") + # ,shape=21
  scale_color_gradient(low="black",high="grey") + labs(color="negative log-likelihood") +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*50) + # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  xlab("") + ylab("reported hospitalisations") + theme_bw() + standard_theme # + theme(legend.text=element_text(size=11))
# save
ggsave(paste0("extras/LHS_trajs_loglkls_blackred_2agegroups_5e3.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/LHS_trajs_R0.png"),width=28,height=16,units="cm")

# plot simul hosps of last 2 yrs of "unsettled" trajectories (to see if they've reached steady state)
# all_hosp %>% mutate(wk_n=gsub("2018-|2019-|2020-","",wk_n)) %>% filter(!wk_n %in% "W16") %>% group_by(k_par,wk_n,broad_age) %>%
#   mutate(interyr_diff=simul_hosp_rate_100k[year %in% "2018-19"]-simul_hosp_rate_100k[year %in% "2019-20"]) %>%
#   group_by(k_par) %>% mutate(sum_interyr_diff=sum(abs(interyr_diff))) %>% filter(sum_interyr_diff>2e2) %>%
#   ggplot(aes(x=wk_n,group=interaction(k_par,year,broad_age))) + 
#   geom_line(aes(y=simul_hosp_rate_100k,color=factor(year),linetype=factor(broad_age)),alpha=1/2) + # simulations
#   facet_wrap(~k_par) + theme_bw() + standard_theme

# plot params showing dissimilar patterns 2018/19 vs 19/20
# all_hosp %>% filter(k_par %in% c(35,38,100)) %>% 
#   ggplot() + geom_line(aes(x=date,y=simul_hosp_sum_full,group=k_par,color=factor(k_par))) + 
#   facet_grid(broad_age~year,scales="free") + # scale_color_gradient2(low="blue",high="red",midpoint=2.1) +
#   scale_x_date(date_breaks="1 month") + theme_bw() + standard_theme

# LogLklhoods vs parameters
# left_join(all_attack_rates,all_hosp %>% select(c(k_par,sum_neg_llh))) %>% 
#   group_by(k_par) %>% summarise(sum_neg_llh=unique(sum_neg_llh),R0=unique(R0),expdep=unique(expdep),agedep=unique(agedep)) %>%
# ggplot() + geom_point(aes(x=expdep,y=sum_neg_llh,color=agedep)) + 
#   scale_color_gradient(low="grey",high="black") + theme_bw() + standard_theme

# filtering criteria: 1) 2.5x within attack rate estimates 2) 85% w42-w9 [inclus.] 3) cumul diff between 2 last years <20%
# crit 1: attack rates
all_crit_table = left_join(left_join(
left_join(all_attack_rates, estim_attack_rates_kenya,by="agegroup_name") %>% 
  mutate(crit1=attack_rate_stat_pol>=sympt_attack_rate/2.5 & attack_rate_stat_pol<=sympt_attack_rate*2.5) %>% 
  group_by(k_par) %>% 
  summarise(sum_crit1=sum(crit1),crit1=(sum_crit1>=9)), # %>% filter(sum_crit1==11)
# crit 2: seasonal concentration (of 5-year olds)
all_hosp %>% filter(broad_age %in% "<5y") %>% mutate(week=as.numeric(gsub("2018-W|2019-W|2020-W","",wk_n))) %>% 
  group_by(k_par) %>%
  summarise(all_cumul=round(sum(simul_cases)),
            in_season_cumul=round(sum(simul_cases[week>=46 | week<=2])),
            out_season_cumul=round(all_cumul-in_season_cumul),
            in_season_share=round(in_season_cumul/all_cumul,3), crit2=in_season_share>=0.6)),
# crit 3: regularity
all_hosp %>% mutate(week=as.numeric(gsub("2018-W|2019-W|2020-W","",wk_n))) %>% 
  group_by(k_par,week,broad_age) %>% filter(!week %in% 16) %>% 
  summarise(interyr_diff=simul_cases[year %in% "2018-19"]-simul_cases[year %in% "2019-20"],
            simul_cases=sum(simul_cases)/2) %>%
  group_by(k_par) %>% summarise(sum_interyr_diff=round(sum(abs(interyr_diff))),cumul_cases=round(sum(simul_cases)),
                                relat_diff=round(sum_interyr_diff/cumul_cases,4),crit3=relat_diff<=0.2)
) %>% mutate(all_crit=crit1&crit2&crit3) %>% mutate(crit_num=crit1+crit2+crit3, crit_text=paste0(ifelse(all_crit,"retain",
            paste0("reject",ifelse(!crit1,"_1",""), ifelse(!crit2,"_2",""), ifelse(!crit3,"_3","")) )) )

# plot neg-LogLH for accepted and rejected parameter sets
# summary_stats_by_crit = left_join(all_attack_rates %>% group_by(k_par,broad_age) %>% select(c(k_par,broad_age,sum_neg_llh)) %>% 
#                 distinct() %>% summarise(sum_neg_llh=sum(sum_neg_llh)),all_crit_table) %>% group_by(all_crit) %>% 
#                 summarise(median=median(sum_neg_llh),mean=mean(sum_neg_llh),
#                           CI50_low=quantile(sum_neg_llh,0.25),CI50_high=quantile(sum_neg_llh,0.75))
#

# plot distrib of 

# plot SUM of neg LLH (all components)
left_join(all_attack_rates %>% select(c(k_par,sum_neg_llh_5y,sum_neg_llh_65y,ALL_neg_log_LLH)) %>% distinct(),
          all_crit_table) %>% 
  ggplot(aes(x=all_crit,y=ALL_neg_log_LLH,color=all_crit)) + 
  geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA) + 
  geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/2,size=1.5,shape=21) +
  scale_color_manual(values=c("black","red")) + labs(color="accepted parameterisations") +
  scale_y_log10() + xlab("accepted parameterisations") + ylab("negative log-likelihood") + theme_bw() + standard_theme 
  # geom_text(aes(label=ifelse(ALL_neg_log_LLH<3e3,k_par,"")),size=3.5,
  #  position=position_jitterdodge(seed=1,dodge.width=0.9,jitter.height=1/9))
# save
ggsave(paste0("extras/jitter_llh_allcomps_filtering_5e3.png"),width=28,height=16,units="cm")

# plot LLH from <5y and attack rates (binomial)
# left_join(all_attack_rates %>% select(c(k_par,sum_neg_llh_5y,sum_neg_llh_65y,AR_dbinom_neg_LLH,ALL_neg_log_LLH)) %>% distinct(),
#           all_crit_table) %>% mutate(sum_neg_llh_5y=sum_neg_llh_5y+AR_dbinom_neg_LLH) %>%
#   ggplot(aes(x=all_crit,y=sum_neg_llh_5y,color=all_crit)) + scale_color_manual(values=c("black","red")) +
#   geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA) + 
#   geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/2,size=3) +
#   scale_y_continuous(breaks=(0:30)*5e2) + labs(color="accepted parameterisations") +
#   xlab("accepted by filtering") + ylab("negative log-likelihood") + theme_bw() + standard_theme
# # save
# ggsave(paste0("extras/jitter_llh_filtering_under5_AR_LLH.png"),width=28,height=16,units="cm")

# plot LogLLH separately by categories
left_join(all_attack_rates %>% select(c(k_par,sum_neg_llh_5y,sum_neg_llh_65y,AR_dbinom_neg_LLH,ALL_neg_log_LLH)) %>% distinct(),
          all_crit_table %>% select(c(k_par,all_crit)) ) %>% pivot_longer(!c(k_par,all_crit)) %>%
  ggplot(aes(x=all_crit,y=value,color=all_crit)) + facet_wrap(~name,scales="free") + 
  geom_boxplot(fill=NA,size=1/2,width=1/2,outlier.colour=NA) + 
  # geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/3,size=3) +
  geom_point(position=position_jitterdodge(seed=1,dodge.width=0.9),alpha=1/2,size=1.5,shape=21) +
  scale_y_log10() + # scale_y_continuous(breaks=(0:30)*5e2) + 
  labs(color="accepted parameterisations") + scale_color_manual(values=c("black","red")) +
  xlab("accepted by filtering") + ylab("negative log-likelihood") + theme_bw() + standard_theme
# save
ggsave(paste0("extras/jitter_llh_all_sep_filtering_separate.png"),width=28,height=16,units="cm")

# look at excluded parsets with likelihood lower (=better) than worst accepted parameter sets
top_parsets_rejected = left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), 
                        all_crit_table %>% select(c(k_par,all_crit)) ) %>% # group_by(all_crit) %>% 
  mutate(worst_fit_accepted=max(ALL_neg_log_LLH[all_crit],na.rm=T)) %>% # ALL_neg_log_LLH[ALL_neg_log_LLH==max(ALL_neg_log_LLH)]
  filter(ALL_neg_log_LLH<=worst_fit_accepted & !all_crit )

sel_pars <- list( top_parsets_rejected$k_par, all_crit_table$k_par[all_crit_table$all_crit] )[[1]]

# plot time courses of NUMBER of hospitalisations of good fits, color coded by accepted/rejected
left_join(left_join(all_hosp %>% mutate(iso_week=isoweek(date)),
                    all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct() ),all_crit_table) %>% 
  filter(all_crit & (iso_week>39|iso_week<=12) & ALL_neg_log_LLH<=1e3) %>% # 
ggplot(aes(x=wk_n,group=k_par)) + # 
  facet_grid(broad_age~year,scales="free") + # ,labeller=labeller(all_crit=label_both)
  geom_line(aes(y=simul_hosp_scaled,color=ALL_neg_log_LLH),alpha=1/2) + # ,linetype=crit2 # crit1 alpha=sum_crit1>=9 ,size=crit3
  # geom_ribbon(aes(ymin=ci50_low,ymax=ci50_up,fill=factor(all_crit)),color=NA) + 
  # scale_color_manual(values=c("black","blue")) + # scale_fill_manual(values=c("black","blue")) +
  scale_color_gradient(low="black",high="grey") + 
  scale_linetype_manual(values=c("dashed","solid")) + # scale_alpha_manual(values=c(1/2,1)) + scale_size_manual(values=c(1,1/2)) +
  # data
  geom_line(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                           SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)) %>%
              mutate(iso_week=isoweek(date)) %>% filter(iso_week>39|iso_week<=12),
            aes(y=casenumber),group=1,color="red",size=1/2,linetype="dashed") +
  geom_point(data=bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total),
                            SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal)) %>%
               mutate(iso_week=isoweek(date)) %>% filter(iso_week>39|iso_week<=12),aes(y=casenumber),group=1,color="red",size=2) +
  scale_y_continuous(expand=expansion(0.01,0)) + scale_x_discrete(expand=expansion(0.01,0)) + 
  labs(color="Negative Log-LLH",linetype="seasonal conc. > 60%") + # # >80% attack rates correct
  xlab("") + ylab("reported hospitalisations") + theme_bw() + standard_theme # + theme(legend.position="top")
# SAVE
# ggsave(paste0("extras/hosp_fit_rejected_goodfits_negLLH_below2e3_5e3samples.png"),width=28,height=20,units="cm")
ggsave(paste0("extras/hosp_fit_accepted_pars_5e3samples.png"),width=28,height=20,units="cm")
ggsave(paste0("extras/hosp_fit_accepted_pars_logLLH_under1e3_5e3samples.png"),width=28,height=20,units="cm")
# ggsave(paste0("extras/hosp_fit_rejected_goodfits_ci50.png"),width=28,height=16,units="cm")
# group_by(all_crit,wk_n,broad_age,year) %>% summarise(ci50_low=quantile(simul_hosp_scaled,probs=0.25),
#                                            ci50_up=quantile(simul_hosp_scaled,probs=0.75),
#                                       bestfit_simul_hosp_scaled=simul_hosp_scaled[ALL_neg_log_LLH<2*min(ALL_neg_log_LLH)],
#                                       simul_hosp_scaled=median(simul_hosp_scaled)) %>%

# plot time courses of RATE of hospitalisations of good fits, color coded by accepted/rejected
left_join(left_join(all_hosp %>% mutate(iso_week=isoweek(date)),
                    all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct() ),all_crit_table) %>% 
  filter((iso_week>39|iso_week<=12) & ALL_neg_log_LLH<=1.25e3) %>% # all_crit & 
  ggplot(aes(x=wk_n,group=k_par)) + facet_grid(broad_age~year,scales="free") + # ,labeller=labeller(all_crit=label_both)
  geom_line(aes(y=ifelse(broad_age %in% "<5y",simul_hosp_rate_100k*under_report_factor_under5,
                         simul_hosp_rate_100k*under_report_factor_over65y),alpha=ALL_neg_log_LLH,color=all_crit)) + 
  scale_alpha_continuous(trans="reverse")+scale_color_manual(values=c("black","blue")) +#scale_color_gradient(low="black",high="grey")+
  # data
  geom_line(data=bind_rows(
    SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total,rate=rate_under5yrs),
    SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal,rate=rate_65yplus)) %>%
              mutate(iso_week=isoweek(date)) %>% filter(iso_week>39|iso_week<=12),
            aes(y=rate),group=1,color="red",size=1/2,linetype="dashed") +
  geom_point(data=bind_rows(
    SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total,rate=rate_under5yrs),
    SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal,rate=rate_65yplus)) %>%
               mutate(iso_week=isoweek(date)) %>% filter(iso_week>39|iso_week<=12),aes(y=rate),group=1,color="red",size=2) +
  scale_y_continuous(expand=expansion(0.01,0)) + scale_x_discrete(expand=expansion(0.01,0)) + 
  labs(color="accepted params.",alpha="Negative log-LLH",linetype="seasonal conc. > 60%") + # # >80% attack rates correct
  xlab("") + ylab("hospitalisations per 100K population") + theme_bw() + standard_theme + theme(legend.position="top")
# save
ggsave(paste0("extras/hosp_rates_rejected_goodfits_negLLH_below1200_5e3samples.png"),width=30,height=20,units="cm")

# param distrib of accepted params
left_join(all_attack_rates[,8:ncol(all_attack_rates)] %>% distinct(),all_crit_table %>% select(k_par,all_crit)) %>% filter(R0<5) %>%
  mutate(waning=1/waning) %>%
  pivot_longer(!c(k_par,AR_dbinom_neg_LLH,sum_neg_llh_5y,sum_neg_llh_65y,ALL_neg_log_LLH,all_crit)) %>% 
ggplot(aes(x=all_crit,y=value,color=all_crit)) + 
  geom_jitter(alpha=1/2,size=1,shape=21) + geom_boxplot(size=0.4,width=0.85,outlier.colour=NA,fill=NA,color="black") + 
  scale_color_manual(values=c("black","red")) +
  facet_wrap(~name,scales="free") + labs(color="accepted parameterisations") + theme_bw() + standard_theme
# save
ggsave(paste0("extras/param_distrib_5e3samples.png"),width=28,height=20,units="cm")

# plot distrib of accepted parsets only
left_join(all_attack_rates[,8:ncol(all_attack_rates)] %>% distinct(),all_crit_table %>% select(k_par,all_crit)) %>% 
  filter(all_crit) %>% mutate(waning=1/waning) %>%
  pivot_longer(!c(k_par,AR_dbinom_neg_LLH,sum_neg_llh_5y,sum_neg_llh_65y,ALL_neg_log_LLH,all_crit)) %>% 
  ggplot(aes(x=all_crit,y=value,color=all_crit)) + 
  geom_jitter(alpha=1/2,size=1,shape=21) + geom_boxplot(size=0.4,width=0.85,outlier.colour=NA,fill=NA,color="black") + 
  scale_color_manual(values=c("black","red")) +
  facet_wrap(~name,scales="free") + labs(color="accepted parameterisations") + theme_bw() + standard_theme

# 
table((left_join(all_attack_rates[,8:ncol(all_attack_rates)] %>% distinct(),all_crit_table %>% select(k_par,all_crit)) %>% 
  filter(all_crit))$peakweek)

# table of summ statistic
param_eval_lhs = left_join(all_attack_rates[,8:ncol(all_attack_rates)] %>% distinct(),all_crit_table %>% select(k_par,all_crit)) %>% 
  mutate(waning=1/waning) %>% pivot_longer(!c(k_par,AR_dbinom_neg_LLH,sum_neg_llh_5y,sum_neg_llh_65y,ALL_neg_log_LLH,all_crit)) %>% 
  group_by(name,all_crit) %>% summarise(minval=min(value),maxval=max(value),medianval=median(value),
            ci95_l=quantile(value,probs=0.025),ci95_up=quantile(value,probs=0.975))
write_csv(param_eval_lhs,file = "extras/param_eval_lhs.csv")

# attack rates of params color-coded by filtered out/kept
# crit 1: AR, 2: seasonal conc., 3: regularity
worst_accepted_fit_llh=max(all_attack_rates$ALL_neg_log_LLH[all_attack_rates$k_par %in% 
                                                              all_crit_table$k_par[all_crit_table$all_crit]],na.rm=T)
# plot rejected parameter sets with Neg-LogLLH LOWER than WORST accepted param
left_join(all_attack_rates, all_crit_table ) %>% filter(ALL_neg_log_LLH<=worst_accepted_fit_llh ) %>%
ggplot() + 
  geom_point(aes(x=ALL_neg_log_LLH,y=attack_rate_stat_pol*100,color=all_crit),size=1.5,alpha=1/2,shape=21) +
  facet_wrap(~agegroup_name,scales="free_y") +
  geom_hline(data=estim_attack_rates_kenya,aes(yintercept=1e2*sympt_attack_rate),size=1/2) + 
  scale_y_log10() + # scale_x_log10(expand=expansion(0.02,0)) + 
  geom_rect(data=estim_attack_rates_kenya,xmin=0,xmax=Inf, # (min(all_attack_rates$sum_neg_llh)*0.9)
        aes(ymin=1e2*sympt_attack_rate/2.5,ymax=1e2*sympt_attack_rate*2.5),alpha=1/3,fill="grey") + 
  scale_color_manual(values=c("black","blue")) +
  # geom_text(aes(x=ALL_neg_log_LLH+1e2,y=attack_rate_stat_pol*100,label=ifelse(!all_crit & !crit1,k_par,"")),size=2.5) +
  xlab("Negative log-likelihood") + ylab("% attack rate") + labs(color="accepted parameterisations") +
  theme_bw() + standard_theme + theme(legend.position="top")
# save
# ggsave(paste0("extras/LHS_attack_rate_loglkls_2agegroups_filtercolorcode.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/LHS_attack_rate_loglkls_crit1.png"),width=28,height=16,units="cm")
ggsave(paste0("extras/LHS_attack_rate_loglkls.png"),width=28,height=16,units="cm")

# plot seasonal conc / regularity as well
left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
  filter(ALL_neg_log_LLH <= worst_accepted_fit_llh ) %>%
ggplot() + 
  geom_point(aes(x=in_season_share,y=ALL_neg_log_LLH,fill=ifelse(all_crit,"accepted","rejected"),size=sum_crit1),
             alpha=1/2,shape=21,stroke=0) +
  # x if attack rates <80% correct
  geom_point(data=left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
               filter(ALL_neg_log_LLH <= worst_accepted_fit_llh & (!all_crit & !crit1)), 
             aes(x=in_season_share,y=ALL_neg_log_LLH,shape="AR <80% correct",size=sum_crit1*0.6),color="black") + # show.legend=F
  # red circle if seasons irregular
  geom_point(data=left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
               filter(ALL_neg_log_LLH <= worst_accepted_fit_llh & (!crit3)), 
             aes(x=in_season_share,y=ALL_neg_log_LLH,color="biennial/irregular seasons",size=sum_crit1*0.7),shape=21,fill=NA,stroke=2/3) +
  scale_fill_manual(values=c("accepted"="blue","rejected"="darkgrey"),
                    guide=guide_legend(nrow=2,byrow=TRUE,override.aes=list(color=NA,stroke=NA))) +
  scale_color_manual(values=c("biennial/irregular seasons"="red"),guide=guide_legend(override.aes=list(size=4,stroke=1.2))) + 
  scale_shape_manual(values=c("AR <80% correct"=4),guide=guide_legend(override.aes=list(size=3))) +
  scale_radius(breaks=1+(0:5)*2,range=c(0.5,5),guide=guide_legend(nrow=2,byrow=TRUE,override.aes=list(color="black",fill=NA,shape=21))) + 
  geom_vline(xintercept=0.6,size=1/2,linetype="dashed") + xlab("seasonal concentration") + ylab("Negative log-likelihood") + 
  labs(fill="accepted \nparameterisations",color="",shape="",size="correct attack rates") + theme_bw() + standard_theme +
  theme(legend.position="top",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.text=element_text(size=12))
#
# ggsave(paste0("extras/params_phase_plot_accepted_rej_good_fits.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/params_phase_plot_all_fits_x_seasconc_y_LogLLH.png"),width=28,height=16,units="cm")
# ggsave(paste0("extras/params_phase_plot_all_fits_x_seasconc_y_LogLLH_shape_annualseas.png"),width=30,height=20,units="cm")
ggsave(paste0("extras/params_phase_plot_goodfits_x_seasconc_y_LogLLH_shape_annualseas.png"),width=30,height=20,units="cm")

# on separate facets by # of age groups with correct AR
left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
  filter(ALL_neg_log_LLH <= worst_accepted_fit_llh ) %>% mutate(`correct attack rates`=sum_crit1) %>%
ggplot() + 
  facet_wrap(~`correct attack rates`,labeller=labeller(`correct attack rates`=label_both)) +
  geom_point(aes(x=in_season_share,y=ALL_neg_log_LLH,fill=ifelse(all_crit,"accepted","rejected")),alpha=1/2,shape=21,size=3,stroke=0) +
  # x if attack rates <80% correct
  geom_point(data=left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
    filter(ALL_neg_log_LLH <= worst_accepted_fit_llh & (!all_crit & !crit1)) %>% mutate(`correct attack rates`=sum_crit1), 
    aes(x=in_season_share,y=ALL_neg_log_LLH,shape="AR <80% correct"),color="black") + # show.legend=F
  # red circle if seasons irregular
  geom_point(data=left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
               filter(ALL_neg_log_LLH <= worst_accepted_fit_llh & (!crit3)) %>% mutate(`correct attack rates`=sum_crit1), 
             aes(x=in_season_share,y=ALL_neg_log_LLH,color="biennial/irregular seasons"),shape=21,fill=NA,size=3,stroke=2/3) +
  scale_fill_manual(values=c("accepted"="blue","rejected"="darkgrey"),
                    guide=guide_legend(nrow=2,byrow=TRUE,override.aes=list(color=NA,stroke=NA))) +
  scale_color_manual(values=c("biennial/irregular seasons"="red")) + 
  scale_shape_manual(values=c("AR <80% correct"=4),guide=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=4)) ) +
  geom_vline(xintercept=0.6,size=1/2,linetype="dashed") + xlab("seasonal concentration") + ylab("Negative log-likelihood") + 
  labs(fill="accepted \nparameterisations",color="",shape="") + theme_bw() + standard_theme + theme(legend.position="top",
    axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),legend.text=element_text(size=12),strip.text=element_text(size=12))
# faceted by agegroups where AR correct
ggsave(paste0("extras/params_phase_plot_goodfits_x_seasconc_y_LogLLH_shape_annualseas_ARfaceted.png"),width=30,height=20,units="cm")

# rejected params with AR correct>=9: are they truly biennial or simulation has not settled?
# plot time courses
irregular_seas_par = (left_join(all_crit_table, all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct()) %>% 
  filter(ALL_neg_log_LLH<=worst_accepted_fit_llh & !all_crit & crit1 & crit2))$k_par
# plot
bind_rows(list_hosp_with_data_weekly[which(unique(all_hosp$k_par) %in% irregular_seas_par)]) %>%
  filter(k_par %in% irregular_seas_par & (isoweek(date)<=12|isoweek(date)>=40) & (broad_age %in% c("<5y",">65y")[1]) ) %>%
  mutate(wk_n=gsub("202.-","",gsub("201.-","",as.character(gsub("-","-W",year_week)))),
         wk_n=factor(wk_n,levels=paste0("W",c(40:53,paste0("0",1:9),11:12))),
         year=ifelse(isoweek(date)>=40,paste0(isoyear(date),"-",isoyear(date)+1),paste0(isoyear(date)-1,"-",isoyear(date))) ) %>% 
  group_by(date,broad_age) %>% mutate(k_par_broad=floor((rank(k_par)-1)/3)+1) %>% group_by(k_par_broad,date,broad_age) %>% 
  mutate(n_par=rank(k_par)) %>% filter(year %in% c("2018-2019","2020-2021") & !is.na(wk_n)) %>%
ggplot(aes(x=wk_n,group=interaction(k_par,year),color=factor(n_par),linetype=year)) + 
           facet_wrap(~k_par_broad,scales="free") + geom_line(aes(y=simul_hosp_rate_100k)) + # *under_report_factor_over65y
  scale_y_continuous(expand=expansion(0.01,0)) + scale_x_discrete(expand=expansion(0.01,0)) + 
  xlab("") + ylab("hospitalisations per 100K population") + theme_bw() + standard_theme + theme(legend.position="top")
# save
# ggsave(paste0("extras/irregular_params_good_AR_under5.png"),width=34,height=20,units="cm")
ggsave(paste0("extras/irregular_params_good_AR_over65.png"),width=34,height=20,units="cm")

# these params all look biennial. lets calculate relative (cumulative) difference btwn 2018-19 and 2020-21
relat_interyr_diff = bind_rows(list_hosp_with_data_weekly[which(unique(all_hosp$k_par) %in% irregular_seas_par)]) %>% 
  filter(k_par %in% irregular_seas_par & (isoweek(date)<=12|isoweek(date)>=40) ) %>%
  mutate(wk_n=gsub("202.-","",gsub("201.-","",as.character(gsub("-","-W",year_week)))),
         wk_n=factor(wk_n,levels=paste0("W",c(40:53,paste0("0",1:9),11:12))),
         year=ifelse(isoweek(date)>=40,paste0(isoyear(date),"-",isoyear(date)+1),paste0(isoyear(date)-1,"-",isoyear(date))) ) %>%
         filter(year %in% c("2018-2019","2020-2021") & !is.na(wk_n)) %>%
  group_by(wk_n,broad_age,k_par) %>% summarise(rate_diff=abs(diff(simul_hosp_rate_100k)),mean_inf=mean(simul_hosp_rate_100k)) %>%
  group_by(broad_age,k_par) %>% summarise(cumul_diff=sum(rate_diff),cumul_mean_inf=sum(mean_inf),
                                          cumul_relat_diff=sum(rate_diff)/sum(mean_inf))
  
# plot x: seasonal conc, y: age groups where AR correct
left_join(all_attack_rates %>% select(c(k_par,ALL_neg_log_LLH)) %>% distinct(), all_crit_table) %>%
  filter(ALL_neg_log_LLH <= worst_accepted_fit_llh ) %>%
  ggplot() + 
  geom_jitter(aes(x=in_season_share,y=sum_crit1,size=ALL_neg_log_LLH,color=all_crit,shape=crit3),alpha=1/2,
              position=position_jitter(w=0,h=0.1)) +
  geom_vline(xintercept=0.6,size=1/2,linetype="dashed") + # geom_hline(yintercept=(3:11),size=1/2,color="grey") + 
  scale_color_manual(values=c("black","blue")) + scale_shape_manual(values=c("square","circle")) + scale_size(range=c(6,2.5)) + 
  scale_y_continuous(breaks=1:11) + xlab("seasonal concentration") + ylab("age groups with correct AR") +
  labs(color="accepted parameterisations",size="Neg LogLLH",shape="annual seasons") +
  guides(color=guide_legend(nrow=2,byrow=TRUE),size=guide_legend(nrow=2,byrow=TRUE,reverse=TRUE),
         shape=guide_legend(nrow=2,byrow=TRUE)) + theme_bw() + standard_theme + 
  theme(legend.position="top",axis.text.x=element_text(size=14),axis.text.y=element_text(size=14))
# faceted by agegroups where AR correct
ggsave(paste0("extras/params_phase_plot_goodfits_x_seasconc_y_AR.png"),width=30,height=20,units="cm")

### ### ### ### ###
# seasonal concentration of hospitalisation data
bind_rows(SARIwatch_RSVhosp_under5_2018_2020_weekly_counts %>% mutate(broad_age="<5y",casenumber=casesunder5total,rate=rate_under5yrs),
          SARIwatch_RSVhosp_over65_2018_2020_weekly_counts %>% mutate(broad_age=">65y",casenumber=cases65plustotal,rate=rate_65yplus)) %>%
  mutate(week=as.numeric(gsub("2018-W|2019-W|2020-W","",wk_n))) %>% 
  group_by(year,broad_age) %>% summarise(all_cumul=sum(rate),in_season_cumul=sum(rate[week>=46|week<=2]),
                               out_season_cumul=all_cumul-in_season_cumul,in_season_share=in_season_cumul/all_cumul)
# in DataMart case data
datamart_0_5y_2012_2022 %>% mutate(epi_year=ifelse(week_number>=40,year(date),year(date)-1)) %>% group_by(epi_year) %>%
  summarise(all_cumul=sum(num_positives_0_5y,na.rm=T),in_season_cumul=sum(num_positives_0_5y[week_number>=46|week_number<=2],na.rm=T),
            out_season_cumul=all_cumul-in_season_cumul,in_season_share=in_season_cumul/all_cumul)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# hospitalisation data: RATES (publicly available)
# this is for ALL hospitalisations in SARI watch
SARI_watch_all_hosp <- read_csv(here::here("repo_data/SARI_watch_all_hosp.csv")) %>% 
  mutate(year_week=paste0(cal_year,"-",week_number),
         date=as.Date(paste0(cal_year,ifelse(week_number<10,paste0("0",week_number),week_number),"1"), "%Y%U%u"))

# ggplot(SARI_watch_all_hosp,aes(x=date,y=rate_per_100000)) + geom_point(shape=21,size=1) + geom_line() +
#   geom_vline(data=SARI_watch_all_hosp %>% filter(week_number %in% c(41,9)),
#              aes(xintercept=date),size=1/2,color="red",linetype="dashed") + 
#   scale_x_date(date_breaks="2 month") + theme_bw() + standard_theme

# hospitalisations amongst children UNDER 5YR: these are rates per 100K population
library(ISOweek) # install.packages("ISOweek"); 
SARI_watch_under5y_hosp <- read_csv(here::here("repo_data/SARI_watch_under5y_hosp.csv")) %>%
  mutate(year_week=factor(year_week,levels=unique(year_week)), 
         date=ISOweek2date(paste0(gsub("-","-W",year_week),"-1")), # as.Date(paste0(gsub("-","-W",year_week),"1"),"%Y%U%u"),
         epi_year=ifelse(week_number>=40,isoyear(date),isoyear(date)-1),
         number_hospitalisations_under5y=round(rate_under5yrs*sum(rsv_age_groups$value[1:7])/1e5)) %>% filter(!is.na(date)) %>%
  relocate(date,.before=year_week) %>% relocate(rate_under5yrs,.after=epi_year)

# fill in missing dates
plot_SARI_watch_under5y_hosp = bind_rows(SARI_watch_under5y_hosp,
                        data.frame(date=seq.Date(from=min(SARI_watch_under5y_hosp$date,na.rm=T),
                        to=max(SARI_watch_under5y_hosp$date,na.rm=T),by="week"),year_week=NA,rate_under5yrs=NA,week_number=NA) %>%
                        filter(!date %in% SARI_watch_under5y_hosp$date)) %>%
  arrange(date) %>% relocate(rate_under5yrs,.after=last_col())
# plot
varname=c("number_hospitalisations_under5y","rate_under5yrs")[2]
ggplot(plot_SARI_watch_under5y_hosp) + 
  geom_point(aes(x=date,y=get(varname)),size=3/2) + # shape=21
  geom_line(aes(x=date,y=get(varname),group=1),color="darkgrey",linetype="dashed") +
  # scale_x_discrete(breaks=every_nth(n=4)) +
  geom_rect(data=plot_SARI_watch_under5y_hosp %>% filter(week_number %in% c(40,8)) %>% group_by(epi_year) %>%
              summarise(start=max(date),end=min(date)),aes(xmin=end,xmax=start,ymin=0,ymax=Inf),color=NA,fill="red",alpha=0.1) + 
  xlab("") + ylab("rate of hospitalisations <5y (/100K)") + theme_bw() + standard_theme + scale_x_date(date_breaks="1 month") +
  scale_y_continuous(breaks=list((0:15)*250,(0:10)*10)[[ifelse(grepl("rate",varname),2,1)]],expand=expansion(0.01,0))
# geom_vline(data=SARI_watch_under5y_hosp %>% filter(week_number %in% c(40,10)),
#               aes(xintercept=date),size=1/2,color="red",linetype="dashed") +
# save
ggsave(paste0("data/SARI_watch_hosp_under5y",ifelse(grepl("rate",varname),"_rate_per_100k",""),".png"),
       width=28,height=16,units="cm")

# sum by season
# ref [25]: 13862+12862+2436=29160
# ref [26]: 16202+7108+10251=33561
SARI_watch_under5y_hosp %>% group_by(epi_year) %>% summarise(sum(number_hospitalisations_under5y))
# ~19K, this is about 60% of estimates in the literature: this is under-ascertainment (?)

### ### ### ### ### ### ### ### ### ### ### ###
# compare SARI-Watch data to simulation
ggplot() + 
  # geom_bar(data=simul_hosp_rate,aes(x=year_week,y=hosp_rate_per_100k,fill=agegroup_name),
  #          color="black",size=1/10,stat="identity") + #
  geom_bar(data=simul_hosp_rate %>% group_by(year_week) %>% summarise(hosp_rate_per_100k=sum(hosp_rate_per_100k)),
           aes(x=year_week,y=hosp_rate_per_100k),fill=NA,color="black",stat="identity") + # 
  # geom_line(data=simul_hosp_rate %>% group_by(year_week) %>% summarise(hosp_rate_per_100k=sum(hosp_rate_per_100k)),
  #                     aes(x=year_week,y=hosp_rate_per_100k,group=1)) + 
  geom_point(data=SARI_watch_under5y_hosp,aes(x=year_week,y=rate_under5yrs),color="red") + # shape=21,
  scale_x_discrete(breaks=every_nth(n=4)) + scale_y_continuous(expand=expansion(0.01,0)) + 
  xlab("") + ylab("weekly hospitalisations per 100k (<5yr)") + theme_bw() + standard_theme


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# DATAMART RSV data

datamart_0_5y_2012_2022 <- read_csv("data/datamart_0_5y_2012_2022.csv") # %>% mutate(date=dmy(week_commencing))
# plot 
datamart_0_5y_2012_2022 %>% select(!num_tests_0_5y) %>% pivot_longer(!c(week_commencing,week_number,date)) %>%
ggplot() + 
  geom_point(aes(x=date,y=value),shape=21,size=1) + geom_line(aes(x=date,y=value),linetype="dashed",color="darkgrey") +
  facet_wrap(~name,scales="free_y",nrow=2) + scale_x_date(date_breaks="3 month",expand=expansion(0.01,0)) + 
  geom_rect(data=datamart_0_5y_2012_2022 %>% filter(week_number %in% c(40,8)) %>% 
          mutate(epi_year=ifelse(week_number>=40,year(date),year(date)-1)) %>% group_by(epi_year) %>% 
          summarise(start=max(date),end=min(date)),aes(xmin=end,xmax=start,ymin=0,ymax=Inf),color=NA,fill="red",alpha=0.1) +
  xlab("") + theme_bw() + standard_theme
# save
ggsave("data/datamart_under5_num_perc_positvs.png",width=28,height=16,units="cm")

# as bar plot showing # of tests and positives
datamart_0_5y_2012_2022 %>% mutate(num_positives_0_5y=ifelse(positivity_0_5y==0,0,num_positives_0_5y),
                                   num_neg_tests_0_5y=num_tests_0_5y-num_positives_0_5y) %>%
  select(!c(week_commencing,num_tests_0_5y,positivity_0_5y)) %>% 
  rename(`number of negative RSV tests`=num_neg_tests_0_5y,`number of RSV positives`=num_positives_0_5y) %>%
  pivot_longer(!c(week_number,date)) %>%
ggplot() + geom_bar(aes(x=date,y=value,fill=name),stat="identity",color="black",size=1/7) +
  scale_x_date(date_breaks="3 month",expand=expansion(0.01/3,0))+scale_y_continuous(expand=expansion(0.01,0),breaks=(1:8)*200)+
  scale_fill_manual(values = c("white","red")) +
  geom_rect(data=datamart_0_5y_2012_2022 %>% filter(week_number %in% c(40,8)) %>% 
          mutate(epi_year=ifelse(week_number>=40,year(date),year(date)-1)) %>% group_by(epi_year) %>% 
          summarise(start=max(date),end=min(date)),aes(xmin=end,xmax=start,ymin=0,ymax=Inf),color=NA,fill="red",alpha=0.1) +
  xlab("") + ylab("") + theme_bw() + standard_theme + 
  theme(legend.position="top",legend.title=element_blank(),legend.text=element_text(size=15))
# save
ggsave("data/datamart_under5yr_num_perc_positvs_BAR.png",width=28,height=16,units="cm")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 'Respiratory viral detections by any method'
# from https://www.gov.uk/government/publications/respiratory-infections-laboratory-reports-2020

# RSV_resp_viral_detections_by_age <- read_csv("data/Respiratory viral detections by any method UK Ages.csv")  %>%
#   select(c(RSV,Year,startweek,Age)) %>% filter(Age %in% c("Under 1","1-4 Y")) %>%
#   mutate(year_week=paste0(Year,"-W",ifelse(nchar(startweek)<2,paste0("0",startweek),startweek)))
# write_csv(RSV_resp_viral_detections_by_age,file="data/RSV_resp_viral_detections_under5y.csv")
RSV_resp_viral_detections_by_age <- read_csv(file="data/RSV_resp_viral_detections_under5y.csv") %>% 
  mutate(date=ISOweek2date(paste0(year_week,"-1")),year_week=factor(year_week,levels=unique(year_week))) %>% 
  group_by(Year,Age) %>% mutate(nr_yr_sample=row_number())
# plot
ggplot(RSV_resp_viral_detections_by_age) + geom_bar(aes(x=year_week,y=RSV,fill=Age),stat="identity",color="black",size=1/3) +
  scale_x_discrete(breaks=every_nth(n=2)) + scale_y_continuous(breaks=(0:25)*5e2,expand=expansion(0.01/2,0)) + 
  xlab("") + ylab("# positive tests") + theme_bw() + standard_theme
# save
ggsave(paste0("data/RSV_resp_viral_detections_under5yr.png"), width=28,height=16,units="cm")

# compare to datamart # year_week
merge_datamart_resplab <- left_join(datamart_0_5y_2012_2022 %>% select(!week_commencing),
          RSV_resp_viral_detections_by_age %>% ungroup() %>% select(c(date,Age,nr_yr_sample,year_week,RSV)) %>%
            group_by(year_week,date) %>% 
            summarise(resp_lab_value=sum(RSV),nr_yr_sample=unique(nr_yr_sample)), by=c("date")) %>% 
  relocate(date,.before=week_number) %>% 
  filter(date>=min(date[!is.na(nr_yr_sample)]) & date<=max(date[!is.na(nr_yr_sample)]) ) %>% 
  fill(nr_yr_sample,year_week) %>% mutate(year_sample=substr(year_week,1,4)) %>% group_by(nr_yr_sample,year_sample) %>% 
  summarise(min_date=min(date),max_date=max(date),year_week=unique(year_week),
            week_number_span=paste0(week_number[date==min(date)],"-",week_number[date==max(date)]),
            datamart_positives_0_5y=sum(num_positives_0_5y,na.rm=T),num_tests_0_5y=sum(num_tests_0_5y,na.rm=T),
            resp_lab_positives_0_5y=sum(resp_lab_value,na.rm=T)) %>% arrange(year_week)
# plot
ggplot(merge_datamart_resplab) +
  # geom_bar(aes(x=year_week,y=resp_lab_positives_0_5y,color="resp. labs"),stat="identity",fill=NA) +
  geom_segment(aes(x=as.numeric(year_week)-1/2,xend=as.numeric(year_week)+1/2,
                    y=resp_lab_positives_0_5y,yend=resp_lab_positives_0_5y,color="resp. labs"),size=1) + 
  geom_point(aes(x=year_week,y=datamart_positives_0_5y,color="DataMart"),size=2.5) +
  geom_vline(aes(xintercept=as.numeric(year_week)-1/2),color="grey",size=1/2) + 
  scale_color_manual(name="data source",values=c("resp. labs"="black","DataMart"="red")) +
  scale_x_discrete(breaks=every_nth(n=2)) + scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*500) + 
  xlab("") + ylab("DataMart vs resp. labs") + theme_bw() + standard_theme  + theme(panel.grid.major.x = element_blank())
# save
ggsave(paste0("data/datamart_resp_lab_under5yr_comparison.png"), width=28,height=16,units="cm")

# ratio of values from two data sources
merge_datamart_resplab %>% mutate(datamart_resp_lab_ratio=datamart_positives_0_5y/resp_lab_positives_0_5y) %>%
ggplot(aes(x=year_week,y=datamart_resp_lab_ratio)) + geom_point() + geom_line(group=1) + 
  scale_y_continuous(breaks=(0:15)/10) + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot how susceptib. changes
# df_plot_suscept <- ode_solution_tidy %>% mutate(compartment=as.character(compartment)) %>%
#   filter((agegroup %in% c(1,2,4)) & date >= as.Date("2017-12-01") & date <= as.Date("2019-02-01")) %>%
#   relocate(c(value,value_fract),.after=last_col()) %>% relocate(date,.after=t) %>% relocate(name,.after=agegroup)
# # plot
# ggplot(df_plot_suscept %>% filter(!compartment %in% "R"), aes(x=date,y=value_fract*1e2,fill=infection)) + # geom_line() + 
#   geom_area(size=1/4,color="black",position="stack",alpha=1/2) + 
#   facet_grid(compartment~agegroup_name,scales="free_y") + # facet_wrap(~agegroup) + # ,scales="free_y"
#   geom_vline(data=df_plot_suscept %>% filter(!compartment %in% "R" & yday(date)==275),aes(xintercept=date),size=1/2,linetype="dashed") + 
#   scale_x_date(date_breaks="1 month",expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)) + 
#   xlab("") + ylab("% of age group") + theme_bw() + standard_theme 
# 
# # show only the season
# ggplot(df_plot_suscept %>% filter(date >= as.Date("2018-10-01") & date <= as.Date("2019-03-01") & !(compartment %in% c("I","R"))), 
#        aes(x=date,y=value_fract*100,fill=infection)) + geom_bar(size=1/5,color="black",stat="identity",alpha=1/2) + 
#   facet_grid(compartment~agegroup_name,scales="free_y") + # facet_wrap(~agegroup) + # ,scales="free_y"
#   scale_x_date(date_breaks="1 week",expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.03,0)) + 
#   xlab("") + ylab("% age group") + theme_bw() + standard_theme
