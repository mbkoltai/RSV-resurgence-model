rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# library(contactdata); library(fitdistrplus);  library(bbmle); library(Rcpp); library(GillespieSSA)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","plotly","Rcpp","zoo","lubridate","tsibble","pracma","qs","ungeviz"),
       library,character.only=TRUE)
source('fcns/RSV_model_functions.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SET PARAMETERS --------------------------------------------------------
### ### ### ### ### ### ### ###
# constant parameters
# selected country
country_sel="United Kingdom"
# time resolution (in days)
elem_time_step=1
# population data
standard_age_groups <- fun_cntr_agestr(country_sel,i_year="2020",seq(0,75,5),c(seq(4,74,5),99))
popul_struct=fcn_cntr_fullpop(n_year="2020",country_sel)
# RSV age groups (population data from wpp2019)
rsv_age_groups<-fun_rsv_agegroups(standard_age_groups,popul_struct,rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,15, 45, 65),
                                  rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3), 9, 29, 19, 34))
# number of age groups, reinfections and variables (S,I,R)
rsv_age_groups$value=rsv_age_groups$value*67e6/sum(rsv_age_groups$value)
# DEATHS (2019: 530841 deaths [England and Wales!]) # "uk_death_rate_byage_rsv_agegroups.csv" is for 1000 population!
# read_csv("data/uk_death_rate_byage_rsv_agegroups.csv")
# i slightly adjusted the age-specific deaths rates to get a stationary population at the 2019 total and age struct
uk_death_rate=c(rep(1e-5,2)*3,rep(1e-6,5),rep(0,2),1e-6,1.79e-4) 
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
daily_births=2314; birth_rates=matrix(c(daily_births,rep(0,dim_sys-1)),dim_sys,1)
# we want population to be stationary (at 2019 or 2020 value), so deaths = births
if (!any(grepl("death",colnames(rsv_age_groups)))){
  rsv_age_groups <- rsv_age_groups %>% mutate(deaths_per_person_per_day=uk_death_rate,
        stationary_popul=fcn_calc_stat_popul(rsv_age_groups,rsv_age_groups$duration,daily_births,uk_death_rate,output_type="") ) }
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3

# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc=fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds=l_inf_susc[[1]]; susc_vars_inds=l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
# if UK -> England's contact matrix # check: cm_parameters_SEI3R(cm_uk_locations("UK", 1))$pop[[1]]$matrices 
list_contmatrs=fun_covidm_contactmatrix(country_sel,currentdir_path,cm_path=cm_path) 
# make matrix reciprocal
C_m_polymod=Reduce('+',list_contmatrs) # fun_recipr_contmatr(Reduce('+',list_contmatrs),age_group_sizes=standard_age_groups$values)
# create for our age groups
C_m_merged_nonrecipr=fun_create_red_C_m(C_m_polymod,rsv_age_groups,
                                        orig_age_groups_duration=standard_age_groups$duration,orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m_merged_nonrecipr,age_group_sizes=rsv_age_groups$stationary_popul)
# bc of reinfections we need to input contact matrix repeatedly
contmatr_rowvector=t(do.call(cbind, lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m=fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,agegroup_durations=rsv_age_groups$duration)
# estimated attack rates
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
agedep_fact=1; delta_primary=c(0.09,0.07,0.05)/3 # c(0.09,0.07,0.05)/2.5
# c(0.27,0.03,0.01)/5 # c(0.21,0.11,0.01)/5 # c(0.09,0.07,0.05)/1.6 # rep(0.15,3) # rep(0.09,3)
delta_susc <- sapply(1:n_age, function(x) {delta_primary/((agedep_fact^(x-1))*rsv_age_groups$stationary_popul[x])})
delta_susc_prop <- delta_susc*matrix(rep(rsv_age_groups$stationary_popul,3),nrow=3,byrow=T)
dep_subfolder_name<-fun_subfld(delta_primary,delta_susc_prop)
### PLOT susceptibility: 
# fcn_plot_suscept_table(fcn_suscept_agedeptable(rsv_age_groups,delta_susc,n_inf)) 
# ggsave(paste0("simul_output/suscept_age_dep/",agedep_fact,"delta_prim",unique(delta_primary),".png"),width=32,height=22,units="cm")
# calculate R0 (at max seasonal forcing=1)
R0_calc_SIRS(C_m,delta_susc_prop,rho,n_inf)
####
# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) | npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-04-01")); seaspeakval=1/3; seasforc_width_wks=3
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,years_pre_post_npi=c(5,3),
            season_width_wks=seasforc_width_wks,init_mt_day="06-01",peak_week=44,forcing_above_baseline=seaspeakval,npireduc_strength=0.5)
# SAVE: ggsave(paste0("simul_output/NPI_y",npi_year,"_on",preseas_npi_on,"w_off",postseas_npi_off,"w.png"),units="cm",height=10,width=20)
# estimated attack rates
estim_attack_rates <- data.frame(agegroup_name=paste0("age=",rsv_age_groups$agegroup_name,"yr"),
                                 median_est=c(rep(65,4),rep(40,4),10,8,5)) %>% mutate(min_est=median_est*0.25,max_est=median_est*1.75)
