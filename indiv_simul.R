# Run individual simulations
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
rm(list=ls())
if (!any(row.names(installed.packages()) %in% "here")) {install.packages("here")}; library(here)
# load libraries
x1<-c("tidyverse","deSolve","lubridate","gtools")
# "rstudioapi","devtools","wpp2019","Rcpp","tsibble","pracma","qs","ungeviz","zoo","RcppRoll"
# currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# if (!any(grepl("ungeviz",row.names(installed.packages())))) { devtools::install_github("wilkelab/ungeviz") }
x2 <- x1 %in% row.names(installed.packages()); if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
# Load all packages    
lapply(x1,library,character.only=TRUE) # as.Date <- zoo::as.Date
# select <- dplyr::select; # row_number <- dplyr::row_number; summarise <- dplyr::summarise
source('fcns/essential_fcns.R')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# set plotting theme (remove Calibri if it throws an error)
standard_theme=theme(plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90,vjust=1/2),
  axis.text.y=element_text(size=9),axis.title=element_text(size=14), panel.grid.minor=element_blank()) 
# panel.grid=element_line(linetype="solid",colour="black",size=0.1) # , text=element_text(family="Calibri")
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
popul_struct <- fcn_cntr_fullpop(n_year="2020",country_sel)
# RSV age groups (population data from wpp2019)
rsv_age_groups <- fun_rsv_agegroups(standard_age_groups,popul_struct,
                                  rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,15, 45, 65),
                                  rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3), 9, 29, 19, 34))
# ONS popul estimates
ons_2020_midyear_estimates_uk <- read_csv(here::here("repo_data/ons_2020_midyear_estimates_uk.csv")) %>% 
  mutate(age_num=as.numeric(gsub("\\+","",age)))
low_inds <- findInterval(rsv_age_groups$age_low,ons_2020_midyear_estimates_uk$age_num)
high_inds <- findInterval(rsv_age_groups$age_low+rsv_age_groups$duration-0.1,ons_2020_midyear_estimates_uk$age_num)
rsv_age_groups$value <- unlist(lapply(1:length(low_inds), 
            function(x) sum(ons_2020_midyear_estimates_uk$value[low_inds[x]:high_inds[x]])*ifelse(rsv_age_groups$duration[x]<1,
            rsv_age_groups$duration[x],1) ))

# DEATHS (2019: 530841 deaths [England and Wales!]) # "uk_death_rate_byage_rsv_agegroups.csv" is for 1000 population!
# read_csv("data/uk_death_rate_byage_rsv_agegroups.csv")
# slightly adjusted age-specific death rates to get stationary population as close as possible to 2019 total & age struct
uk_death_rate=c(rep(1e-5,2)*3,rep(1e-6,5),rep(0,2),1e-6,1.79e-4) 
# number of age groups, reinfections and variables (S,I,R)
n_age=nrow(rsv_age_groups); varname_list=c('S','I','R'); n_compartment=length(varname_list); n_inf=3
dim_sys=n_age*n_compartment*n_inf; n_days_year=365
# BIRTH RATE into S_1_1 (Germany 2019: 778e3 births)
daily_births=2314; birth_rates=matrix(c(daily_births,rep(0,dim_sys-1)),dim_sys,1)
# we want population to be stationary (at 2019 value), so deaths = births
if (!any(grepl("death",colnames(rsv_age_groups)))){
  rsv_age_groups <- rsv_age_groups %>% mutate(deaths_per_person_per_day=uk_death_rate,
      stationary_popul=fcn_calc_stat_popul(rsv_age_groups,rsv_age_groups$duration,daily_births,uk_death_rate,output_type="")) 
}
# query variables: fun_sub2ind(1:3,11,"R",varname_list,n_age,n_inf)
# force of infection terms
# linear indices of the I & S variables
l_inf_susc<-fun_inf_susc_index_lists(n_age,n_inf,varname_list); inf_vars_inds<-l_inf_susc[[1]]; susc_vars_inds<-l_inf_susc[[2]]
# CONTACT MATRIX
# contact matrix from covidm ("home","work","school","other")
# if UK -> England's contact matrix
C_m_polymod <- readRDS(here::here("repo_data/UK_contact_matrix_sum.RDS"))
# create for our age groups
C_m_merged_nonrecipr <- fun_create_red_C_m(C_m_full=C_m_polymod,rsv_agegroups=rsv_age_groups,
                          orig_age_groups_duration=standard_age_groups$duration,
                          orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m <- fun_recipr_contmatr(C_m_merged_nonrecipr,age_group_sizes=rsv_age_groups$stationary_popul)
# bc of reinfections we need to input contact matrix repeatedly for the ODE system, 
# normalisation by population is to construct the force of infection terms
contmatr_rowvector=t(do.call(cbind, 
    lapply(1:nrow(C_m), function(x){diag(C_m[x,]) %*% matrix(1,n_age,n_inf)})))/rsv_age_groups$stationary_popul[col(C_m)]
# build kinetic matrix
# WANING (immunity) terms: R_i_j -> S_min(i+1,n_inf)_j
omega=1/350 # 1/runif(1,60,200)
# RECOVERY
rho=1/7 # 1/rho=rweibull(1, shape=4.1,scale=8.3)
# KINETIC MATRIX (aging terms need to be scaled by duration of age groups!)
K_m <- fun_K_m_sirs_multiage(dim_sys,n_age,n_inf,n_compartment,rho,omega,varname_list,
                              agegroup_durations=rsv_age_groups$duration)
# agegroup indices for maternal immunity
mat_imm_flag <- TRUE; mat_imm_inds<-list(fun_sub2ind(i_inf=1,j_age=1,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"R",c("S","I","R"),n_age,3),
                                         fun_sub2ind(i_inf=c(1,2,3),j_age=9,"S",c("S","I","R"),n_age,3))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# set variable/control parameters
# SUSCEPTIBILITY (normalised by age group sizes, for infection terms ~ delta*(I1+I2+...+In)*S_i/N_i)
# agedep_fact determines strength of age dependence, if agedep_fact>1, decreasing susceptibility with age
exp_dep <- 1.750; age_dep <- 0.3125 # partable$exp_dep[k_par]
const_delta <- 2.549890; delta_primary <- const_delta*exp(-exp_dep*(1:3)) # partable$const_delta[k_par]
# matrix of susceptibility parameters
delta_susc <- sapply(1:n_age, function(x) {delta_primary/(exp(age_dep*x))})

# calculate R0
R0_calc_SIRS(C_m,delta_susc,rho,n_inf)

### ### ### ### ### ### ### ### ### ### ### ### 
# DURATION of SIMULATION
# seasonal forcing (baseline level=1, forcing_strength=2 means 200% above baseline) 
# npi_reduc_strength: reduction from baseline 
# set seas lims from UK data: peak is weeks 49/50, on/off is 41,11
npi_dates=as.Date(c("2020-03-26","2021-05-17"))
g(n_years,timesteps,simul_start_end,forcing_vector_npi) %=% fun_shutdown_seasforc(npi_dates,
        years_pre_post_npi=c(6,3),season_width_wks=5,init_mt_day="06-01",
        peak_week=48,forcing_above_baseline=1,npireduc_strength=0.5)
# plot seasonal forcing
fcn_plot_seas_forc(simul_start_end,forcing_vector_npi,seas_lims_wks=c(7,42),npi_dates,date_resol="3 month")
# interpolation fcns for seas forcing & extern introds
approx_seas_forc <- approxfun(data.frame(t=timesteps,seas_force=forcing_vector_npi))
# how many introductions every 30 days?
approx_introd <- approxfun(data.frame(t=timesteps,as.numeric(timesteps %% 30==0)*5))

### ### ### ### ### ### ### ### ### ### ### ### 
# set INITIAL VALUES
# init_set can be "from scratch" (all susceptible) OR last state of a previous simulation
initvals_sirs_model <- fcn_set_initconds(rsv_age_groups$stationary_popul,
      init_set=c("previous","fromscratch")[1],init_cond_src=c("output","file")[1], # 
      input_from_prev_simul=ode_solution,init_seed=10,seed_vars="all",filename="")

# all parameters collected as a list (sent to fcn running ODE solver)
params <- list(list(birth_rates, matrix(unlist(lapply(uk_death_rate,function(x) rep(x,n_inf*n_compartment))))),
                  K_m,contmatr_rowvector,inf_vars_inds,susc_vars_inds,delta_susc,mat_imm_inds)

# SOLVE ODE SYSTEM
# can use 'lsoda' (faster) or 'lsodes'
tm <- proc.time(); ode_sol <- lsoda(initvals_sirs_model,timesteps,
                                func=sirs_seasonal_forc_mat_immun,parms=params); round(proc.time()-tm,2)
# PROCESS OUTPUT
### ### ### ### ### ### ### ### ### ### ### ### 
# this step only needed if want to investigate other variables or want to use it as input for next simulation (initial conds)
# reshape data: get all variables in `ode_solution_tidy` 
g(final_pop,ode_solution,ode_solution_tidy) %=% fun_process_simul_output(ode_sol,varname_list,
                              incidvar="newinf",incid_only=F,init_date=simul_start_end[1],n_age,n_inf,
                              rsv_age_groups,neg_thresh=-0.01)
### ### ### ### ### ### ### ### ### ### ### ### 
# Extract incident infections (do this for plots below!)
df_cases_infs <- fcn_process_odesol_incid(ode_sol,n_age,n_inf,n_compartment,simul_start_end) %>% mutate(value=round(value,2))

# select the weeks of RSV season (limits)
sel_weeks <- df_cases_infs %>% mutate(week=week(date),year=year(date)) %>% 
  filter(week %in% c(9,42,49)) %>% group_by(year,agegroup,week) %>% 
  filter(date==min(date) & infection==1)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT RESULTS
# disaggregate by age, summed over infection type (1/2/3)
df_cases_infs %>% 
  filter(agegroup<=9 & date>max(c(as.Date("2016-09-01"),min(df_cases_infs$date))) & date<as.Date("2023-05-01")) %>%
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
  unique(rsv_age_groups$agegroup_name))) %>% 
  group_by(t,date,agegroup_name) %>% summarise(value=sum(value)) %>%
  ggplot() + geom_line(aes(x=date,y=value)) + facet_wrap(~agegroup_name,scales="free_y") + 
  geom_vline(data=sel_weeks %>% filter(agegroup<=9),
             aes(xintercept=date,linetype=ifelse(week==unique(sel_weeks$week)[2],"solid","dashed"),
                 size=ifelse(week==49,"a","b")),show.legend=F) + 
  scale_size_manual(breaks=c("a","b"),values=c(1/4,1/8)) +
  scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + xlab("") + ylab("# cases") +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=1/100) + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT incident infections, disaggregated by age AND 1/2/3rd infection
df_cases_infs %>% 
  filter(t %% 7==0 & agegroup<=9 & date>as.Date("2018-09-01") & date<as.Date("2022-04-01")) %>% 
  mutate(agegroup_name=factor(rsv_age_groups$agegroup_name[agegroup],
                              unique(rsv_age_groups$agegroup_name))) %>%
  ggplot() + geom_area(aes(x=date,y=value,fill=factor(infection)),position=position_stack(reverse=T),
                       alpha=0.3,color="black",size=0.25) + 
  facet_wrap(~agegroup_name,scales="free_y") + xlab("") + 
  geom_vline(data=sel_weeks %>% filter(agegroup<=9),
             aes(xintercept=date,linetype=ifelse(week==unique(sel_weeks$week)[2],"solid","dashed"),
                 size=ifelse(week==49,"a","b")),show.legend=F) + 
  scale_size_manual(breaks=c("a","b"),values=c(1/4,1/8)) +
  scale_x_date(date_breaks="2 month",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) +
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=1/100) +
  theme_bw() + standard_theme + ylab("# cases") # ifelse(varname=="value","# cases","% age group")

### ### ### ### ### ### ### ### ### ### ### ### ### 
# sum across age groups
df_cases_infs %>% filter(agegroup<=9 & date>as.Date("2018-09-01") & # t %% 7==0 & 
      date<as.Date("2023-04-01")) %>% group_by(date,infection) %>% summarise(value=sum(value)) %>%
  ggplot(aes(x=date,y=value,fill=factor(infection),group=infection)) + 
  geom_area(position=position_stack(reverse=T),alpha=0.3,color="black",size=0.25) + 
  scale_x_date(date_breaks="3 months",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0)) + 
  geom_rect(xmin=npi_dates[1],xmax=npi_dates[2],ymin=-Inf,ymax=Inf,fill="pink",alpha=0.01) + 
  geom_vline(data=sel_weeks %>% filter(week %in% c(42,9)),
             aes(xintercept=date,linetype="solid"),size=1/4,show.legend=F,linetype="dashed") + 
  ylab("# new infections") + theme_bw() + standard_theme
